//
// Created by 赵展 on 2021/3/10.
//
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

typedef std::chrono::high_resolution_clock Clock;
#ifdef CYCLING
#define TDEF(x_) static unsigned long long int x_##_t0, x_##_t1;
    #define TSTART(x_) x_##_t0 = __rdtsc();
    #define TEND(x_) x_##_t1 = __rdtsc();
    #define TPRINT(x_, str) printf("%-20s \t%.6f\t M cycles\n", str, (double)(x_##_t1 - x_##_t0)/1e6);
#elif defined TIMING
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
#else
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
#endif
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int n_thread=1;
long long NUM_N[500]={0};
long long NUM_M[500]={0};

void read_pack(BGZF *fp,BamBlock *block){
    pair<bam_block *,int> b;
    b=block->getEmpty();
    int count=0;
    while(read_block(fp,b.first)==0){
        block->inputblock(b.second);
        //printf("read block is %d\n",++count);
        b=block->getEmpty();
    }
    block->ReadComplete();
}

void write_pack(Buffer *buffer){
    while(!buffer->is_complete()){
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}

void consumer_pack(BamBlock *block,Buffer *buffer,int id){

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    pair<bam_block *,int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *)malloc(sizeof(bam_block));
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    pair<char *,int> thread_buffer;
    thread_buffer=buffer->getCap();
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp,b,0)>=0) {
            NUM_N[id]++;
            if (b->core.flag&2048) continue;
            NUM_M[id]++;
            seq=bam_get_seq(b);
            long long len = strlen(bam_get_qname(b));
            if (pos>buffer->config->Maxn)
            {
                buffer->initoutput(thread_buffer.second,pos);
                thread_buffer=buffer->getCap();
                pos=0;
            }
            thread_buffer.first[pos++]='@';
            memcpy(thread_buffer.first+pos,bam_get_qname(b),len);
            pos+=len;
            thread_buffer.first[pos++]='\n';
            if (b->core.flag&16){
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+b->core.l_qseq-1-i] = BaseRever[bam_seqi(seq,i)];
                pos+=b->core.l_qseq;
                thread_buffer.first[pos++]='\n';thread_buffer.first[pos++]='+';thread_buffer.first[pos++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ thread_buffer.first[pos+i]=(char)quality[b->core.l_qseq-i-1]+33;}
                pos+=b->core.l_qseq;
            }else{
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+i] = Base[bam_seqi(seq,i)];
                pos+=b->core.l_qseq;
                thread_buffer.first[pos++]='\n';thread_buffer.first[pos++]='+';thread_buffer.first[pos++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ thread_buffer.first[pos+i]=(char)quality[i]+33;}
                pos+=b->core.l_qseq;
            }
            thread_buffer.first[pos++]='\n';
        }
    }
    buffer->initoutput(thread_buffer.second,pos);
    buffer->complete_thread();
}

int main(int argc,char* argv[]){
    TDEF(fq)
    TSTART(fq)
    printf("Starting Running\n");
    n_thread=atoi(argv[1]);
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    ofstream fout;
    fout.open("./bam.fq");
    // /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
    // ../data/NC_T_1.sorted.bam
    if ((sin=sam_open("../data/NC_N_1.sorted.bam", "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    BufferConfig config(150,n_thread,10000000);
    Buffer buffer(&config,&fout);
    BamBlockConfig bamconfig(20000);
    BamBlock block(&bamconfig);
    thread **Bam = new thread *[n_thread+2];
    Bam[0]=new thread(&read_pack,sin->fp.bgzf,&block);
    for (int i=1;i<=n_thread;i++)
        Bam[i]=new thread(&consumer_pack,&block,&buffer,i);
    Bam[n_thread+1]=new thread(&write_pack,&buffer);
    for (int i=0;i<n_thread+2;i++)
        Bam[i]->join();
    long long N=0,M=0;
    for (int i=0;i<=n_thread;i++) N+=NUM_N[i];
    for (int i=0;i<=n_thread;i++) M+=NUM_M[i];
    printf("total read is %lld\n",N);
    printf("totol process is %lld\n",M);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}

/*
 * in 金牌
 *  thread_number |    time    |   lib
 *         1      |   1674     |   zlib
 *         2      |    874     |   zlib
 *         4      |    482     |   zlib
 *         8      |    279     |   zlib
 *        16      |    166     |   zlib
 *        24      |    132     |   zlib
 *        32      |    111     |   zlib
 *
 *  thread_number |    time    |   lib
 *         1      |    865     |   libdeflate
 *         2      |    494     |   libdeflate
 *         4      |    304     |   libdeflate
 *         8      |    190     |   libdeflate
 *        16      |    138     |   libdeflate
 *        24      |    111     |   libdeflate
 *        30      |    104     |   libdeflate
 *
 *
 *
 * in AMD
 *
 *  thread_number |    time    |   lib
 *         1      |    876     |   libdeflate
 *         2      |    508     |   libdeflate
 *         4      |    298     |   libdeflate
 *         8      |    200     |   libdeflate
 *        16      |    115     |   libdeflate
 *        20      |     91     |   libdeflate
 *        24      |     75     |   libdeflate
 *        28      |     63     |   libdeflate
 *        32      |     53     |   libdeflate
 *        40      |     40     |   libdeflate
 *        94      |     28     |   libdeflate
 */

