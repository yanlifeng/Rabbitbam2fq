//
// Created by 赵展 on 2021/2/22.
//

#include "Buffer.h"
#include <htslib/sam.h>
#include <chrono>
#include <fstream>
#include <thread>
#include <mutex>
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
using namespace std;
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};


ofstream fout;
int chro_number=0;
const int chro_number_all=25;
char chromosome[][50]= {
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
        "chrX", "chrY", "*"
};
int n_thread=1;
int num[50] = {0};
char *path;
char *index_path;
void write_pack(Buffer *buffer){
    while(!buffer->is_complete()){
        std::this_thread::sleep_for(chrono::milliseconds(10));
        //cout << thread_have << endl;
        buffer->output();
    }
}
/*
 * number into change is 6583191
 * 6583494
change time is :     	0.963243	 sec
 */


void consumer_pack(Buffer *buffer,int id){
    htsFile *hin;
    sam_hdr_t *hdr;
    hts_idx_t *idx = NULL;
    bam1_t *b;


    if ((hin=hts_open(path, "r"))==NULL){
        printf("Can`t open this file!\n");
    }
    if ((idx = sam_index_load(hin, index_path)) == 0) {
        fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
    }
    if ((hdr = sam_hdr_read(hin)) == NULL) {
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    hts_itr_t *iter;
    int fg=1;
    int thread_chro_num=0;
    pair<char *,int> thread_buffer;
    thread_buffer=buffer->getCap();
    while (fg){
        buffer->mtx.lock();
        if (chro_number>=chro_number_all) { fg=0;}
        else thread_chro_num=chro_number++;
        buffer->mtx.unlock();
        if (!fg) break;
        if ((iter = sam_itr_querys(idx, hdr, chromosome[thread_chro_num])) == 0) {continue;}
        while (sam_itr_next(hin, iter, b)>=0) {
            if (b->core.flag&2048) continue;
            num[id]++;
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
    hts_close(hin);
}
int main(int argc,char* argv[]){
    TDEF(fq)
    TSTART(fq)
//    n_thread = atoi(argv[1]);
//    path = argv[2];
//    index_path = argv[3];
//    fout.open("./out.fq");
    n_thread=15;
    path = "/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam";
    index_path = "/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam.bai";
    fout.open("/Users/zhaozhan/CLionProjects/BAM/out.fq");
    BufferConfig config(30,n_thread,10000000);
    Buffer buffer(&config,&fout);
    thread **Bam = new thread *[n_thread+1];

    for (int i=0;i<n_thread;i++)
        Bam[i]=new thread(&consumer_pack,&buffer,i);
    Bam[n_thread]=new thread(&write_pack,&buffer);
    for (int i=0;i<n_thread+1;i++)
        Bam[i]->join();
    long long num_all=0;
    for (int i=0;i<n_thread;i++) num_all+=num[i];
    cout << "number into change is " << num_all << endl;
    TEND(fq)
    TPRINT(fq,"change time is : ");
}
/*
 * 6512186
 */