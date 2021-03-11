//// Created by 赵展 on 2021/2/18.
////
#include <iostream>
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
const int n_thread=16;
const long long Maxn = 1e8;
char Buffer[n_thread][Maxn+10000];
long long pos[n_thread];
int chro_number=0;
uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};


ofstream fout;
mutex mtx;
char chromosome[][50]= {
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
        "chrX", "chrY", "*"
};
void consumer_pack(int id){
    cout << id << endl;
    htsFile *hin;
    sam_hdr_t *hdr;
    hts_idx_t *idx = NULL;
    const char *path = "/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam";
    const char *index_path = "/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam.bai";

    if ((hin=hts_open(path, "r"))==NULL){
        printf("Can`t open this file!\n");
    }
    if ((idx = sam_index_load(hin, index_path)) == 0) {
        fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
    }
    if ((hdr = sam_hdr_read(hin)) == NULL) {
    }
    pos[id]=0;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    hts_itr_t *iter;
    int fg=1;
    int thread_chro_num=0;
    while (fg){
        mtx.lock();
        if (chro_number>=25) {fg=0;}
        thread_chro_num=chro_number++;
        mtx.unlock();
        if (!fg) break;
        if ((iter = sam_itr_querys(idx, hdr, chromosome[thread_chro_num])) == 0) {continue;}
        while (sam_itr_next(hin, iter, b)>=0) {
            if (b->core.flag&2048) continue;
            seq=bam_get_seq(b);
            long long len = strlen(bam_get_qname(b));
            if (pos[id]>Maxn-1)
            {
                mtx.lock();
                fout.write(Buffer[id],pos[id]);
                mtx.unlock();
                pos[id]=0;
            }
            Buffer[id][pos[id]++]='@';
            memcpy(Buffer[id]+pos[id],bam_get_qname(b),len);
            pos[id]+=len;
            Buffer[id][pos[id]++]='\n';
            if (b->core.flag&16){
                for (int i=0;i<b->core.l_qseq;i++) Buffer[id][pos[id]+b->core.l_qseq-1-i] = BaseRever[bam_seqi(seq,i)];
                pos[id]+=b->core.l_qseq;
                Buffer[id][pos[id]++]='\n';Buffer[id][pos[id]++]='+';Buffer[id][pos[id]++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ Buffer[id][pos[id]+i]=(char)quality[b->core.l_qseq-i-1]+33;}
                pos[id]+=b->core.l_qseq;
            }else{
                for (int i=0;i<b->core.l_qseq;i++) Buffer[id][pos[id]+i] = Base[bam_seqi(seq,i)];
                pos[id]+=b->core.l_qseq;
                Buffer[id][pos[id]++]='\n';Buffer[id][pos[id]++]='+';Buffer[id][pos[id]++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ Buffer[id][pos[id]+i]=(char)quality[i]+33;}
                pos[id]+=b->core.l_qseq;
            }
            Buffer[id][pos[id]++]='\n';
        }
        mtx.lock();
        fout.write(Buffer[id],pos[id]);
        mtx.unlock();
        pos[id]=0;
    }
    hts_close(hin);
}
int main(){
    TDEF(fq)
    TSTART(fq)
    fout.open("/Users/zhaozhan/CLionProjects/BAM/my1.fq");
    thread **Bam = new thread *[n_thread];
    for (int i=0;i<n_thread;i++)
        Bam[i]=new thread(&consumer_pack,i);
    for (int i=0;i<n_thread;i++)
        Bam[i]->join();
    TEND(fq)
    TPRINT(fq,"change time is : ");
}
/*
 * 1. 不确定 tid 到底有多少个，现在采用了25个tid ，但是输出大小为1.5GB，标准答案应该是1.52GB
 *    线程数   |   时间  5.5
 *      2     ｜  2.8-3.0
 *      4     ｜  1.6-1.75
 *      8     ｜  1.15-1.2
 *     16     |   1.0
 */