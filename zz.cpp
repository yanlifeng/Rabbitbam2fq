//
// Created by 赵展 on 2021/1/29.
//
//
// Created by 赵展 on 2021/1/27.
//
#include <iostream>
#include <htslib/sam.h>
#include <chrono>
#include <fstream>
#include <cstring>
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
int8_t seq_comp_table[16] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};

const long long Maxn = 1e12;
char Buffer[Maxn];
long long pos = 0;


int main() {
    // 时间 ： (fprintf) 7.35s  (stream) 7.03s
    const char *path = "/Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam";
    htsFile *in;
    sam_hdr_t *hdr;
    bam1_t *b;
    char moder[8];
    strcpy(moder, "r");
    if ((in = hts_open(path, moder)) == NULL) {
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(in)) == NULL) {
        return 0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    int num = 1, ret = 0;
//    FILE *file = fopen("/Users/zhaozhan/CLionProjects/BAM/my1.txt", "w");
    ofstream fout;
    fout.open("/Users/ylf9811/Desktop/QCQCQC/data/hg19/zz.fq");
    TDEF(fq)
    TSTART(fq)
    while ((ret = sam_read1(in, hdr, b)) >= 0 && num++) {

        char *seq = bam_get_qname(b);
        long long len = strlen(bam_get_qname(b));
        if (pos + len + 2 * b->core.l_qseq + 4 < Maxn - 1) {
            memcpy(Buffer + pos, bam_get_qname(b), len);
            pos += len;
            Buffer[pos] = '\n';
            pos++;
            if (b->core.flag & 16) {
                for (int i = 0; i < b->core.l_qseq; i++)
                    Buffer[pos + b->core.l_qseq - 1 - i] = seq_nt16_str[seq_comp_table[bam_seqi(seq, i)]];
                pos += b->core.l_qseq;
                Buffer[pos] = '\n';
                pos++;
                memcpy(Buffer + pos, (char *) bam_get_qual(b), b->core.l_qseq);
                for (int i = 0; i < (b->core.l_qseq + 1) >> 1; i++) {
                    swap(Buffer[pos + b->core.l_qseq - 1 - i], Buffer[pos + i]);
                    Buffer[pos + b->core.l_qseq - 1 - i] += 33;
                    if (b->core.l_qseq != 2 * i + 1) Buffer[pos + i] += 33;
                }
                pos += b->core.l_qseq;
                Buffer[pos] = '\n';
                pos++;
            } else {
                for (int i = 0; i < b->core.l_qseq; i++) Buffer[pos + i] = seq_nt16_str[bam_seqi(seq, i)];
                pos += b->core.l_qseq;
                Buffer[pos] = '\n';
                pos++;
                memcpy(Buffer + pos, (char *) bam_get_qual(b), b->core.l_qseq);
                for (int i = 0; i < b->core.l_qseq; i++) { Buffer[pos + i] = Buffer[pos + i] + 33; }
                pos += b->core.l_qseq;
                Buffer[pos] = '\n';
                pos++;
            }
        } else {
//            Buffer[pos]='\0';
//            fprintf(file,Buffer);
            fout.write(Buffer, pos);
            pos = 0;
        }
    }
    if (pos != 0) {
//        Buffer[pos]='\0';
//        fprintf(file,Buffer);
        fout.write(Buffer, pos);
        pos = 0;
    }
    hts_close(in);
    TEND(fq)
    TPRINT(fq, "change time is : ");
}