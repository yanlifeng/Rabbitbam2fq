//
// Created by 赵展 on 2021/1/30.
//

#include <cstdio>
#include <cstdlib>
#include <htslib/sam.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <chrono>
#include <fstream>

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
#ifdef _OPENMP
#include <omp.h>
#endif

#define MX 100000000
using namespace std;


uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int main(int argc, char **argv) {
    bam_hdr_t *header;
    bam1_t *aln = bam_init1();

    samFile *in = sam_open("/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam", "r");
    char *outFileName = "/Users/zhaozhan/CLionProjects/BAM/ylf.txt";
    ofstream outStream;
    outStream.open(outFileName, ifstream::out);
    header = sam_hdr_read(in);

    uint8_t *seq;
    uint8_t *qul;
    char *data = new char[MX + 10000];
    int32_t lseq;
    char *qname;
    int N = 0, M = 0, pos, nameLen;
#ifdef _OPENMP
    printf("openmp is open\n");
    printf("use %d threads\n", omp_get_num_procs());
#endif
    TDEF(change)
    TSTART(change)
    while (sam_read1(in, header, aln) >= 0) {
        N++;
        if (aln->core.flag & 2048)continue;
        M++;
        seq = bam_get_seq(aln);
        qul = bam_get_qual(aln);
        qname = bam_get_qname(aln);
        nameLen = strlen(qname);
        if (pos > MX) {
            outStream.write(data, pos);
            pos = 0;
        }
        data[pos++] = '@';
        memcpy(data + pos, qname, nameLen);
        pos += nameLen;
        data[pos++] = '\n';
        lseq = aln->core.l_qseq;
        if (aln->core.flag & 16) {
            for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                data[pos + j] = BaseRever[bam_seqi(seq, i)];
            }
            pos += lseq;
            data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
            for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                data[pos + j] = char(qul[i] + 33);
            }
            pos += lseq;
        } else {
            for (int i = 0; i < lseq; ++i) {
                data[pos + i] = Base[bam_seqi(seq, i)];
            }
            pos += lseq;
            data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
            for (int i = 0; i < lseq; ++i) {
                data[pos + i] = char(qul[i] + 33);
            }
            pos += lseq;
        }
        data[pos++] = '\n';
    }
    outStream.write(data, pos);
    printf("total process %d reads\n", N);
    printf("total print %d reads\n", M);
    sam_close(in);
    TEND(change)
    TPRINT(change,"time is ")
    return 0;
}