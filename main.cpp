#include <cstdio>
#include <cstdlib>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>

#include "mytimer.hpp"

#ifdef _OPENMP

#include <omp.h>

#endif

#define MX 100000000
using namespace std;


uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int main(int argc, char **argv) {
    double t0 = get_wall_time();

    bam_hdr_t *header;
    bam1_t *aln = bam_init1();
    char *inFileName = argv[1];
    char *outFileName = argv[2];
    printf("in %s\n", inFileName);
    printf("out %s\n", outFileName);
    samFile *in = sam_open(inFileName, "r");
    ofstream outStream;
    outStream.open(outFileName, ifstream::out);
    header = sam_hdr_read(in);
    uint8_t *seq;
    uint8_t *qul;
    char *data = new char[MX + 10000];
    int32_t lseq;
    char *qname;
    int N = 0, M = 0;
    int cnt = 0;
    int nameLen, pos = 0;
    printf("cost %.5f\n", get_wall_time() - t0);
    t0 = get_wall_time();
    double c1 = 0, c2 = 0, c3 = 0;
    while (true) {
//        double tt = get_wall_time();
        int r = bam_read1(in->fp.bgzf, aln);

//        if (header && r >= 0) {
//            if (aln->core.tid >= header->n_targets || aln->core.tid < -1 ||
//                aln->core.mtid >= header->n_targets || aln->core.mtid < -1) {
//get_wall_time()                errno = ERANGE;
//                break;
//            }
//        }
        if (r < 0)break;
        N++;
        if (aln->core.flag & 2048)continue;
        M++;
//        c1 += get_wall_time() - tt;
//        tt = get_wall_time();
        qul = bam_get_qual(aln);
        qname = bam_get_qname(aln);
        seq = bam_get_seq(aln);
        nameLen = aln->core.l_qname - aln->core.l_extranul - 1;

//        c2 += get_wall_time() - tt;
//        tt = get_wall_time();

        if (pos > MX) {
            cnt++;
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
//        c3 += get_wall_time() - tt;
    }
    outStream.write(data, pos);
    printf("cost %.5f\n", get_wall_time() - t0);
//    printf("cost1 %.5f\n", c1);
//    printf("cost2 %.5f\n", c2);
//    printf("cost3 %.5f\n", c3);
    printf("total process %d reads\n", N);
    printf("total print %d reads\n", M);
    printf("total write count %d\n", cnt);
    sam_close(in);
    return 0;
}

