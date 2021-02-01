#include <cstdio>
#include <cstdlib>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>


#ifdef _OPENMP

#include <omp.h>

#endif

#define MX 100000000
using namespace std;


uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int main(int argc, char **argv) {
    clock_t t0 = clock();

    bam_hdr_t *header;
    bam1_t *aln = bam_init1();

    samFile *in = sam_open(argv[1], "r");
    char *outFileName = argv[2];
    ofstream outStream;
    outStream.open(outFileName, ifstream::out);
    header = sam_hdr_read(in);

    uint8_t *seq;
    uint8_t *qul;
    char *data = new char[MX + 10000];
    int32_t lseq;
    char *qname;
    int N = 0, M = 0, pos = 0, nameLen;
#ifdef _OPENMP
    printf("openmp is open\n");
    printf("use %d threads\n", omp_get_num_procs());
#endif
    printf("cost %.5f\n", (clock() - t0) * 1e-6);
    t0 = clock();
    while (true) {
        struct BGZF *fp = in->fp.bgzf;
        int r = bam_read1(in->fp.bgzf, aln);
        if (header && r >= 0) {
            if (aln->core.tid >= header->n_targets || aln->core.tid < -1 ||
                aln->core.mtid >= header->n_targets || aln->core.mtid < -1) {
                errno = ERANGE;
                break;
            }
        }
        //        printf("block_length %d\n", fp->block_length);
//        printf("block_offset %d\n", fp->block_offset);
//        printf("res %d\n", r);
        if (r < 0)break;
        N++;
        if (aln->core.flag & 2048)continue;
        M++;
//        seq = bam_get_seq(aln);
//        qul = bam_get_qual(aln);
//        qname = bam_get_qname(aln);
//        nameLen = strlen(qname);
//        if (pos > MX) {
//            outStream.write(data, pos);
//            pos = 0;
//        }
//        data[pos++] = '@';
//        memcpy(data + pos, qname, nameLen);
//        pos += nameLen;
//        data[pos++] = '\n';
//        lseq = aln->core.l_qseq;
//        if (aln->core.flag & 16) {
//            for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
//                data[pos + j] = BaseRever[bam_seqi(seq, i)];
//            }
//            pos += lseq;
//            data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
//            for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
//                data[pos + j] = char(qul[i] + 33);
//            }
//            pos += lseq;
//        } else {
//            for (int i = 0; i < lseq; ++i) {
//                data[pos + i] = Base[bam_seqi(seq, i)];
//            }
//            pos += lseq;
//            data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
//            for (int i = 0; i < lseq; ++i) {
//                data[pos + i] = char(qul[i] + 33);
//            }
//            pos += lseq;
//        }
//        data[pos++] = '\n';
    }
    printf("cost %.5f\n", (clock() - t0) * 1e-6);
    outStream.write(data, pos);
    printf("total process %d reads\n", N);
    printf("total print %d reads\n", M);
    sam_close(in);
    return 0;
}

