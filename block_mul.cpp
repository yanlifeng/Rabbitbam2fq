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
#include <cassert>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#include "threadconfig.h"


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

//A: 65  G: 71  C: 67  T: 84  N:78
int benckmark = 0;
int is_write = 0;

uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

uint8_t Idx[16] = {0, 0, 2, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 4};
uint8_t IdxRever[16] = {0, 3, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4};

int n_thread = 1;


void read_pack(BGZF *fp, BamBlock *block) {
    pair<bam_block *, int> b;
    b = block->getEmpty();
    int count = 0;
    while (read_block(fp, b.first) == 0) {
        block->inputblock(b.second);
        //printf("read block is %d\n",++count);
        b = block->getEmpty();
    }
    block->ReadComplete();
}

void write_pack(Buffer *buffer) {
    while (!buffer->is_complete()) {
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}

void consumer_pack(BamBlock *block, Buffer *buffer, ThreadConfig *config) {

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    pair<bam_block *, int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *) malloc(sizeof(bam_block));
    int pos = 0;
    uint8_t *quality;
    uint8_t *seq;
    uint8_t *qname;
    pair<char *, int> thread_buffer;
    thread_buffer = buffer->getCap();
    while (1) {
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp = block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second < 0) {
            //printf("%d is Over\n",id);
            break;
        }
        // bam_decode_func
        block_decode_func(comp.first, un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp, b, 0) >= 0) {
            config->TotalReads++;
            if (b->core.flag & 2048) continue;
            config->ReadsQC++;
            if (benckmark)continue;

            if (b->core.flag & 4)config->ReadsAligned--;

            config->TotalBases += b->core.l_qseq;
            seq = bam_get_seq(b);
            if (b->core.flag & 16) {
                for (int i = 0; i < b->core.l_qseq; i++) {
                    config->Bases[IdxRever[bam_seqi(seq, i)]]++;
                }

            } else {
                for (int i = 0; i < b->core.l_qseq; i++)
                    config->Bases[Idx[bam_seqi(seq, i)]]++;
            }


            if (is_write) {
                seq = bam_get_seq(b);
                long long len = strlen(bam_get_qname(b));
                if (pos > buffer->config->Maxn) {
                    buffer->initoutput(thread_buffer.second, pos);
                    thread_buffer = buffer->getCap();
                    pos = 0;
                }
                thread_buffer.first[pos++] = '@';
                memcpy(thread_buffer.first + pos, bam_get_qname(b), len);
                pos += len;
                thread_buffer.first[pos++] = '\n';
                if (b->core.flag & 16) {
                    for (int i = 0; i < b->core.l_qseq; i++)
                        thread_buffer.first[pos + b->core.l_qseq - 1 - i] = BaseRever[bam_seqi(seq, i)];
                    pos += b->core.l_qseq;
                    thread_buffer.first[pos++] = '\n';
                    thread_buffer.first[pos++] = '+';
                    thread_buffer.first[pos++] = '\n';
                    quality = bam_get_qual(b);
                    for (int i = 0; i < b->core.l_qseq; i++) {
                        thread_buffer.first[pos + i] = (char) quality[b->core.l_qseq - i - 1] + 33;
                    }
                    pos += b->core.l_qseq;
                } else {
                    for (int i = 0; i < b->core.l_qseq; i++)
                        thread_buffer.first[pos + i] = Base[bam_seqi(seq, i)];
                    pos += b->core.l_qseq;
                    thread_buffer.first[pos++] = '\n';
                    thread_buffer.first[pos++] = '+';
                    thread_buffer.first[pos++] = '\n';
                    quality = bam_get_qual(b);
                    for (int i = 0; i < b->core.l_qseq; i++) { thread_buffer.first[pos + i] = (char) quality[i] + 33; }
                    pos += b->core.l_qseq;
                }
                thread_buffer.first[pos++] = '\n';
            }


        }
    }
    if (is_write) {
        buffer->initoutput(thread_buffer.second, pos);
        buffer->complete_thread();
    }
}

int main(int argc, char *argv[]) {
    TDEF(fq)
    TSTART(fq)
    printf("Starting Running\n");
    char *in_file = argv[1];
    char *out_file = argv[2];
    n_thread = atoi(argv[3]);
    if (argc > 4)benckmark = 1;
    if (strlen(out_file) > 2)is_write = 1;
    if (benckmark)is_write = 0;
    printf("in : %s\n", in_file);
    printf("out : %s\n", out_file);
    printf("thread : %d\n", n_thread);
    printf("benchmakr %d \n", benckmark);
    printf("is_write %d \n", is_write);
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    ofstream fout;
    fout.open(out_file);
    // /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
    if ((sin = sam_open(in_file, "r")) == NULL) {
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return 0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    BufferConfig config(30, n_thread, 10000000);
    Buffer buffer(&config, &fout);
    BamBlockConfig bamconfig(5000);
    BamBlock block(&bamconfig);

    ThreadConfig **configs = new ThreadConfig *[n_thread];
    for (int t = 0; t < n_thread; t++) {
        configs[t] = new ThreadConfig(t);
    }

    thread producer(&read_pack, sin->fp.bgzf, &block);

    thread **Bam = new thread *[n_thread];
    for (int i = 0; i < n_thread; i++)
        Bam[i] = new thread(&consumer_pack, &block, &buffer, configs[i]);

    thread *writer = NULL;
    if (is_write)writer = new thread(&write_pack, &buffer);

    producer.join();
    for (int i = 0; i < n_thread; i++)
        Bam[i]->join();
    if (is_write)writer->join();

    long long MTotalReads = 0, MReadsQC = 0, MTotalBases = 0, MReadsAligned = 0, MBases[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < n_thread; i++)
        MTotalReads += configs[i]->TotalReads;
    for (int i = 0; i < n_thread; i++)
        MReadsQC += configs[i]->ReadsQC;
    for (int i = 0; i < n_thread; i++)
        MTotalBases += configs[i]->TotalBases;
    for (int i = 0; i < n_thread; i++)
        MReadsAligned += configs[i]->ReadsAligned;
    MReadsAligned += MTotalReads;
    for (int i = 0; i < n_thread; i++)
        for (int j = 0; j < 5; j++) {
            MBases[j] += configs[i]->Bases[j];
        }
    if (benckmark) {
        printf("total read is %lld\n", MTotalReads);
        printf("totol process is %lld\n", MReadsQC);
    } else {
        printf("total read is %lld\n", MTotalReads);
        printf("total process is %lld\n", MReadsQC);
        printf("total base is %lld\n", MTotalBases);
        printf("total reads aligned is %lld\n", MReadsAligned);
        printf("A bases is %lld\n", MBases[0]);
        printf("G bases is %lld\n", MBases[1]);
        printf("C bases is %lld\n", MBases[2]);
        printf("T bases is %lld\n", MBases[3]);
        printf("N bases is %lld\n", MBases[4]);
    }


    sam_close(sin);
    TEND(fq)
    TPRINT(fq, "change time is : ");
    return 0;
}

