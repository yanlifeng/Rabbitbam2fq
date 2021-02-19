#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>

#include <cram/cram.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/hts_log.h>

#ifdef _OPENMP

#include <omp.h>

#endif

uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};
const int MX = 10000000;

struct opts {
    char *fn_ref;
    int flag;
    int clevel;
    int ignore_sam_err;
    int nreads;
    int extra_hdr_nuls;
    int benchmark;
    int nthreads;
    int multi_reg;
    char *index;
    int min_shift;
};

enum test_op {
    READ_COMPRESSED = 1,
    WRITE_BINARY_COMP = 2, // eg bam, bcf
    READ_CRAM = 4,
    WRITE_CRAM = 8,
    WRITE_UNCOMPRESSED = 16,
    WRITE_COMPRESSED = 32, // eg vcf.gz, sam.gz
};

int sam_loop(int argc, char **argv, int optind, struct opts *opts, htsFile *in, htsFile *out, char *moder) {
    int r = 0;
    sam_hdr_t *h = NULL;
    hts_idx_t *idx = NULL;
    bam1_t *b = NULL;

    h = sam_hdr_read(in);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    h->ignore_sam_err = opts->ignore_sam_err;
    if (opts->extra_hdr_nuls > 0) {
        char *new_text = realloc(h->text, h->l_text + opts->extra_hdr_nuls);
        if (new_text == NULL) {
            fprintf(stderr, "Error reallocing header text\n");
            goto fail;
        }
        h->text = new_text;
        memset(&h->text[h->l_text], 0, opts->extra_hdr_nuls);
        h->l_text += opts->extra_hdr_nuls;
    }

    b = bam_init1();
    if (b == NULL) {
        fprintf(stderr, "Out of memory allocating BAM struct\n");
        goto fail;
    }

    /* CRAM output */
    if ((opts->flag & WRITE_CRAM) && opts->fn_ref) {
        // Create CRAM references arrays
        int ret = hts_set_fai_filename(out, opts->fn_ref);

        if (ret != 0)
            goto fail;
    }

    if (!opts->benchmark && sam_hdr_write(out, h) < 0) {
        fprintf(stderr, "Error writing output header.\n");
        goto fail;
    }

    if (opts->index) {
        if (sam_idx_init(out, h, opts->min_shift, opts->index) < 0) {
            fprintf(stderr, "Failed to initialise index\n");
            goto fail;
        }
    }


    clock_t t0 = clock();
    char *outFileName = out->fn;
    FILE *outStream;
    if ((outStream = fopen(outFileName, "w+t")) == NULL) {
        printf("error occur when open %s\n", outFileName);
        return -1;
    }
    int N = 0, M = 0;
    printf("cost %.5f\n", (clock() - t0) * 1e-6);
    t0 = clock();

    int tag = 1;
    if (tag) {
        if ((idx = sam_index_load(in, argv[optind])) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            goto fail;
        }
        printf("cost %.5f\n", (clock() - t0) * 1e-6);
        t0 = clock();
        printf("single region\n");
        printf("optind %d argc %d\n", optind, argc);

//        double cost1 = 0, cost2 = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(1)
#endif
        for (int chr = -2; chr <= 22; chr++) {
            htsFile *ini = hts_open(argv[optind], moder);
            if (ini == NULL) {
                fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
                continue;
            }
            uint8_t *seq;
            uint8_t *qul;
            char *data = (char *) malloc(MX + 1000);
            int32_t lseq;
            char *qname;
            int Ni = 0, Mi = 0;
            int nameLen, pos = 0;
            hts_itr_t *iter;
            char regi[10];
            if (chr == -2)sprintf(regi, "*");
            else if (chr == -1)sprintf(regi, "chrX");
            else if (chr == 0)sprintf(regi, "chrY");
            else
                sprintf(regi, "chr%d", chr);
            if ((iter = sam_itr_querys(idx, h, regi)) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, regi);
                continue;
            }
            int res;
            bam1_t *bi = bam_init1();
            while ((res = sam_itr_next(ini, iter, bi)) >= 0) {
                Ni++;
                if (bi->core.flag & 2048)continue;
                Mi++;
                seq = bam_get_seq(bi);
                qul = bam_get_qual(bi);
                qname = bam_get_qname(bi);
                nameLen = strlen(qname);
                if (pos > MX) {
//                    fwrite(data, sizeof(char), pos, outStream);
                    pos = 0;
                }
                data[pos++] = '@';
                memcpy(data + pos, qname, nameLen);
                pos += nameLen;
                data[pos++] = '\n';
                lseq = bi->core.l_qseq;
                if (bi->core.flag & 16) {
                    for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                        data[pos + j] = BaseRever[bam_seqi(seq, i)];
                    }
                    pos += lseq;
                    data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
                    for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                        data[pos + j] = qul[i] + 33;
                    }
                    pos += lseq;
                } else {
                    for (int i = 0; i < lseq; ++i) {
                        data[pos + i] = Base[bam_seqi(seq, i)];
                    }
                    pos += lseq;
                    data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
                    for (int i = 0; i < lseq; ++i) {
                        data[pos + i] = qul[i] + 33;
                    }
                    pos += lseq;
                }
                data[pos++] = '\n';
            }
//            fwrite(data, sizeof(char), pos, outStream);
            hts_itr_destroy(iter);
            bam_destroy1(bi);
            if (res < -1) {
                fprintf(stderr, "Error reading input.\n");
                continue;
            }
            N += Ni, M += Mi;
            free(data);
        }
        printf("cost %.5f\n", (clock() - t0) * 1e-6);
        t0 = clock();
        hts_idx_destroy(idx);
        idx = NULL;
    } else {

        uint8_t *seq;
        uint8_t *qul;
        char *data = (char *) malloc(MX + 1000);
        int32_t lseq;
        char *qname;
        int nameLen, pos = 0;
        printf("no region\n");
        printf("in %d\n", in->format.format);
        printf("out %d\n", out->format.format);
//        int cnt[1000];
//        for (int i = 0; i < 1000; i++)cnt[i] = 0;
        while (1) {
            struct BGZF *fp = in->fp.bgzf;
            int r = bam_read1(in->fp.bgzf, b);
            if (h && r >= 0) {
                if (b->core.tid >= h->n_targets || b->core.tid < -1 ||
                    b->core.mtid >= h->n_targets || b->core.mtid < -1) {
                    errno = ERANGE;
                    break;
                }
            }
            if (r < 0)break;
            //            cnt[b->core.tid + 500]++;
            N++;
            if (b->core.flag & 2048)continue;
            M++;
            seq = bam_get_seq(b);
            qul = bam_get_qual(b);
            qname = bam_get_qname(b);
            nameLen = strlen(qname);
            if (pos > MX) {
                fwrite(data, sizeof(char), pos, outStream);
                pos = 0;
            }
            data[pos++] = '@';
            memcpy(data + pos, qname, nameLen);
            pos += nameLen;
            data[pos++] = '\n';
            lseq = b->core.l_qseq;
            if (b->core.flag & 16) {
                for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                    data[pos + j] = BaseRever[bam_seqi(seq, i)];
                }
                pos += lseq;
                data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
                for (int i = lseq - 1, j = 0; i >= 0; i--, j++) {
                    data[pos + j] = qul[i] + 33;
                }
                pos += lseq;
            } else {
                for (int i = 0; i < lseq; ++i) {
                    data[pos + i] = Base[bam_seqi(seq, i)];
                }
                pos += lseq;
                data[pos++] = '\n', data[pos++] = '+', data[pos++] = '\n';
                for (int i = 0; i < lseq; ++i) {
                    data[pos + i] = qul[i] + 33;
                }
                pos += lseq;
            }
            data[pos++] = '\n';
        }
        fwrite(data, sizeof(char), pos, outStream);

//        for (int i = 0; i < 1000; i++)
//            if (cnt[i])printf("%d %d\n", i - 500, cnt[i]);
    }

    printf("cost %.5f\n", (clock() - t0) * 1e-6);
    printf("total process %d reads\n", N);
    printf("total print %d reads\n", M);
    if (r < -1) {
        fprintf(stderr, "Error parsing input.\n");
        goto fail;
    }

    if (opts->index) {
        if (sam_idx_save(out) < 0) {
            fprintf(stderr, "Error saving index\n");
            goto fail;
        }
    }

    bam_destroy1(b);
    sam_hdr_destroy(h);

    return 0;
    fail:
    if (b) bam_destroy1(b);
    if (h) sam_hdr_destroy(h);
    if (idx) hts_idx_destroy(idx);

    return 1;
}

int main(int argc, char *argv[]) {
    htsFile *in, *out;
    char moder[8];
    char modew[800];
    int c, exit_code = EXIT_SUCCESS;
    hts_opt *in_opts = NULL, *out_opts = NULL;
    char *out_fn = "-";

    struct opts opts;
    opts.fn_ref = NULL;
    opts.flag = 0;
    opts.clevel = -1;
    opts.ignore_sam_err = 0;
    opts.nreads = 0;
    opts.extra_hdr_nuls = 0;
    opts.benchmark = 0;
    opts.nthreads = 0; // shared pool
    opts.multi_reg = 0;
    opts.index = NULL;
    opts.min_shift = 0;

    while ((c = getopt(argc, argv, "DSIt:i:bzCul:o:N:BZ:@:Mx:m:p:v")) >= 0) {
        switch (c) {
            case 'D':
                opts.flag |= READ_CRAM;
                break;
            case 'S':
                opts.flag |= READ_COMPRESSED;
                break;
            case 'I':
                opts.ignore_sam_err = 1;
                break;
            case 't':
                opts.fn_ref = optarg;
                break;
            case 'i':
                if (hts_opt_add(&in_opts, optarg)) return 1;
                break;
            case 'b':
                opts.flag |= WRITE_BINARY_COMP;
                break;
            case 'z':
                opts.flag |= WRITE_COMPRESSED;
                break;
            case 'C':
                opts.flag |= WRITE_CRAM;
                break;
            case 'u':
                opts.flag |= WRITE_UNCOMPRESSED;
                break; // eg u-BAM not SAM
            case 'l':
                opts.clevel = atoi(optarg);
                break;
            case 'o':
                if (hts_opt_add(&out_opts, optarg)) return 1;
                break;
            case 'N':
                opts.nreads = atoi(optarg);
                break;
            case 'B':
                opts.benchmark = 1;
                break;
            case 'Z':
                opts.extra_hdr_nuls = atoi(optarg);
                break;
            case 'M':
                opts.multi_reg = 1;
                break;
            case '@':
                opts.nthreads = atoi(optarg);
                break;
            case 'x':
                opts.index = optarg;
                break;
            case 'm':
                opts.min_shift = atoi(optarg);
                break;
            case 'p':
                out_fn = optarg;
                break;
            case 'v':
                hts_verbose++;
                break;
        }
    }
    if (argc == optind) {
        fprintf(stderr,
                "Usage: test_view [-DSI] [-t fn_ref] [-i option=value] [-bC] [-l level] [-o option=value] [-N num_reads] [-B] [-Z hdr_nuls] [-@ num_threads] [-x index_fn] [-m min_shift] [-p out] [-v] <in.bam>|<in.sam>|<in.cram> [region]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-D: read CRAM format (mode 'c')\n");
        fprintf(stderr, "-S: read compressed BCF, BAM, FAI (mode 'b')\n");
        fprintf(stderr, "-I: ignore SAM parsing errors\n");
        fprintf(stderr,
                "-t: fn_ref: load CRAM references from the specified fasta file instead of @SQ headers when writing a CRAM file\n");
        fprintf(stderr, "-i: option=value: set an option for CRAM input\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-b: write binary compressed BCF, BAM, FAI (mode 'b')\n");
        fprintf(stderr, "-z: write text compressed VCF.gz, SAM.gz (mode 'z')\n");
        fprintf(stderr, "-C: write CRAM format (mode 'c')\n");
        fprintf(stderr, "-l 0-9: set zlib compression level\n");
        fprintf(stderr, "-o option=value: set an option for CRAM output\n");
        fprintf(stderr, "-N: num_reads: limit the output to the first num_reads reads\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-B: enable benchmarking\n");
        fprintf(stderr, "-M: use hts_itr_multi iterator\n");
        fprintf(stderr, "-Z hdr_nuls: append specified number of null bytes to the SAM header\n");
        fprintf(stderr, "-@ num_threads: use thread pool with specified number of threads\n\n");
        fprintf(stderr, "-x fn: write index to fn\n");
        fprintf(stderr, "-m min_shift: specifies BAI/CSI bin size; 0 is BAI(BAM) or TBI(VCF), 14 is CSI default\n");
        fprintf(stderr, "-p out_fn: output to out_fn instead of stdout\n");
        fprintf(stderr, "-v: increase verbosity\n");
        fprintf(stderr,
                "The region list entries should be specified as 'reg:beg-end', with intervals of a region being disjunct and sorted by the starting coordinate.\n");
        return 1;
    }
    strcpy(moder, "r");
    if (opts.flag & READ_CRAM) strcat(moder, "c");
    else if ((opts.flag & READ_COMPRESSED) == 0) strcat(moder, "b");

    in = hts_open(argv[optind], moder);
    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    strcpy(modew, "w");
    if (opts.clevel >= 0 && opts.clevel <= 9) sprintf(modew + 1, "%d", opts.clevel);
    if (opts.flag & WRITE_CRAM) strcat(modew, "c");
    else if (opts.flag & WRITE_BINARY_COMP) strcat(modew, "b");
    else if (opts.flag & WRITE_COMPRESSED) strcat(modew, "z");
    else if (opts.flag & WRITE_UNCOMPRESSED) strcat(modew, "bu");
    out = hts_open(out_fn, modew);
    if (out == NULL) {
        fprintf(stderr, "Error opening standard output\n");
        return EXIT_FAILURE;
    }

    // Process any options; currently cram only.
    if (hts_opt_apply(in, in_opts))
        return EXIT_FAILURE;
    hts_opt_free(in_opts);

    if (hts_opt_apply(out, out_opts))
        return EXIT_FAILURE;
    hts_opt_free(out_opts);

    // Create and share the thread pool
    htsThreadPool p = {NULL, 0};
    if (opts.nthreads > 0) {
        p.pool = hts_tpool_init(opts.nthreads);
        if (!p.pool) {
            fprintf(stderr, "Error creating thread pool\n");
            exit_code = 1;
        } else {
            hts_set_opt(in, HTS_OPT_THREAD_POOL, &p);
            hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
        }
    }

    int ret;
    switch (hts_get_format(in)->category) {
        case sequence_data:
            ret = sam_loop(argc, argv, optind, &opts, in, out, moder);
            break;


        default:
            fprintf(stderr, "Unsupported or unknown category of data in input file\n");
            return EXIT_FAILURE;
    }

    if (ret != 0)
        exit_code = EXIT_FAILURE;

    ret = hts_close(out);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        exit_code = EXIT_FAILURE;
    }
    ret = hts_close(in);
    if (ret < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = EXIT_FAILURE;
    }

    if (p.pool)
        hts_tpool_destroy(p.pool);

    return exit_code;
}
