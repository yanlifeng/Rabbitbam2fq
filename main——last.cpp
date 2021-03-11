#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <htslib/sam.h>
#include <iostream>

enum filetype {
	FBAM = 1,   // BAM file
	FSAM = 2,   // SAM file
};

int ftype;

int sam_test_extract(int argc, char **argv, int optind, htsFile *in, htsFile *out) {
		sam_hdr_t *hdr;
		bam1_t *b;
		hts_idx_t *idx = NULL;
		hts_itr_t *iter = NULL;
		int ret;

		if ((hdr = sam_hdr_read(in)) == NULL) {
			fprintf(stderr, "[E::%s] couldn't read header for '%s'\n", __func__, argv[optind]);
			return  -1;
		}
		if ((b = bam_init1()) == NULL) {
			fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
			goto fail;
		}
		if (ftype == FBAM && optind + 2 <= argc) { // BAM input and has a region.
			if ((idx = sam_index_load(in, argv[optind])) == 0) {
				fprintf(stderr, "[E::%s] fail to load the index for '%s'\n", __func__, argv[optind]);
				goto fail;
			}
			if ((iter = sam_itr_querys(idx, hdr, argv[optind + 1])) == 0) {
				fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[optind + 1]);
				goto fail;
			}
			while ((ret = sam_itr_next(in, iter, b)) >= 0) {
				if (sam_write1(out, hdr, b) < 0) {
					fprintf(stderr, "[E::%s] Error writing output.\n", __func__);
					goto fail;
				}
			}
			if (ret < -1) {
				fprintf(stderr, "[E::%s] Error reading input.\n", __func__);
				goto fail;
			}
			hts_itr_destroy(iter);
			iter = NULL;
			hts_idx_destroy(idx);
			idx = NULL;
		} else if (optind + 2 > argc) {
		    printf("It is in this\n");
		    int num=0;
			while ((ret = sam_read1(in, hdr, b)) >= 0 && num++<10) {
//			    printf("%s\n",b->data); //名字
                printf("%s\n",hdr->sdict);
				if (sam_write1(out, hdr, b) < 0) {
					fprintf(stderr,"[E::%s] Error writing alignments.\n", __func__);
					goto fail;
				}
			}
			if (ret < -1) {
				fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
				goto fail;
			}
		} else { // SAM input and has a region.
			fprintf(stderr, "[E::%s] couldn't extract alignments directly from raw sam file.\n", __func__);
			goto fail;
		}
		bam_destroy1(b);
		sam_hdr_destroy(hdr);
		return  0;

	fail:
		if (iter) sam_itr_destroy(iter);
		if (b) bam_destroy1(b);
		if (idx) hts_idx_destroy(idx);
		if (hdr) sam_hdr_destroy(hdr);
		return  1;
}

int main(int argc, char **argv) {

	htsFile *in, *out;
	int c, ret, exit_code;
	char moder[8];
	//char modew[800];
	char *outfn = "-";

	ftype = FSAM;
	exit_code = 0;
	strcpy(moder, "r");
	while ((c = getopt(argc, argv, "bo:")) >= 0) {
		switch (c) {
            case 'b': strcat(moder, "b"); ftype = FBAM; break;
			case 'o': outfn = optarg; break;
			          		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: %s [-b] [-o out.sam] <in.bam>|<in.sam> [region]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "\t-b:\tUse BAM as input if this option is set, otherwise use SAM as input.\n");
		fprintf(stderr, "\t-o:\tPath to the output file. Output to stdout if this option is not set.\n");
		return  -1;
	}
	if ((in = hts_open(argv[optind], moder)) == NULL) {
		fprintf(stderr, "Error opening '%s'\n", argv[1]);
		return  -3;
	}
	if ((out = hts_open(outfn, "w")) == NULL) {
		fprintf(stderr, "Error opening '%s'\n", argv[2]);
		return  -3;
	}


	if ((ret = sam_test_extract(argc, argv, optind, in, out)) != 0) {
		fprintf(stderr, "Error extracting alignment from '%s'\n", argv[optind]);
		exit_code = -5;
	}


	if ((ret = hts_close(out)) < 0) {
		fprintf(stderr, "Error closing output.\n");
		exit_code = -3;
	}
	if ((ret = hts_close(in)) < 0) {
		fprintf(stderr, "Error closing input.\n");
		exit_code = -3;
	}
	return exit_code;
}