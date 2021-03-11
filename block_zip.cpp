//
// Created by 赵展 on 2021/3/5.
//
#include <iostream>
#include <htslib/sam.h>
#include <chrono>
#include <fstream>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include "htslib/khash.h"
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
using namespace std;
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const long long Maxn = 1e8;
char Buffer[Maxn+10000];
long long pos=0;
uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};
struct ZZ_block{
    int errcode=0;
    int count=0;
    unsigned char comp_data[BGZF_MAX_BLOCK_SIZE];
    unsigned int comp_len;
    unsigned char uncomp_data[BGZF_MAX_BLOCK_SIZE];
    unsigned int uncomp_len;
    int64_t block_address;
    int hit_eof;
};
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t)
struct bgzf_cache_t {
    khash_t(cache) *h;
    khint_t last_pos;
};
static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;

    khash_t(cache) *h = fp->cache->h;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, p->size);
    if ( hseek(fp->fp, p->end_offset, SEEK_SET) < 0 )
    {
        // todo: move the error up
        hts_log_error("Could not hseek to %" PRId64, p->end_offset);
        exit(1);
    }
    return p->size;
}
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}
//static int bgzf_uncompress(uint8_t *dst, unsigned int *dlen,const uint8_t *src, unsigned int slen,uint32_t expected_crc)
//{
//        z_stream zs = {
//                .next_in = (Bytef*)src,
//                .avail_in = slen,
//                .total_in = BGZF_MAX_BLOCK_SIZE,
//                .next_out = (Bytef*)dst,
//                .avail_out = *dlen,
//                .total_out = BGZF_MAX_BLOCK_SIZE
//        };
//
//    int ret = inflateInit2(&zs, -15);
//    if (ret != Z_OK) { return -1;}
//    if ((ret = inflate(&zs, Z_FINISH)) != Z_STREAM_END) {
//        if ((ret = inflateEnd(&zs)) != Z_OK) {}
//        return -1;
//    }
//    if ((ret = inflateEnd(&zs)) != Z_OK) {
//        return -1;
//    }
//    *dlen = *dlen - zs.avail_out;
//
//    uint32_t crc = crc32(crc32(0L, NULL, 0L), (unsigned char *)dst, *dlen);
//    if (crc != expected_crc) {
//        hts_log_error("CRC32 checksum mismatch");
//        return -2;
//    }
//    return 0;
//}
//static void *bgzf_decode_func(ZZ_block *arg) {
//    ZZ_block *j = arg;
//
//    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
//    uint32_t crc = le_to_u32((uint8_t *)j->comp_data + j->comp_len-8);
//    int ret = bgzf_uncompress(j->uncomp_data, &j->uncomp_len,
//                              j->comp_data+18, j->comp_len-18, crc);
//    if (ret != 0)
//        j->errcode |= BGZF_ERR_ZLIB;
//
//    return arg;
//}
int read_block(BGZF *fp, ZZ_block *j)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    size = count;
    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }

    compressed_block = (uint8_t*)j->comp_data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    j->count=size;
    j->comp_len = block_length;
    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    j->block_address = block_address;
    j->errcode = 0;

    return 0;
}



int main(){
    TDEF(fq)
    TSTART(fq)
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    if ((sin=sam_open("/Users/zhaozhan/CLionProjects/BAM/SRR_sort.bam", "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    uint8_t* quality;
    uint8_t * seq;
    uint8_t* qname;
    int num=1,ret=0;
    FILE *file=fopen("/Users/zhaozhan/CLionProjects/BAM/my1.txt","w");
    ofstream fout;
    fout.open("./my1.fq");
    int count=0;
    int size=0;
    ZZ_block *j = new ZZ_block;
    while (read_block(sin->fp.bgzf,j) == 0) {
        count++;
        size+=j->count;
       // bgzf_decode_func(j);
        uint32_t x[8];
        memcpy(x,j->uncomp_data+4,32);
        //printf("%d\n",x[0]);
        //break;
    }
    printf("count is %d\n",count);
    printf("size is %d\n",size);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}


/*
 * 代码重构
 * struct{
 * int errcode
 * unsigned char data[]
 * unsigned int data_len
 * int pos
 * }block
 * read——block ： return block_comp
 * bgzf_decode_func: return block_uncomp
 * changeToBam_init1 : return bam1_t
 * bam1_t2fq: return fq
 *
 */