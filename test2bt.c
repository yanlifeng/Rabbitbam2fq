#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef unsigned char uint8_t;
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

static inline void solve(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
            "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
            "A=AAACAMAGARASAVATAWAYAHAKADABAN"
            "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
            "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
            "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
            "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
            "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
            "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
            "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
            "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
            "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
            "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
            "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
            "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
            "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
            "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len / 2;
    seq[0] = 0;

    for (i = 0; i < len2; i++) {
        printf("%c%c", code2base[(size_t) nib[i] * 2], code2base[(size_t) nib[i] * 2 + 1]);
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i * 2], &code2base[(size_t) nib[i] * 2], 2);
    }
    printf("\n");

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

static inline void solveRE(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
            "==T=G=M=C=R=S=V=A=W=Y=H=K=D=B=N="
            "=TTTGTMTCTRTSTVTATWTYTHTKTDTBTNT"
            "=GTGGGMGCGRGSGVGAGWGYGHGKGDGBGNG"
            "=MTMGMMMCMRMSMVMAMWMYMHMKMDMBMNM"
            "=CTCGCMCCCRCSCVCACWCYCHCKCDCBCNC"
            "=RTRGRMRCRRRSRVRARWRYRHRKRDRBRNR"
            "=STSGSMSCSRSSSVSASWSYSHSKSDSBSNS"
            "=VTVGVMVCVRVSVVVAVWVYVHVKVDVBVNV"
            "=ATAGAMACARASAVAAAWAYAHAKADABANA"
            "=WTWGWMWCWRWSWVWAWWWYWHWKWDWBWNW"
            "=YTYGYMYCYRYSYVYAYWYYYHYKYDYBYNY"
            "=HTHGHMHCHRHSHVHAHWHYHHHKHDHBHNH"
            "=KTKGKMKCKRKSKVKAKWKYKHKKKDKBKNK"
            "=DTDGDMDCDRDSDVDADWDYDHDKDDDBDND"
            "=BTBGBMBCBRBSBVBABWBYBHBKBDBBBNB"
            "=NTNGNMNCNRNSNVNANWNYNHNKNDNBNNN";
    int i, len2 = len / 2;
    seq[0] = 0;

    for (i = 0; i < len2; i++) {
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i * 2 + (len & 1)], &code2base[(size_t) nib[len2 - 1 - i] * 2], 2);
    }

    if ((i *= 2) < len)
        seq[0] = seq_nt16_str[bam_seqi(nib, i)];
}

int main() {

    int srci[] = {68, 40, 132, 33, 34, 40, 40, 65, 20, 33, 72, 65, 40, 68, 20, 136, 72, 17, 40, 68, 72, 33, 24, 136, 65,
                  66, 129, 33, 34, 132, 65, 66, 132, 17, 66, 20, 132, 66, 33, 34, 130, 34, 34, 17, 17, 17, 20, 66, 34,
                  33, 34, 34, 34, 34, 34, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
                  38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38,
                  38};
    int stdi[] = {71, 71, 67, 84, 84, 71, 67, 65, 67, 67, 67, 84, 67, 84, 71, 65, 65, 71, 67, 65, 71, 84, 71, 65, 67,
                  84, 71, 71, 65, 71, 84, 84, 71, 84, 65, 65, 67, 84, 71, 71, 71, 84, 67, 65, 65, 84, 84, 84, 71, 65,
                  71, 67, 84, 65, 67, 65, 67, 67, 84, 71, 71, 65, 71, 67, 84, 71, 65, 65, 71, 67, 65, 71, 84, 71, 71,
                  67, 67, 65, 67, 67, 84, 67, 67, 67, 67, 67, 65, 65, 65, 65, 65, 65, 65, 71, 71, 67, 67, 67, 67, 65};
    int len = 99;
    char *src = malloc(len);
    for (int i = 0; i < len; i++)
        src[i] = srci[i];
    src[len] = 0;
//    printf("%s\n", src);
    char *tar = malloc(len);
    solveRE(src, tar, len);
    for (int i = 0; i < len; i++)
        printf("%c", tar[len - 1 - i]);
    printf("\n");
    for (int i = 0; i < len; i++)
        printf("%c", stdi[i]);
    printf("\n");
    int ok = 1;
    for (int i = 0; i < len; i++) {
        if (tar[i] != stdi[len - 1 - i]) {
            ok = 0;
            break;
        }
    }
    if (ok)printf("ac\n");
    else printf("wa\n");
    return 0;
}