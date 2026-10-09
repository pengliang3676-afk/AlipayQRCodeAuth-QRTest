#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// ===== 从 ALPQRCode.m 原样抄过来的部分 =====
static uint8_t gExp[512];
static uint8_t gLog[256];
static int gInit = 0;

static void GFInit(void) {
    if (gInit) return;
    int x = 1;
    for (int i = 0; i < 255; i++) {
        gExp[i] = (uint8_t)x;
        gLog[x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x100) x ^= 0x11D;
    }
    for (int i = 255; i < 512; i++) gExp[i] = gExp[i - 255];
    gInit = 1;
}

static uint8_t GFMul(uint8_t a, uint8_t b) {
    if (a == 0 || b == 0) return 0;
    return gExp[gLog[a] + gLog[b]];
}

static void RSGenPoly(int degree, uint8_t *gen) {
    memset(gen, 0, degree + 1);
    gen[0] = 1;
    for (int i = 0; i < degree; i++) {
        for (int j = i + 1; j >= 1; j--) {
            gen[j] = gen[j - 1] ^ GFMul(gen[j], gExp[i]);
        }
        gen[0] = GFMul(gen[0], gExp[i]);
    }
}

static void RSEcc(const uint8_t *data, int len, int degree, uint8_t *ecc) {
    uint8_t gen[64];
    RSGenPoly(degree, gen);
    uint8_t rem[64];
    memset(rem, 0, degree);
    for (int i = 0; i < len; i++) {
        uint8_t factor = data[i] ^ rem[0];
        if (degree > 1) memmove(rem, rem + 1, degree - 1);
        rem[degree - 1] = 0;
        for (int j = 0; j < degree; j++) {
            rem[j] ^= GFMul(gen[j + 1], factor);
        }
    }
    memcpy(ecc, rem, degree);
}

typedef struct { int total; int data; int ecc; int blocks; } QRVerInfo;

static const QRVerInfo kVerL[41] = {
    {0,0,0,0},
    {26,19,7,1},{44,34,10,1},{70,55,15,1},{100,80,20,1},{134,108,26,1},
    {172,136,18,2},{196,156,20,2},{242,194,24,2},{292,232,30,2},{346,274,18,4},
    {404,324,20,4},{466,370,24,4},{532,428,26,4},{581,461,30,4},{655,523,22,6},
    {733,589,24,6},{815,647,28,6},{901,721,30,6},{991,795,28,7},{1085,861,28,8},
    {1156,932,28,8},{1258,1006,30,9},{1364,1094,30,9},{1474,1174,30,10},{1588,1276,26,12},
    {1706,1370,28,12},{1828,1468,30,12},{1921,1531,30,13},{2051,1631,30,14},{2185,1735,30,15},
    {2323,1843,30,16},{2465,1955,30,17},{2611,2071,30,18},{2761,2191,30,19},{2876,2306,30,19},
    {3034,2434,30,20},{3196,2566,30,21},{3362,2702,30,22},{3532,2812,30,24},{3706,2956,30,25},
};

static const int kAlign[41][8] = {
    {0},{0},{6,18},{6,22},{6,26},{6,30},{6,34},
    {6,22,38},{6,24,42},{6,26,46},{6,28,50},{6,30,54},{6,32,58},{6,34,62},
    {6,26,46,66},{6,26,48,70},{6,26,50,74},{6,30,54,78},{6,30,56,82},{6,30,58,86},
    {6,34,62,90},{6,28,50,72,94},{6,26,50,74,98},{6,30,54,78,102},{6,28,54,80,106},
    {6,32,58,84,110},{6,30,58,86,114},{6,34,62,90,118},{6,26,50,74,98,122},
    {6,30,54,78,102,126},{6,26,52,78,104,130},{6,30,56,82,108,134},
    {6,34,60,86,112,138},{6,30,58,86,114,142},{6,34,62,90,118,146},
    {6,30,54,78,102,126,150},{6,24,50,76,102,128,154},{6,28,54,80,106,132,158},
    {6,32,58,84,110,136,162},{6,26,54,82,110,138,166},{6,30,58,86,114,142,170},
};
static const int kAlignCount[41] = {
    0,0,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7
};

#define QRNAX 177

typedef struct {
    int size;
    int ver;
    uint8_t mod[QRNAX][QRNAX];
    uint8_t fn[QRNAX][QRNAX];
    uint8_t data[QRNAX][QRNAX];
} QRMat;

static void buildFunction(QRMat *m) {
    int size = m->size;
    memset(m->fn, 0, sizeof(m->fn));
    memset(m->mod, 0, sizeof(m->mod));

    int corners[3][2] = {{0, 0}, {0, size - 7}, {size - 7, 0}};
    for (int k = 0; k < 3; k++) {
        int r0 = corners[k][0], c0 = corners[k][1];
        for (int dr = -1; dr <= 7; dr++) {
            for (int dc = -1; dc <= 7; dc++) {
                int r = r0 + dr, c = c0 + dc;
                if (r < 0 || c < 0 || r >= size || c >= size) continue;
                int v = 0;
                if (dr >= 0 && dr <= 6 && dc >= 0 && dc <= 6) {
                    if (dr == 0 || dr == 6 || dc == 0 || dc == 6) v = 1;
                    else if (dr >= 2 && dr <= 4 && dc >= 2 && dc <= 4) v = 1;
                }
                m->fn[r][c] = 1;
                m->mod[r][c] = v;
            }
        }
    }

    int n = kAlignCount[m->ver];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int r = kAlign[m->ver][i], c = kAlign[m->ver][j];
            if (m->fn[r][c]) continue;
            for (int dr = -2; dr <= 2; dr++) {
                for (int dc = -2; dc <= 2; dc++) {
                    int v = (abs(dr) == 2 || abs(dc) == 2 || (dr == 0 && dc == 0)) ? 1 : 0;
                    m->fn[r + dr][c + dc] = 1;
                    m->mod[r + dr][c + dc] = v;
                }
            }
        }
    }

    for (int i = 8; i < size - 8; i++) {
        int v = (i % 2 == 0) ? 1 : 0;
        if (!m->fn[6][i]) { m->fn[6][i] = 1; m->mod[6][i] = v; }
        if (!m->fn[i][6]) { m->fn[i][6] = 1; m->mod[i][6] = v; }
    }

    for (int i = 0; i < 9; i++) {
        if (!m->fn[8][i]) { m->fn[8][i] = 1; m->mod[8][i] = 0; }
        if (!m->fn[i][8]) { m->fn[i][8] = 1; m->mod[i][8] = 0; }
    }
    for (int i = 0; i < 8; i++) {
        m->fn[8][size - 1 - i] = 1; m->mod[8][size - 1 - i] = 0;
        m->fn[size - 1 - i][8] = 1; m->mod[size - 1 - i][8] = 0;
    }
    m->fn[size - 8][8] = 1;
    m->mod[size - 8][8] = 1;
}

static void placeData(QRMat *m, const uint8_t *bits, long bitCount) {
    int size = m->size;
    memset(m->data, 0, sizeof(m->data));
    int dir = -1;
    int row = size - 1;
    long idx = 0;
    for (int col = size - 1; col > 0; col -= 2) {
        if (col == 6) col--;
        while (1) {
            for (int c = 0; c < 2; c++) {
                int cc = col - c;
                if (m->fn[row][cc]) continue;
                int v = 0;
                if (idx < bitCount) v = bits[idx++];
                m->data[row][cc] = v;
            }
            row += dir;
            if (row < 0 || row >= size) {
                row -= dir;
                dir = -dir;
                break;
            }
        }
    }
}

static int maskFn(int mask, int r, int c) {
    switch (mask) {
        case 0: return (r + c) % 2 == 0;
        case 1: return r % 2 == 0;
        case 2: return c % 3 == 0;
        case 3: return (r + c) % 3 == 0;
        case 4: return (r / 2 + c / 3) % 2 == 0;
        case 5: return (r * c) % 2 + (r * c) % 3 == 0;
        case 6: return ((r * c) % 2 + (r * c) % 3) % 2 == 0;
        case 7: return ((r + c) % 2 + (r * c) % 3) % 2 == 0;
    }
    return 0;
}

static int fmtBits(int mask) {
    int d = (1 << 3) | mask;
    int rem = d;
    for (int i = 0; i < 10; i++) {
        rem = (rem << 1) ^ ((rem >> 9) * 0x537);
    }
    return ((d << 10) | rem) ^ 0x5412;
}

static void applyMask(QRMat *m, int mask) {
    int size = m->size;
    for (int r = 0; r < size; r++) {
        for (int c = 0; c < size; c++) {
            if (m->fn[r][c]) continue;
            m->mod[r][c] = m->data[r][c] ^ (maskFn(mask, r, c) ? 1 : 0);
        }
    }
    int bv = fmtBits(mask);
    for (int i = 0; i < 6; i++)  m->mod[8][i] = (bv >> (14 - i)) & 1;
    m->mod[8][7] = (bv >> 8) & 1;
    m->mod[8][8] = (bv >> 7) & 1;
    m->mod[7][8] = (bv >> 6) & 1;
    for (int i = 9; i < 15; i++) m->mod[14 - i][8] = (bv >> (14 - i)) & 1;
    for (int i = 0; i < 8; i++)  m->mod[size - 1 - i][8] = (bv >> (14 - i)) & 1;
    for (int i = 8; i < 15; i++) m->mod[8][size - 15 + i] = (bv >> (14 - i)) & 1;
    m->mod[size - 8][8] = 1;
}

static int readBack(QRMat *m, uint8_t *outCw, long maxCw, long *outLen) {
    int size = m->size;
    int got[15];
    for (int i = 0; i < 6; i++) got[i] = m->mod[8][i];
    got[6] = m->mod[8][7];
    got[7] = m->mod[8][8];
    got[8] = m->mod[7][8];
    for (int i = 0; i < 6; i++) got[9 + i] = m->mod[5 - i][8];
    int val = 0;
    for (int i = 0; i < 15; i++) val = (val << 1) | got[i];
    val ^= 0x5412;
    int mask = val & 7;

    static uint8_t bits[QRNAX * QRNAX];
    long bitCount = 0;
    int dir = -1;
    int row = size - 1;
    for (int col = size - 1; col > 0; col -= 2) {
        if (col == 6) col--;
        while (1) {
            for (int c = 0; c < 2; c++) {
                int cc = col - c;
                if (m->fn[row][cc]) continue;
                int v = m->mod[row][cc];
                if (maskFn(mask, row, cc)) v ^= 1;
                bits[bitCount++] = v;
            }
            row += dir;
            if (row < 0 || row >= size) {
                row -= dir;
                dir = -dir;
                break;
            }
        }
    }
    long n = 0;
    for (long i = 0; i + 7 < bitCount && n < maxCw; i += 8) {
        int v = 0;
        for (int j = 0; j < 8; j++) v = (v << 1) | bits[i + j];
        outCw[n++] = v;
    }
    *outLen = n;
    return n > 0;
}

static long deinterleave(const uint8_t *cw, long cwLen, int ver, uint8_t *outData) {
    const QRVerInfo *vi = &kVerL[ver];
    int blocks = vi->blocks;
    int dbytes = vi->data;
    int shortLen = dbytes / blocks;
    int numLong = dbytes % blocks;
    int numShort = blocks - numLong;
    int maxD = shortLen + (numLong > 0 ? 1 : 0);

    static uint8_t blkbuf[64][3000];
    int blkLen[64];
    for (int b = 0; b < blocks; b++) {
        blkLen[b] = shortLen + (b >= numShort ? 1 : 0);
    }
    long idx = 0;
    for (int i = 0; i < maxD; i++) {
        for (int b = 0; b < blocks; b++) {
            if (i < blkLen[b] && idx < cwLen) {
                blkbuf[b][i] = cw[idx++];
            }
        }
    }
    long off = 0;
    for (int b = 0; b < blocks; b++) {
        memcpy(outData + off, blkbuf[b], blkLen[b]);
        off += blkLen[b];
    }
    return off;
}

// ===== 主流程：和 ALPQRCode.m 的 matrixWithText 一致，但输出矩阵 =====
static int buildMatrix(const char *text, uint8_t outMod[QRNAX][QRNAX], int *outSize, int *outVer, int *outMask, int verbose) {
    GFInit();
    int n = (int)strlen(text);
    const uint8_t *src = (const uint8_t *)text;

    int ver = 0;
    for (int v = 1; v <= 40; v++) {
        int cap = kVerL[v].data * 8;
        int overhead = 4 + (v < 10 ? 8 : 16);
        if (n * 8 + overhead <= cap) { ver = v; break; }
    }
    if (ver == 0) return 0;
    const QRVerInfo *vi = &kVerL[ver];

    static uint8_t bits[QRNAX * QRNAX];
    long bc = 0;
#define PUT(val, cnt) for (int _i = (cnt) - 1; _i >= 0; _i--) bits[bc++] = ((val) >> _i) & 1
    PUT(4, 4);
    PUT(n, ver < 10 ? 8 : 16);
    for (int i = 0; i < n; i++) PUT(src[i], 8);
    long totalDataBits = (long)vi->data * 8;
    for (int i = 0; i < 4 && bc < totalDataBits; i++) bits[bc++] = 0;
    while (bc % 8 != 0) bits[bc++] = 0;
    {
        int p = 0;
        uint8_t pad[2] = {0xEC, 0x11};
        while (bc < totalDataBits) { PUT(pad[p % 2], 8); p++; }
    }
#undef PUT

    int totalDataBytes = vi->data;
    uint8_t *dcw = calloc(totalDataBytes, 1);
    for (long i = 0; i < totalDataBits; i++) {
        if (bits[i]) dcw[i / 8] |= (1 << (7 - (i % 8)));
    }

    int blocks = vi->blocks;
    int eccLen = vi->ecc;
    int shortLen = totalDataBytes / blocks;
    int numLong = totalDataBytes % blocks;
    int numShort = blocks - numLong;
    int maxD = shortLen + (numLong > 0 ? 1 : 0);

    static uint8_t dblk[64][3000], eblk[64][64];
    int dlen[64];
    int off = 0;
    for (int b = 0; b < blocks; b++) {
        dlen[b] = shortLen + (b >= numShort ? 1 : 0);
        memcpy(dblk[b], dcw + off, dlen[b]);
        off += dlen[b];
        RSEcc(dblk[b], dlen[b], eccLen, eblk[b]);
    }

    uint8_t *inter = calloc(vi->total, 1);
    long ip = 0;
    for (int i = 0; i < maxD; i++) {
        for (int b = 0; b < blocks; b++) {
            if (i < dlen[b]) inter[ip++] = dblk[b][i];
        }
    }
    for (int i = 0; i < eccLen; i++) {
        for (int b = 0; b < blocks; b++) inter[ip++] = eblk[b][i];
    }
    free(dcw);

    uint8_t *ibit = calloc(vi->total * 8, 1);
    long ibc = 0;
    for (long i = 0; i < ip; i++) {
        for (int k = 7; k >= 0; k--) ibit[ibc++] = (inter[i] >> k) & 1;
    }

    QRMat m;
    memset(&m, 0, sizeof(m));
    m.size = 17 + 4 * ver;
    m.ver = ver;
    buildFunction(&m);
    placeData(&m, ibit, ibc);

    static uint8_t backCw[QRNAX * QRNAX];
    static uint8_t backData[3000];
    int chosen = -1;
    for (int mask = 0; mask < 8; mask++) {
        /* 0 = 没有选到 */
        QRMat *t = malloc(sizeof(QRMat));
        memcpy(t, &m, sizeof(QRMat));
        applyMask(t, mask);
        long cwLen = 0;
        if (!readBack(t, backCw, sizeof(backCw), &cwLen)) { free(t); continue; }
        long dLen = deinterleave(backCw, cwLen, ver, backData);
        if (dLen < totalDataBytes) { free(t); continue; }
        uint8_t *bs = calloc((size_t)dLen * 8 + 8, 1);
        long bsp = 0;
        for (long i = 0; i < dLen; i++) {
            for (int k = 7; k >= 0; k--) bs[bsp++] = (backData[i] >> k) & 1;
        }
        int mode = 0;
        for (int i = 0; i < 4; i++) mode = (mode << 1) | bs[i];
        int cntLen = (ver < 10) ? 8 : 16;
        int cnt = 0;
        for (int i = 4; i < 4 + cntLen; i++) cnt = (cnt << 1) | bs[i];
        int ok = (mode == 4 && cnt == n);
        if (ok) {
            for (int i = 0; i < n && ok; i++) {
                int v = 0;
                long base = 4 + cntLen + (long)i * 8;
                for (int j = 0; j < 8; j++) v = (v << 1) | bs[base + j];
                if (v != src[i]) ok = 0;
            }
        }
        free(bs);
        if (ok) {
            if (verbose) printf("     掩模 %d 自检通过\n", mask);
            chosen = mask;
            memcpy(outMod, t->mod, sizeof(t->mod));
            *outSize = t->size;
            *outVer = ver;
            *outMask = mask;
            free(t);
            break;
        }
        free(t);
    }
    free(inter);
    free(ibit);
    if (chosen < 0) {
        if (verbose) printf("     !! 8 个掩模自检全部失败\n");
        return 0;
    }
    return 1;
}


int main(void) {
    // 内容：真实结构，866 字符
    static char text[4096];
    const char *prefix =
        "https://wappaygw.alipay.com/home/exterfaceAssign.htm?"
        "pid=2088631028484361&app_id=2021002165619311&timestamp=2026-10-09+22%3A22%3A20"
        "&charset=utf-8&format=json&sign_type=RSA2"
        "&alipay_root_cert_sn=687b59193f3f462dd5336e5abf83c5d8_02941eef3187dddf3d3b83462e1dfcf6"
        "&app_cert_sn=408b30c8c37959ecb86a25281e607ccb"
        "&apiname=com.alipay.account.auth&method=alipay.open.auth.sdk.code.get"
        "&app_name=mc&biz_type=openservice&product_id=APP_FAST_LOGIN&scope=auth_user"
        "&target_id=459b9c43-f1f6-4a15-b5d1-666f9305f396&auth_type=AUTHACCOUNT&sign=";
    int plen = (int)strlen(prefix);
    memcpy(text, prefix, plen);
    for (int i = 0; i < 344; i++) text[plen + i] = 'A';
    text[plen + 344] = 0;

    printf("内容长度: %d\n", (int)strlen(text));
    fflush(stdout);

    static uint8_t mod[QRNAX][QRNAX];
    int size = 0, ver = 0, mask = -1;
    int ok = buildMatrix(text, mod, &size, &ver, &mask, 1);
    printf("结果: %s  版本=%d 尺寸=%d 掩模=%d\n", ok ? "成功" : "失败", ver, size, mask);
    fflush(stdout);
    if (!ok) return 1;

    FILE *f = fopen("objc_qr.pbm", "wb");
    int quiet = 4;
    int total = size + quiet * 2;
    fprintf(f, "P1\n%d %d\n", total, total);
    for (int r = 0; r < total; r++) {
        for (int c = 0; c < total; c++) {
            int rr = r - quiet, cc = c - quiet;
            int v = 0;
            if (rr >= 0 && cc >= 0 && rr < size && cc < size) v = mod[rr][cc];
            fprintf(f, "%d ", v);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    printf("已写出 objc_qr.pbm (%dx%d)\n", total, total);

    f = fopen("objc_qr_content.txt", "wb");
    fwrite(text, 1, strlen(text), f);
    fclose(f);
    printf("已写出 objc_qr_content.txt\n");
    fflush(stdout);
    return 0;
}
