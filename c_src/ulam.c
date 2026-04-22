/*
 * Ulam number generator - Gibbs residue-based algorithm.
 *
 * Reference article that the algorithm is based on: Philip Gibbs, "An efficient method for computing Ulam numbers" arXiv:1508.0085v2
 * 
 * Instead of storing each Ulam number as a uint64 (8 bytes), it is stored as the
 * deviation from the linear trend:
 *
 *   a_off[i] = ulam(i+1) - (int64_t)((i+1) * STEP) [int32, 4 bytes]
 *
 * where STEP = 13.5178 is the mean spacing between consecutive Ulam numbers.
 * The Gibbs quasi-periodicity theorem guarantees these offsets are always
 * small - a few thousand at most - so int32_t is sufficient for all n.
 *
 * For n = 20 billion this reduces the sequence-array memory from around 12 GB
 * (uint64) to around 6 GB (int32), making the computation feasible on 32 GB RAM.
 *
 * Save format (VERSION 4, little-endian):
 *   magic(uint32)  version(uint32=4)  max_computed(uint64)  count(uint64)
 *   a_off[0..count-1]                         (int32 * count)
 *   k_cap(uint64)  k[0..k_cap-1]             (uint32 * k_cap)
 *   track_pairs(int32)
 *   [pairs_count(uint64) pairs[0..pairs_count*2-1] (uint64 * pairs_count*2)]
 * 
 */

#include "ulam.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#define CHUNK_ELEMS(esz) ((size_t)(64 * 1024 * 1024 / (esz)))

static int fwrite_all(const void *buf, size_t esz, size_t n, FILE *f) {
    const char *p = (const char *)buf;
    size_t chunk = CHUNK_ELEMS(esz);
    while (n > 0) {
        size_t batch = n < chunk ? n : chunk;
        if (fwrite(p, esz, batch, f) != batch) return 0;
        p += batch * esz;
        n -= batch;
    }
    return 1;
}

#define MAGIC    0x554C414DU
#define VERSION  4U
#define VERSION3 3U
#define VERSION2 2U

#define LAMBDA  2.44344296778474
#define STEP    13.517831473

#define INIT_ULAM_CAP  4096
#define INIT_PAIR_CAP  512

static const uint32_t POW2_30[30] = {
    1u<< 0, 1u<< 1, 1u<< 2, 1u<< 3, 1u<< 4, 1u<< 5,
    1u<< 6, 1u<< 7, 1u<< 8, 1u<< 9, 1u<<10, 1u<<11,
    1u<<12, 1u<<13, 1u<<14, 1u<<15, 1u<<16, 1u<<17,
    1u<<18, 1u<<19, 1u<<20, 1u<<21, 1u<<22, 1u<<23,
    1u<<24, 1u<<25, 1u<<26, 1u<<27, 1u<<28, 1u<<29
};

struct UlamState {
    int32_t  *a_off;
    size_t    count;
    size_t    ulam_cap;
    uint32_t *k;
    size_t    k_cap;
    uint32_t *nx;
    uint32_t *pv;
    uint32_t *res_index;
    size_t    nindex;

    uint64_t  max_computed;
    uint64_t  prev_ulam;
    FILE     *adj_fp;

    int       track_pairs;
    uint64_t *pairs;
    size_t    pairs_count;
    size_t    pairs_cap;
};

static double dmod(uint64_t x, double m) {
    double z = (double)x / m;
    long long iz = (long long)z;
    z = (z - (double)iz) * m;
    if (z < 0.0) z += m;
    return z;
}

static inline uint64_t get_ulam(const UlamState *s, uint32_t seqpos) {
    if (seqpos == 0) return 0ULL;
    int64_t ground = (int64_t)((double)seqpos * STEP);
    return (uint64_t)(ground + (int64_t)s->a_off[seqpos - 1]);
}

static inline int is_ulam_bit(const UlamState *s, uint64_t v) {
    size_t w = (size_t)(v / 30);
    if (w >= s->k_cap) return 0;
    return (s->k[w] & POW2_30[v % 30]) != 0;
}

static inline void set_ulam_bit(UlamState *s, uint64_t v) {
    s->k[v / 30] |= POW2_30[v % 30];
}

static inline void store_ulam(UlamState *s, size_t slot, uint64_t val) {
    int64_t ground = (int64_t)((double)(slot + 1) * STEP);
    s->a_off[slot] = (int32_t)((int64_t)val - ground);
}

static int grow_k(UlamState *s, uint64_t max_val) {
    size_t needed = (size_t)(max_val / 30) + 2;
    if (needed <= s->k_cap) return 1;
    size_t nc = s->k_cap ? s->k_cap : 256;
    while (nc < needed) nc *= 2;
    uint32_t *p = realloc(s->k, nc * sizeof(uint32_t));
    if (!p) return 0;
    memset(p + s->k_cap, 0, (nc - s->k_cap) * sizeof(uint32_t));
    s->k     = p;
    s->k_cap = nc;
    return 1;
}

static void rebuild_res_index(UlamState *s) {
    memset(s->res_index, 0, s->nindex * sizeof(uint32_t));
    uint32_t last = 0;
    uint32_t i    = s->nx[0];
    for (int j = 0; j < (int)s->nindex; j++) {
        double thr = (double)j / (double)s->nindex;
        while (i != 0 && dmod(get_ulam(s, i), LAMBDA) / LAMBDA < thr) {
            last = i;
            i    = s->nx[i];
        }
        s->res_index[j] = last;
    }
}

static int grow_seq(UlamState *s, size_t needed) {
    if (needed <= s->ulam_cap) return 1;
    size_t nc = s->ulam_cap ? s->ulam_cap : INIT_ULAM_CAP;
    while (nc < needed) nc *= 2;

    int32_t *pa = realloc(s->a_off, nc * sizeof(int32_t));
    if (!pa) return 0;
    s->a_off = pa;

    uint32_t *pnx = realloc(s->nx, (nc + 1) * sizeof(uint32_t));
    if (!pnx) return 0;
    s->nx = pnx;

    uint32_t *ppv = realloc(s->pv, (nc + 1) * sizeof(uint32_t));
    if (!ppv) return 0;
    s->pv = ppv;

    size_t new_ni = nc / 100 + 1;
    int nindex_grew = (new_ni > s->nindex);
    if (nindex_grew) {
        uint32_t *pi = realloc(s->res_index, new_ni * sizeof(uint32_t));
        if (!pi) return 0;
        s->res_index = pi;
        s->nindex    = new_ni;
    }

    s->ulam_cap = nc;

    if (nindex_grew)
        rebuild_res_index(s);

    return 1;
}

static int pairs_push(UlamState *s, uint64_t a, uint64_t b) {
    if (s->pairs_count >= s->pairs_cap) {
        size_t nc = s->pairs_cap ? s->pairs_cap * 2 : INIT_PAIR_CAP;
        uint64_t *p = realloc(s->pairs, nc * 2 * sizeof(uint64_t));
        if (!p) return 0;
        s->pairs     = p;
        s->pairs_cap = nc;
    }
    s->pairs[s->pairs_count * 2]     = a;
    s->pairs[s->pairs_count * 2 + 1] = b;
    s->pairs_count++;
    return 1;
}

static void setlinks(UlamState *s, uint32_t n) {
    double rdn = dmod(get_ulam(s, n), LAMBDA) / LAMBDA;

    int j = (int)((double)s->nindex * rdn);
    if (j >= (int)s->nindex) j = (int)s->nindex - 1;
    if (j < 0)               j = 0;

    uint32_t pvi = s->res_index[j];
    for (;;) {
        uint32_t nxt = s->nx[pvi];
        if (nxt == 0) break;
        if (dmod(get_ulam(s, nxt), LAMBDA) / LAMBDA >= rdn) break;
        pvi = nxt;
    }

    uint32_t nxi = s->nx[pvi];
    s->pv[n]   = pvi;
    s->nx[pvi] = n;
    s->nx[n]   = nxi;
    s->pv[nxi] = n;

    j++;
    while (j < (int)s->nindex) {
        uint32_t cur = s->res_index[j];
        if (cur != 0 && dmod(get_ulam(s, cur), LAMBDA) / LAMBDA >= rdn) break;
        s->res_index[j] = n;
        j++;
    }
}

static int ulam_register(UlamState *s, uint64_t val) {
    if (!grow_seq(s, s->count + 2)) return 0;

    if (s->track_pairs && s->count > 0)
        if (!pairs_push(s, s->prev_ulam, val)) return 0;

    if (s->adj_fp && s->prev_ulam + 1 == val)
        fprintf(s->adj_fp, "%" PRIu64 ",%" PRIu64 "\n", s->prev_ulam, val);

    s->prev_ulam = val;
    store_ulam(s, s->count, val);
    s->count++;
    set_ulam_bit(s, val);
    setlinks(s, (uint32_t)s->count);
    return 1;
}

UlamState *ulam_init(int track_pairs) {
    UlamState *s = calloc(1, sizeof *s);
    if (!s) return NULL;

    s->ulam_cap = INIT_ULAM_CAP;
    s->a_off    = malloc(s->ulam_cap * sizeof(int32_t));
    if (!s->a_off) { free(s); return NULL; }

    s->nx = calloc(s->ulam_cap + 1, sizeof(uint32_t));
    s->pv = calloc(s->ulam_cap + 1, sizeof(uint32_t));
    if (!s->nx || !s->pv) { ulam_free(s); return NULL; }

    s->nindex    = s->ulam_cap / 100 + 1;
    s->res_index = calloc(s->nindex, sizeof(uint32_t));
    if (!s->res_index) { ulam_free(s); return NULL; }

    s->k_cap = 256;
    s->k     = calloc(s->k_cap, sizeof(uint32_t));
    if (!s->k) { ulam_free(s); return NULL; }

    s->track_pairs = track_pairs;
    s->adj_fp      = NULL;
    s->prev_ulam   = 0;

    store_ulam(s, 0, 1);
    s->count = 1;
    set_ulam_bit(s, 1);
    setlinks(s, 1);
    s->prev_ulam = 1;

    if (!ulam_register(s, 2)) { ulam_free(s); return NULL; }

    s->max_computed = 2;
    return s;
}

void ulam_free(UlamState *s) {
    if (!s) return;
    if (s->adj_fp) { fclose(s->adj_fp); s->adj_fp = NULL; }
    free(s->a_off);
    free(s->nx);
    free(s->pv);
    free(s->res_index);
    free(s->k);
    free(s->pairs);
    free(s);
}

int ulam_compute_up_to(UlamState *s, uint64_t n) {
    if (n <= s->max_computed) return 1;

    if (!grow_k(s, n)) return 0;

    for (uint64_t a0 = s->max_computed + 1; a0 <= n; a0++) {
        double   rd0  = dmod(a0, LAMBDA) / LAMBDA;
        int      ulam = 0;
        uint64_t a1x  = 0;

        if (rd0 < 0.24 || rd0 > 0.80) {
            uint32_t j  = (uint32_t)s->count;
            uint64_t aj = get_ulam(s, j);
            while (2 * aj > a0 && j > 0) {
                uint64_t a2 = a0 - aj;
                if (is_ulam_bit(s, a2)) {
                    if (ulam) { ulam = 0; break; }
                    ulam = 1;
                    a1x  = a2;
                }
                j--;
                aj = get_ulam(s, j);
            }
        } else {
            int more = 1;

            uint32_t i  = s->nx[0];
            uint64_t ai = get_ulam(s, i);
            double   rdi = (i != 0) ? dmod(ai, LAMBDA) / LAMBDA : 0.0;

            while (more && i != 0 && 2.0 * rdi <= rd0 + 0.0002) {
                uint64_t a2 = a0 - ai;
                if (is_ulam_bit(s, a2) && ai != a2 && a2 != a1x) {
                    if (ulam) { ulam = 0; more = 0; }
                    else      { ulam = 1; a1x = ai; }
                }
                i   = s->nx[i];
                ai  = get_ulam(s, i);
                rdi = (i != 0) ? dmod(ai, LAMBDA) / LAMBDA : 0.0;
            }

            i   = s->pv[0];
            ai  = get_ulam(s, i);
            rdi = (i != 0) ? dmod(ai, LAMBDA) / LAMBDA : 1.0;

            while (more && i != 0 && 2.0 * (1.0 - rdi) <= (1.0 - rd0) + 0.0002) {
                uint64_t a2 = a0 - ai;
                if (is_ulam_bit(s, a2) && ai != a2 && a2 != a1x) {
                    if (ulam) { ulam = 0; more = 0; }
                    else      { ulam = 1; a1x = ai; }
                }
                i   = s->pv[i];
                ai  = get_ulam(s, i);
                rdi = (i != 0) ? dmod(ai, LAMBDA) / LAMBDA : 1.0;
            }
        }

        if (ulam)
            if (!ulam_register(s, a0)) return 0;
    }

    s->max_computed = n;
    return 1;
}

int ulam_is_ulam(UlamState *s, uint64_t n) {
    if (n > s->max_computed) return -1;
    return is_ulam_bit(s, n) ? 1 : 0;
}

uint64_t ulam_max_computed(UlamState *s) { return s->max_computed; }
uint64_t ulam_count(UlamState *s)        { return (uint64_t)s->count; }
int      ulam_track_pairs(UlamState *s)  { return s->track_pairs; }
uint64_t ulam_pairs_count(UlamState *s)  { return (uint64_t)s->pairs_count; }
uint64_t ulam_pairs_ptr_raw(UlamState *s){ return (uint64_t)(uintptr_t)s->pairs; }

uint64_t *ulam_get_ptr(UlamState *s) { (void)s; return NULL; }
uint64_t *ulam_pairs_ptr(UlamState *s) { return s->pairs; }

uint64_t ulam_get(UlamState *s, uint64_t idx) {
    return idx < s->count ? get_ulam(s, (uint32_t)(idx + 1)) : 0;
}

void ulam_get_pair(UlamState *s, uint64_t idx, uint64_t *a, uint64_t *b) {
    if (idx < s->pairs_count) {
        *a = s->pairs[idx * 2];
        *b = s->pairs[idx * 2 + 1];
    } else {
        *a = *b = 0;
    }
}

uint64_t *ulam_find_adjacent_pairs(UlamState *s, uint64_t *out_count) {
    size_t cap = 1024, cnt = 0;
    uint64_t *res = malloc(cap * 2 * sizeof(uint64_t));
    if (!res) { *out_count = 0; return NULL; }

    for (size_t i = 0; i + 1 < s->count; i++) {
        uint64_t a = get_ulam(s, (uint32_t)(i + 1));
        uint64_t b = get_ulam(s, (uint32_t)(i + 2));
        if (b == a + 1) {
            if (cnt >= cap) {
                cap *= 2;
                uint64_t *p = realloc(res, cap * 2 * sizeof(uint64_t));
                if (!p) { free(res); *out_count = 0; return NULL; }
                res = p;
            }
            res[cnt * 2]     = a;
            res[cnt * 2 + 1] = b;
            cnt++;
        }
    }

    *out_count = (uint64_t)cnt;
    if (cnt == 0) { free(res); return NULL; }
    return res;
}

void ulam_free_ptr(void *p) { free(p); }

int ulam_open_adj_stream(UlamState *s, const char *path) {
    if (s->adj_fp) { fclose(s->adj_fp); s->adj_fp = NULL; }
    s->adj_fp = fopen(path, "w");
    if (!s->adj_fp) return 0;
    fprintf(s->adj_fp, "n,n+1\n");
    for (size_t i = 0; i + 1 < s->count; i++) {
        uint64_t a = get_ulam(s, (uint32_t)(i + 1));
        uint64_t b = get_ulam(s, (uint32_t)(i + 2));
        if (b == a + 1)
            fprintf(s->adj_fp, "%" PRIu64 ",%" PRIu64 "\n", a, b);
    }
    return 1;
}

void ulam_close_adj_stream(UlamState *s) {
    if (s->adj_fp) { fclose(s->adj_fp); s->adj_fp = NULL; }
}

int ulam_save(UlamState *s, const char *path) {
    /* Write to a temp file first; rename to path only on success so a
       failed/interrupted save never corrupts the existing file. */
    char tmp[4096];
    if ((size_t)snprintf(tmp, sizeof tmp, "%s.tmp", path) >= sizeof tmp) return 0;

    FILE *f = fopen(tmp, "wb");
    if (!f) return 0;

    uint32_t hdr[2] = { MAGIC, VERSION };
    if (fwrite(hdr,              4, 2, f) != 2)               goto fail;
    if (fwrite(&s->max_computed, 8, 1, f) != 1)               goto fail;
    uint64_t cnt = (uint64_t)s->count;
    if (fwrite(&cnt,             8, 1, f) != 1)               goto fail;
    if (!fwrite_all(s->a_off,    4, s->count, f))             goto fail;

    uint64_t kc = (uint64_t)s->k_cap;
    if (fwrite(&kc,              8, 1, f) != 1)               goto fail;
    if (!fwrite_all(s->k,        4, s->k_cap, f))             goto fail;

    int32_t tp = (int32_t)s->track_pairs;
    if (fwrite(&tp, 4, 1, f) != 1) goto fail;
    if (s->track_pairs) {
        uint64_t pc = (uint64_t)s->pairs_count;
        if (fwrite(&pc, 8, 1, f) != 1) goto fail;
        if (s->pairs_count > 0 &&
            !fwrite_all(s->pairs, 8, s->pairs_count * 2, f))  goto fail;
    }

    if (fflush(f) != 0) goto fail;
    fclose(f);
    remove(path);
    if (rename(tmp, path) != 0) { remove(tmp); return 0; }
    return 1;
fail:
    fclose(f);
    remove(tmp);
    return 0;
}

static void rebuild_links(UlamState *s) {
    s->nx[0] = 0;
    s->pv[0] = 0;
    memset(s->res_index, 0, s->nindex * sizeof(uint32_t));
    for (size_t i = 0; i < s->count; i++)
        setlinks(s, (uint32_t)(i + 1));
}

UlamState *ulam_load(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    uint32_t hdr[2];
    if (fread(hdr, 4, 2, f) != 2 || hdr[0] != MAGIC) goto fail_bare;
    if (hdr[1] != VERSION && hdr[1] != VERSION3 && hdr[1] != VERSION2)
        goto fail_bare;

    int version = (int)hdr[1];

    UlamState *s = calloc(1, sizeof *s);
    if (!s) goto fail_bare;

    if (fread(&s->max_computed, 8, 1, f) != 1) goto fail;

    uint64_t cnt;
    if (fread(&cnt, 8, 1, f) != 1) goto fail;
    s->count    = (size_t)cnt;
    s->ulam_cap = s->count > INIT_ULAM_CAP ? s->count * 2 : INIT_ULAM_CAP;

    s->a_off = malloc(s->ulam_cap * sizeof(int32_t));
    if (!s->a_off) goto fail;

    if (version == (int)VERSION) {
        if (cnt > 0 && fread(s->a_off, 4, s->count, f) != s->count) goto fail;
    } else {
        uint64_t *tmp = malloc(s->count * sizeof(uint64_t));
        if (!tmp) goto fail;
        if (cnt > 0 && fread(tmp, 8, s->count, f) != s->count) {
            free(tmp); goto fail;
        }
        for (size_t i = 0; i < s->count; i++)
            store_ulam(s, i, tmp[i]);
        free(tmp);
    }

    s->prev_ulam = (s->count > 0) ? get_ulam(s, (uint32_t)s->count) : 0;

    if (version == (int)VERSION2) {
        uint64_t sc_cap;
        if (fread(&sc_cap, 8, 1, f) != 1) goto fail;
        if (fseek(f, (long)sc_cap, SEEK_CUR) != 0) goto fail;
        s->k_cap = (size_t)(s->max_computed / 30) + 2;
        s->k     = calloc(s->k_cap, sizeof(uint32_t));
        if (!s->k) goto fail;
        for (size_t i = 0; i < s->count; i++)
            set_ulam_bit(s, get_ulam(s, (uint32_t)(i + 1)));
    } else {
        uint64_t kc;
        if (fread(&kc, 8, 1, f) != 1) goto fail;
        s->k_cap = (size_t)kc;
        s->k     = malloc(s->k_cap * sizeof(uint32_t));
        if (!s->k) goto fail;
        if (s->k_cap > 0 && fread(s->k, 4, s->k_cap, f) != s->k_cap) goto fail;
    }

    int32_t tp;
    if (fread(&tp, 4, 1, f) != 1) goto fail;
    s->track_pairs = (int)tp;
    if (tp) {
        uint64_t pc;
        if (fread(&pc, 8, 1, f) != 1) goto fail;
        s->pairs_count = (size_t)pc;
        s->pairs_cap   = s->pairs_count > INIT_PAIR_CAP
                         ? s->pairs_count * 2 : INIT_PAIR_CAP;
        s->pairs = malloc(s->pairs_cap * 2 * sizeof(uint64_t));
        if (!s->pairs) goto fail;
        if (s->pairs_count > 0 &&
            fread(s->pairs, 8, s->pairs_count * 2, f) != s->pairs_count * 2)
            goto fail;
    }

    fclose(f);

    s->nx = calloc(s->ulam_cap + 1, sizeof(uint32_t));
    s->pv = calloc(s->ulam_cap + 1, sizeof(uint32_t));
    if (!s->nx || !s->pv) { ulam_free(s); return NULL; }
    s->nindex    = s->ulam_cap / 100 + 1;
    s->res_index = calloc(s->nindex, sizeof(uint32_t));
    if (!s->res_index) { ulam_free(s); return NULL; }

    rebuild_links(s);
    return s;

fail:
    ulam_free(s);
fail_bare:
    fclose(f);
    return NULL;
}
