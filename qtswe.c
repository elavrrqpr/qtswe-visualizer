/* qtswe.c -- see qtswe.h */
#define _GNU_SOURCE
#include "qtswe.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ rng */
typedef struct { uint64_t s; } rng_t;
static inline uint64_t rng_next(rng_t *r) {
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline double rng_u(rng_t *r) { return (rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }
static inline int rng_int(rng_t *r, int n) { return (int)(rng_u(r) * n); }

#ifdef _WIN32
#include <windows.h>
static double now_sec(void) {
    LARGE_INTEGER f, c; QueryPerformanceFrequency(&f); QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}
#else
static double now_sec(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
#endif

int qtswe_threads(void) {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* ------------------------------------------------------------- defaults */
void qtswe_defaults(qp_t *p) {
    memset(p, 0, sizeof *p);
    p->rows = 10; p->cols = 10; p->sink_r = -1; p->sink_c = -1;
    p->Tg = 0.3; p->Rcon = 1.414; p->Rmax = 5; p->alpha = 1; p->beta = 1; p->Rt = 5;
    p->use_conn = 1;
    p->pop = 100; p->gens = 500; p->detect = 10; p->local_rounds = 10;
    p->delta = 0.1; p->w1 = 1; p->w2 = 1;
    p->entangle = 1; p->local_search = 1; p->retries = 2; p->local_radius = 2; p->sa_T0 = 2.0; p->sa_Tend = 0.1;
    p->seed = 1;
}

int qtswe_set_param(qp_t *p, const char *k, const char *v) {
    double d = atof(v); long long i = atoll(v);
#define EQ(s) (strcmp(k, s) == 0)
    if (EQ("rows") || EQ("m")) p->rows = (int)i;
    else if (EQ("cols") || EQ("n")) p->cols = (int)i;
    else if (EQ("sink_r") || EQ("sinkr")) p->sink_r = (int)i - 1;   /* user-facing 1-indexed */
    else if (EQ("sink_c") || EQ("sinkc")) p->sink_c = (int)i - 1;
    else if (EQ("tg") || EQ("Tg")) p->Tg = d;
    else if (EQ("rcon")) p->Rcon = d;
    else if (EQ("rmax")) p->Rmax = d;
    else if (EQ("alpha")) p->alpha = d;
    else if (EQ("beta")) p->beta = d;
    else if (EQ("rt")) p->Rt = d;
    else if (EQ("conn")) p->use_conn = (int)i;
    else if (EQ("pop")) p->pop = (int)i;
    else if (EQ("gens")) p->gens = (int)i;
    else if (EQ("detect")) p->detect = (int)i;
    else if (EQ("local_rounds") || EQ("localsize")) p->local_rounds = (int)i;
    else if (EQ("delta")) p->delta = d;
    else if (EQ("w1")) p->w1 = d;
    else if (EQ("w2")) p->w2 = d;
    else if (EQ("entangle")) p->entangle = (int)i;
    else if (EQ("local")) p->local_search = (int)i;
    else if (EQ("local_radius") || EQ("lradius")) p->local_radius = (int)i;
    else if (EQ("sa_ref")) p->sa_ref = (int)i;
    else if (EQ("sa_ls")) p->sa_ls = (int)i;
    else if (EQ("sa_t0")) p->sa_T0 = d;
    else if (EQ("sa_tend")) p->sa_Tend = d;
    else if (EQ("retries")) p->retries = (int)i;
    else if (EQ("k")) p->fixed_k = (int)i;
    else if (EQ("time")) p->time_limit = d;
    else if (EQ("seed")) p->seed = (uint64_t)i;
    else if (EQ("include_pg")) p->include_pg = (int)i;
    else if (EQ("threads")) p->threads = (int)i;
    else return -1;
#undef EQ
    return 0;
}

const char *qtswe_validate(const qp_t *p) {
    if (p->rows < 1 || p->cols < 1) return "rows/cols must be >= 1";
    if ((long long)p->rows * p->cols < 2) return "map too small";
    if ((long long)p->rows * p->cols > 4000000LL) return "map too large (limit 4,000,000 cells)";
    if (p->Tg <= 0 || p->Tg >= 1) return "Tg must be in (0,1)";
    if (p->Rmax <= 0 || p->Rmax > 50) return "Rmax must be in (0,50]";
    if (p->Rcon < 0 || p->Rcon > p->Rmax) return "Rcon must be in [0,Rmax]";
    if (p->alpha <= 0 || p->alpha > 1) return "alpha must be in (0,1]";
    if (p->beta < 0) return "beta must be >= 0";
    if (p->Rt <= 0 || p->Rt > 100) return "Rt must be in (0,100]";
    if (p->pop < 2 || p->pop > 5000) return "pop must be in [2,5000]";
    if (p->gens < 1) return "gens must be >= 1";
    if (p->delta <= 0 || p->delta >= 1) return "delta must be in (0,1)";
    if (p->sa_T0 <= 0 || p->sa_Tend <= 0 || p->sa_Tend > p->sa_T0) return "need sa_t0 >= sa_tend > 0";
    if (p->sink_r >= p->rows || p->sink_c >= p->cols) return "sink outside map";
    if (p->fixed_k >= p->rows * p->cols) return "k must be < rows*cols";
    return NULL;
}

/* -------------------------------------------------------------- context */
typedef struct { float key; int g; } pair_t;
typedef struct {
    double  *L;      /* sum of log(1-ps) per cell      */
    uint8_t *occ, *vis;
    int     *q, *tmp;
    pair_t  *pairs;
} scratch_t;

typedef struct {
    const qp_t *P;
    int R, C, RC, sinkcell;
    int nfp; int *fdr, *fdc; double *fw, *fps;  /* detection footprint            */
    int nrt; int *tdr, *tdc;                    /* connectivity offsets (<=Rt)    */
    double thr;
    rng_t rng;
    scratch_t *scr; int nscr;
    uint64_t evals;
} ctx_t;

static void scratch_init(ctx_t *x, scratch_t *s) {
    s->L = malloc(sizeof(double) * x->RC);
    s->occ = calloc(x->RC, 1); s->vis = calloc(x->RC, 1);
    s->q = malloc(sizeof(int) * x->RC); s->tmp = malloc(sizeof(int) * x->RC);
    s->pairs = malloc(sizeof(pair_t) * x->RC);
}
static void scratch_free(scratch_t *s) {
    free(s->L); free(s->occ); free(s->vis); free(s->q); free(s->tmp); free(s->pairs);
}

static void ctx_init(ctx_t *x, const qp_t *P) {
    memset(x, 0, sizeof *x);
    x->P = P; x->R = P->rows; x->C = P->cols; x->RC = x->R * x->C;
    int sr = P->sink_r >= 0 ? P->sink_r : (x->R - 1) / 2;
    int sc = P->sink_c >= 0 ? P->sink_c : (x->C - 1) / 2;
    x->sinkcell = sr * x->C + sc;
    int m = (int)ceil(P->Rmax), cap = (2 * m + 1) * (2 * m + 1);
    x->fdr = malloc(sizeof(int) * cap); x->fdc = malloc(sizeof(int) * cap);
    x->fw = malloc(sizeof(double) * cap); x->fps = malloc(sizeof(double) * cap);
    for (int dr = -m; dr <= m; dr++) for (int dc = -m; dc <= m; dc++) {
        double d = sqrt((double)(dr * dr + dc * dc)), ps;
        if (d <= P->Rcon) ps = P->alpha;
        else if (d <= P->Rmax + 1e-9) ps = pow(d, -P->beta);
        else ps = 0;
        if (ps <= 0) continue;
        if (ps > 1) ps = 1;
        x->fdr[x->nfp] = dr; x->fdc[x->nfp] = dc; x->fps[x->nfp] = ps;
        x->fw[x->nfp] = log(fmax(1.0 - ps, 1e-12));
        x->nfp++;
    }
    int mt = (int)ceil(P->Rt);
    cap = (2 * mt + 1) * (2 * mt + 1);
    x->tdr = malloc(sizeof(int) * cap); x->tdc = malloc(sizeof(int) * cap);
    for (int dr = -mt; dr <= mt; dr++) for (int dc = -mt; dc <= mt; dc++) {
        if (!dr && !dc) continue;
        if (sqrt((double)(dr * dr + dc * dc)) <= P->Rt + 1e-9) { x->tdr[x->nrt] = dr; x->tdc[x->nrt] = dc; x->nrt++; }
    }
    x->thr = log(1.0 - P->Tg) + 1e-9;
    uint64_t seed = P->seed ? P->seed : (uint64_t)time(NULL) * 2654435761ULL;
    x->rng.s = seed;
    x->nscr = qtswe_threads();
    x->scr = calloc(x->nscr, sizeof(scratch_t));
    for (int i = 0; i < x->nscr; i++) scratch_init(x, &x->scr[i]);
}
static void ctx_free(ctx_t *x) {
    free(x->fdr); free(x->fdc); free(x->fw); free(x->fps); free(x->tdr); free(x->tdc);
    for (int i = 0; i < x->nscr; i++) scratch_free(&x->scr[i]);
    free(x->scr);
}

/* ------------------------------------------------------------- fitness */
static void compute_L(const ctx_t *x, double *L, const int *pos, int k) {
    memset(L, 0, sizeof(double) * x->RC);
    for (int i = 0; i < k; i++) {
        int r = pos[i] / x->C, c = pos[i] % x->C;
        for (int t = 0; t < x->nfp; t++) {
            int r2 = r + x->fdr[t], c2 = c + x->fdc[t];
            if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
            L[r2 * x->C + c2] += x->fw[t];
        }
    }
}
static int count_cov(const ctx_t *x, const double *L) {
    int n = 0; for (int g = 0; g < x->RC; g++) n += (L[g] <= x->thr); return n;
}
/* #sensors reachable from the sink (hop <= Rt); occ must be set */
static int conn_count(const ctx_t *x, scratch_t *s) {
    int *q = s->q, head = 0, tail = 0, cnt = 0;
    q[tail++] = x->sinkcell; s->vis[x->sinkcell] = 1;
    while (head < tail) {
        int g = q[head++], r = g / x->C, c = g % x->C;
        for (int t = 0; t < x->nrt; t++) {
            int r2 = r + x->tdr[t], c2 = c + x->tdc[t];
            if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
            int h = r2 * x->C + c2;
            if (s->occ[h] && !s->vis[h]) { s->vis[h] = 1; q[tail++] = h; cnt++; }
        }
    }
    for (int i = 0; i < tail; i++) s->vis[q[i]] = 0;
    return cnt;
}
static void set_occ(const ctx_t *x, scratch_t *s, const int *pos, int k, int v) {
    for (int i = 0; i < k; i++) s->occ[pos[i]] = (uint8_t)v;
    (void)x;
}
/* full evaluation (thread-safe: only touches the scratch) */
static void eval_sol(const ctx_t *x, scratch_t *s, const int *pos, int k, int *f1, int *f2) {
    compute_L(x, s->L, pos, k);
    *f1 = count_cov(x, s->L);
    if (x->P->use_conn) { set_occ(x, s, pos, k, 1); *f2 = conn_count(x, s); set_occ(x, s, pos, k, 0); }
    else *f2 = k;
}
static inline double fit(const ctx_t *x, int f1, int f2) { return x->P->w1 * f1 + x->P->w2 * f2; }
static inline int is_feasible(const ctx_t *x, int f1, int f2, int k) { return f1 == x->RC && f2 == k; }

/* ------------------------------------------------- Measure (quantum->bits) */
static int cmp_desc(const void *a, const void *b) {
    float x = ((const pair_t *)a)->key, y = ((const pair_t *)b)->key;
    return (x < y) - (x > y);
}
static void fill_random(const ctx_t *x, scratch_t *s, rng_t *r, int *pos, int m, int k) {
    set_occ(x, s, pos, m, 1);
    while (m < k) {
        int g = rng_int(r, x->RC);
        if (g == x->sinkcell || s->occ[g]) continue;
        s->occ[g] = 1; pos[m++] = g;
    }
    set_occ(x, s, pos, k, 0);
}
/* Measure(x): R~U[0,1]; x=1 if R > |beta|^2   (Algorithm 1).  Then repair to k. */
static void measure(const ctx_t *x, scratch_t *s, rng_t *r, const float *b2, int k, int *pos) {
    pair_t *pr = s->pairs; int n = 0;
    for (int g = 0; g < x->RC; g++) {
        if (g == x->sinkcell) continue;
        if (rng_u(r) > b2[g]) { pr[n].key = (float)((1.0 - b2[g]) * (0.5 + rng_u(r))); pr[n].g = g; n++; }
    }
    if (n > k) qsort(pr, n, sizeof(pair_t), cmp_desc);   /* keep best-supported k cells */
    int m = n > k ? k : n;
    for (int i = 0; i < m; i++) pos[i] = pr[i].g;
    if (m < k) fill_random(x, s, r, pos, m, k);
}

/* ------------------------------------------- entanglement-inspired init */
/* Sensors are laid along a chain on the x- (or y-) axis at the "proper"
 * distance, successive chains are staggered (triangular lattice) and every
 * step is randomly perturbed, like chain-mail. */
static void entangled(const ctx_t *x, scratch_t *s, rng_t *r, int k, int *pos) {
    int RC = x->RC;
    int horiz = rng_u(r) < 0.5;
    double step = sqrt((double)RC / (0.866025 * k)); if (step < 1) step = 1;
    double rs = 0.866025 * step;
    int U = horiz ? x->C : x->R, V = horiz ? x->R : x->C;
    double ox = rng_u(r) * step, oy = rng_u(r) * rs;
    int n = 0, *tmp = s->tmp, row = 0;
    memset(s->occ, 0, RC);
    for (double v = oy - rs; v < V; v += rs * (1 + 0.3 * (rng_u(r) - 0.5)), row++) {
        double sh = (row & 1) ? step / 2 : 0;
        for (double u = ox + sh - step; u < U; u += step * (1 + 0.3 * (rng_u(r) - 0.5))) {
            int uu = (int)floor(u + 0.5), vv = (int)floor(v + 0.5);
            if (uu < 0 || uu >= U || vv < 0 || vv >= V) continue;
            int g = horiz ? vv * x->C + uu : uu * x->C + vv;
            if (g == x->sinkcell || s->occ[g]) continue;
            s->occ[g] = 1; tmp[n++] = g;
        }
    }
    memset(s->occ, 0, RC);
    if (n > k) {
        for (int i = 0; i < k; i++) { int j = i + rng_int(r, n - i); int t = tmp[i]; tmp[i] = tmp[j]; tmp[j] = t; }
        n = k;
    }
    memcpy(pos, tmp, sizeof(int) * n);
    if (n < k) fill_random(x, s, r, pos, n, k);
}

/* -------------------------------------------------------- local search */
static int fp_apply(const ctx_t *x, double *L, int cell, int sgn) {
    int r = cell / x->C, c = cell % x->C, d = 0;
    for (int t = 0; t < x->nfp; t++) {
        int r2 = r + x->fdr[t], c2 = c + x->fdc[t];
        if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
        double *p = &L[r2 * x->C + c2];
        int was = *p <= x->thr; *p += sgn * x->fw[t]; d += (*p <= x->thr) - was;
    }
    return d;
}
/* For each sensor (random order) try its neighbours and take the best improving move
 * (equal-fitness moves are taken as sideways steps).  If T > 0 (SA option B) and every
 * neighbour is worse, one random neighbour is accepted with prob exp(-delta/T).
 * Repeats up to local_rounds sweeps.  The best layout seen along the way is returned. */
static void local_search(ctx_t *x, scratch_t *s, rng_t *rg, int *pos, int k, int *of1, int *of2, double T) {
    const qp_t *P = x->P;
    compute_L(x, s->L, pos, k);
    int f1 = count_cov(x, s->L), f2;
    set_occ(x, s, pos, k, 1);
    f2 = P->use_conn ? conn_count(x, s) : k;
    double cur = fit(x, f1, f2);
    int *ord = malloc(sizeof(int) * k), *bpos = malloc(sizeof(int) * k);
    for (int i = 0; i < k; i++) ord[i] = i;
    memcpy(bpos, pos, sizeof(int) * k);
    double bcur = cur; int bf1_ = f1, bf2_ = f2;
    for (int round = 0; round < P->local_rounds; round++) {
        int improved = 0;
        for (int i = k - 1; i > 0; i--) { int j = rng_int(rg, i + 1); int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
        for (int oi = 0; oi < k; oi++) {
            int i = ord[oi], a = pos[i], ar = a / x->C, ac = a % x->C;
            double bestf = cur - 1e-9; int bestb = -1, bf1 = 0, bf2 = 0, nside = 0;
            int anyb = -1, any1 = 0, any2 = 0, nany = 0; double anyf = 0;
            for (int dr = -P->local_radius; dr <= P->local_radius; dr++)
                for (int dc = -P->local_radius; dc <= P->local_radius; dc++) {
                if (!dr && !dc) continue;
                int r2 = ar + dr, c2 = ac + dc;
                if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
                int b = r2 * x->C + c2;
                if (b == x->sinkcell || s->occ[b]) continue;
                int d = fp_apply(x, s->L, a, -1) + fp_apply(x, s->L, b, +1);
                int n1 = f1 + d, n2 = k;
                if (P->use_conn) {
                    s->occ[a] = 0; s->occ[b] = 1;
                    n2 = conn_count(x, s);
                    s->occ[b] = 0; s->occ[a] = 1;
                }
                x->evals++;
                double nf = fit(x, n1, n2);
                if (nf > bestf + 1e-9) { bestf = nf; bestb = b; bf1 = n1; bf2 = n2; nside = 0; }
                else if (nf > cur - 1e-9 && bestf <= cur + 1e-9) {          /* equal fitness: sideways move (reservoir pick) */
                    nside++;
                    if (rng_int(rg, nside) == 0) { bestf = cur; bestb = b; bf1 = n1; bf2 = n2; }
                }
                if (T > 0) { nany++; if (rng_int(rg, nany) == 0) { anyb = b; anyf = nf; any1 = n1; any2 = n2; } }
                fp_apply(x, s->L, b, -1); fp_apply(x, s->L, a, +1);
            }
            if (bestb < 0 && T > 0 && anyb >= 0) {                          /* stuck: SA may accept a worse move */
                double dl = cur - anyf;
                if (dl > 0 && rng_u(rg) < exp(-dl / T)) { bestb = anyb; bestf = anyf; bf1 = any1; bf2 = any2; }
            }
            if (bestb >= 0) {
                fp_apply(x, s->L, a, -1); fp_apply(x, s->L, bestb, +1);
                s->occ[a] = 0; s->occ[bestb] = 1; pos[i] = bestb;
                if (bestf > cur + 1e-9) improved = 1;
                cur = fit(x, bf1, bf2); f1 = bf1; f2 = bf2;
                if (cur > bcur + 1e-9) { bcur = cur; memcpy(bpos, pos, sizeof(int) * k); bf1_ = f1; bf2_ = f2; }
                if (is_feasible(x, f1, f2, k)) goto done;
            }
        }
        if (!improved && T <= 0 && round >= 2) break;
    }
done:
    set_occ(x, s, pos, k, 0);
    if (bcur > cur + 1e-9) { memcpy(pos, bpos, sizeof(int) * k); f1 = bf1_; f2 = bf2_; }
    free(ord); free(bpos);
    *of1 = f1; *of2 = f2;
}

/* --------------------------------------- adapt a solution to a new k */
static void adapt(ctx_t *x, scratch_t *s, rng_t *r, const int *from, int kf, int *out, int kt) {
    memcpy(out, from, sizeof(int) * (kf < kt ? kf : kt));
    int *tmp = malloc(sizeof(int) * (kf > kt ? kf : kt));
    memcpy(tmp, from, sizeof(int) * kf);
    int k = kf;
    compute_L(x, s->L, tmp, k);
    if (kt < kf) {
        while (k > kt) {
            int bi = -1, bl = 1 << 30;
            for (int i = 0; i < k; i++) {
                int rr = tmp[i] / x->C, cc = tmp[i] % x->C, loss = 0;
                for (int t = 0; t < x->nfp; t++) {
                    int r2 = rr + x->fdr[t], c2 = cc + x->fdc[t];
                    if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
                    double l = s->L[r2 * x->C + c2];
                    if (l <= x->thr && l - x->fw[t] > x->thr) loss++;
                }
                if (loss < bl || (loss == bl && rng_u(r) < 0.3)) { bl = loss; bi = i; }
            }
            fp_apply(x, s->L, tmp[bi], -1);
            tmp[bi] = tmp[--k];
        }
    } else if (kt > kf) {
        set_occ(x, s, tmp, k, 1);
        while (k < kt) {
            int bestg = -1, bestgain = -1;
            for (int tr = 0; tr < 12; tr++) {
                int g = -1;
                for (int a = 0; a < 40; a++) {            /* random uncovered cell */
                    int c = rng_int(r, x->RC);
                    if (s->L[c] > x->thr) { g = c; break; }
                }
                if (g < 0) g = rng_int(r, x->RC);
                if (g == x->sinkcell || s->occ[g]) continue;
                int rr = g / x->C, cc = g % x->C, gain = 0;
                for (int t = 0; t < x->nfp; t++) {
                    int r2 = rr + x->fdr[t], c2 = cc + x->fdc[t];
                    if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
                    double l = s->L[r2 * x->C + c2];
                    if (l > x->thr && l + x->fw[t] <= x->thr) gain++;
                }
                if (gain > bestgain) { bestgain = gain; bestg = g; }
            }
            if (bestg < 0) continue;
            fp_apply(x, s->L, bestg, +1); s->occ[bestg] = 1; tmp[k++] = bestg;
        }
        set_occ(x, s, tmp, k, 0);
    }
    memcpy(out, tmp, sizeof(int) * kt);
    free(tmp);
}

/* ------------------------------------------------------------ one attempt */
typedef struct { int feasible, gens, f1, f2; double f; } ares_t;

static ares_t attempt(ctx_t *x, int k, const int *seed_pos, int *best_pos,
                      float **curve, int *ncurve, double deadline) {
    const qp_t *P = x->P; int p = P->pop, RC = x->RC;
    ares_t res = {0, 0, 0, 0, -1};
    float *b2 = malloc(sizeof(float) * RC);
    for (int g = 0; g < RC; g++) b2[g] = 1.0f - (float)k / RC;
    if (seed_pos) for (int i = 0; i < k; i++) b2[seed_pos[i]] = 0.5f;   /* warm start: bias Q towards the known layout */
    int *popPos = malloc(sizeof(int) * (size_t)p * k), *ref_pos = malloc(sizeof(int) * k);
    double *pf = malloc(sizeof(double) * p); int *pf1 = malloc(sizeof(int) * p), *pf2 = malloc(sizeof(int) * p);
    uint64_t *seeds = malloc(sizeof(uint64_t) * p);
    uint8_t *bm = malloc(RC), *wm = malloc(RC);
    float *cv = malloc(sizeof(float) * (P->gens + 1)); int nc = 0;
    /* reference solution: attracts Q (may get worse under SA-A).
     * history best (best_pos): best layout ever seen; this is what is returned. */
    double rbf = -1, hbf = -1; int rb1 = 0, rb2 = 0, hb1 = 0, hb2 = 0;
    int stall = seed_pos ? P->detect : 0;
    int use_sa = P->sa_ref || P->sa_ls;
    scratch_t *s0 = &x->scr[0];

    for (int gen = 0; gen < P->gens; gen++) {
        double T = use_sa ? P->sa_T0 * pow(P->sa_Tend / P->sa_T0, P->gens > 1 ? (double)gen / (P->gens - 1) : 1.0) : 0;
        for (int i = 0; i < p; i++) seeds[i] = rng_next(&x->rng);
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < p; i++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            scratch_t *s = &x->scr[tid]; rng_t r = { seeds[i] };
            int *pos = popPos + (size_t)i * k;
            if (gen == 0 && i == 0 && seed_pos) memcpy(pos, seed_pos, sizeof(int) * k);
            else if (gen == 0 && P->entangle) entangled(x, s, &r, k, pos);
            else measure(x, s, &r, b2, k, pos);
            eval_sol(x, s, pos, k, &pf1[i], &pf2[i]);
            pf[i] = fit(x, pf1[i], pf2[i]);
        }
        x->evals += p;
        int bi = 0, wi = 0;
        for (int i = 1; i < p; i++) { if (pf[i] > pf[bi]) bi = i; if (pf[i] < pf[wi]) wi = i; }

        int hist_up = 0;                                   /* history best */
        if (pf[bi] > hbf) {
            hbf = pf[bi]; hb1 = pf1[bi]; hb2 = pf2[bi];
            memcpy(best_pos, popPos + (size_t)bi * k, sizeof(int) * k); hist_up = 1;
        }
        int acc = 0;                                       /* reference solution */
        if (!P->sa_ref) acc = (pf[bi] > rbf);              /* greedy (paper) */
        else {                                             /* SA-A: delta = ref - candidate (maximisation) */
            double dl = rbf - pf[bi];
            acc = (rbf < 0) || dl < 0 || (T > 0 ? rng_u(&x->rng) < exp(-dl / T) : dl <= 0);
        }
        if (acc) { rbf = pf[bi]; rb1 = pf1[bi]; rb2 = pf2[bi]; memcpy(ref_pos, popPos + (size_t)bi * k, sizeof(int) * k); }
        stall = hist_up ? 0 : stall + 1;
        cv[nc++] = (float)hbf; res.gens = gen + 1;
        if (P->on_generation) P->on_generation(k, gen, hbf, P->hook_ud);
        if (is_feasible(x, hb1, hb2, k)) break;

        if (P->local_search && stall >= P->detect) {       /* stagnation -> local search on the reference */
            rng_t lr = { rng_next(&x->rng) };
            local_search(x, s0, &lr, ref_pos, k, &rb1, &rb2, P->sa_ls ? T : 0);
            rbf = fit(x, rb1, rb2); stall = 0;
            if (rbf > hbf || (!P->sa_ref && rbf >= hbf)) { hbf = rbf; hb1 = rb1; hb2 = rb2; memcpy(best_pos, ref_pos, sizeof(int) * k); }
            if (nc) cv[nc - 1] = (float)hbf;
            if (is_feasible(x, hb1, hb2, k)) break;
        }
        /* update Q: attract to the reference, repel from this generation's worst */
        memset(bm, 0, RC); memset(wm, 0, RC);
        for (int i = 0; i < k; i++) { bm[ref_pos[i]] = 1; wm[popPos[(size_t)wi * k + i]] = 1; }
        float dl2 = (float)P->delta;
        for (int g = 0; g < RC; g++) {
            if (bm[g] == wm[g]) continue;
            float v = b2[g] + (bm[g] ? -dl2 : dl2);
            b2[g] = v < 0.001f ? 0.001f : (v > 0.999f ? 0.999f : v);
        }
        if (deadline > 0 && now_sec() > deadline) break;
    }
    res.f = hbf; res.f1 = hb1; res.f2 = hb2; res.feasible = is_feasible(x, hb1, hb2, k);
    *curve = cv; *ncurve = nc;
    free(b2); free(popPos); free(ref_pos); free(pf); free(pf1); free(pf2); free(seeds); free(bm); free(wm);
    return res;
}

/* ---------------------------------------------- proper spacing from Tg */
static double proper_spacing(const qp_t *P) {
    /* triangular lattice; worst point = triangle centroid */
    for (double s = 12.0; s >= 0.5; s -= 0.05) {
        double miss = 1.0, px = s / 2, py = 0.288675 * s;
        int m = (int)(P->Rmax / (0.866 * s)) + 3, mi = (int)(P->Rmax / s) + 4;
        for (int j = -m; j <= m; j++) for (int i = -mi - m; i <= mi + m; i++) {
            double lx = i * s + j * s / 2, ly = j * s * 0.866025;
            double d = hypot(px - lx, py - ly), ps;
            if (d <= P->Rcon) ps = P->alpha; else if (d <= P->Rmax) ps = pow(d, -P->beta); else ps = 0;
            if (ps > 1) ps = 1;
            miss *= (1 - ps);
        }
        if (1 - miss >= P->Tg) return s;
    }
    return 0.5;
}

/* ----------------------------------------------------- final verification */
static void finalize(ctx_t *x, qr_t *r, const int *pos, int k) {
    const qp_t *P = x->P;
    r->k = k; r->pos = malloc(sizeof(int) * k); memcpy(r->pos, pos, sizeof(int) * k);
    double *miss = malloc(sizeof(double) * x->RC);
    for (int g = 0; g < x->RC; g++) miss[g] = 1.0;
    for (int i = 0; i < k; i++) {
        int rr = pos[i] / x->C, cc = pos[i] % x->C;
        for (int t = 0; t < x->nfp; t++) {
            int r2 = rr + x->fdr[t], c2 = cc + x->fdc[t];
            if ((unsigned)r2 >= (unsigned)x->R || (unsigned)c2 >= (unsigned)x->C) continue;
            miss[r2 * x->C + c2] *= (1 - x->fps[t]);
        }
    }
    int cov = 0; double mn = 1, sum = 0;
    for (int g = 0; g < x->RC; g++) {
        double pg = 1 - miss[g];
        if (miss[g] <= (1 - P->Tg) * (1 + 1e-9)) cov++;
        if (pg < mn) mn = pg;
        sum += pg;
        miss[g] = pg;
    }
    scratch_t *s = &x->scr[0];
    set_occ(x, s, pos, k, 1);
    r->f2 = P->use_conn ? conn_count(x, s) : k;
    set_occ(x, s, pos, k, 0);
    r->f1 = cov; r->pg_min = mn; r->pg_avg = sum / x->RC;
    r->feasible = (cov == x->RC && r->f2 == k);
    if (P->include_pg) r->pg = miss; else free(miss);
}

static void add_attempt(qr_t *r, int *cap, int k, ares_t a, double sec) {
    if (r->natt == *cap) { *cap = *cap ? *cap * 2 : 16; r->att = realloc(r->att, sizeof(qattempt_t) * *cap); }
    r->att[r->natt++] = (qattempt_t){ k, a.feasible, a.gens, a.f1, a.f2, a.f, sec };
}

/* ------------------------------------------------------------ main driver */
int qtswe_optimize(const qp_t *P, qr_t *out) {
    memset(out, 0, sizeof *out);
    if (qtswe_validate(P)) return -1;
#ifdef _OPENMP
    if (P->threads > 0) omp_set_num_threads(P->threads);
#endif
    ctx_t x; ctx_init(&x, P);
    double t0 = now_sec(), deadline = P->time_limit > 0 ? t0 + P->time_limit : 0;
    int RC = x.RC, kmax = RC - 1;
    out->rows = x.R; out->cols = x.C; out->sink_r = x.sinkcell / x.C; out->sink_c = x.sinkcell % x.C;
    out->seed_used = x.rng.s;
    out->s_prop = proper_spacing(P);
    double sp = out->s_prop;
    int k0 = (int)ceil((x.R + 0.5 * sp) * (x.C + 0.5 * sp) / (0.866025 * sp * sp));
    if (P->fixed_k > 0) k0 = P->fixed_k;
    if (k0 < 1) k0 = 1;
    if (k0 > kmax) k0 = kmax;
    out->k_start = k0;

    int cap = 0, hi = -1, k = k0, step = k0 / 8 > 0 ? k0 / 8 : 1, fails = 0;
    int *hi_pos = NULL, *last_pos = NULL, last_k = 0, *cur = NULL;
    int *seedbuf = malloc(sizeof(int) * (kmax + 1));
    float *curve = NULL; int ncurve = 0;
    /* Phase 1: gallop upwards from the lattice estimate until a feasible layout is found.
     * Phase 2: remove one sensor at a time (warm start from the best feasible layout)
     *          until the QTSwE search can no longer repair the layout. */
    for (;;) {
        cur = malloc(sizeof(int) * k);
        const int *seed = NULL;
        if (hi_pos && k < hi) { adapt(&x, &x.scr[0], &x.rng, hi_pos, hi, seedbuf, k); seed = seedbuf; }
        else if (last_pos && k > last_k) { adapt(&x, &x.scr[0], &x.rng, last_pos, last_k, seedbuf, k); seed = seedbuf; }
        double ta = now_sec();
        free(curve); curve = NULL;
        ares_t a = attempt(&x, k, seed, cur, &curve, &ncurve, deadline);
        add_attempt(out, &cap, k, a, now_sec() - ta);
        int stop_user = P->on_attempt ? P->on_attempt(k, a.feasible, a.gens, P->hook_ud) : 0;
        if (a.feasible) { hi = k; free(hi_pos); hi_pos = cur; fails = 0; }
        else {
            free(last_pos); last_pos = cur; last_k = k;
            if (hi > 0) fails++;
        }
        if (P->fixed_k > 0 || stop_user) break;
        if (deadline > 0 && now_sec() > deadline) { out->timed_out = 1; break; }
        if (hi < 0) { if (k >= kmax) break; k += step; if (k > kmax) k = kmax; step *= 2; }
        else if (a.feasible) { if (hi <= 1) break; k = hi - 1; }
        else if (fails > P->retries) break;
        else k = hi - 1;                      /* retry the same k with a new random stream */
    }
    if (hi > 0) finalize(&x, out, hi_pos, hi);
    else finalize(&x, out, last_pos, last_k);
    out->curve = curve; out->ncurve = ncurve;
    out->evals = x.evals; out->sec = now_sec() - t0;
    free(hi_pos); free(last_pos); free(seedbuf);
    ctx_free(&x);
    return 0;
}

void qtswe_free(qr_t *r) { free(r->pos); free(r->pg); free(r->att); free(r->curve); memset(r, 0, sizeof *r); }
