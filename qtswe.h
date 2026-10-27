/*
 * qtswe.h  --  QTSwE (Quantum-inspired Tabu Search with Entanglement) core
 *
 * Based on: Kuo, Chou & Chen, "Quantum-inspired algorithm for cyber-physical
 * visual surveillance deployment systems", Computer Networks 117 (2017) 5-18.
 *
 * Problem: on an R x C grid, place the FEWEST sensors such that
 *   (1) every grid cell g has detection probability Pg >= Tg      (Eq.1, Eq.3)
 *   (2) (optional) all sensors are connected to the sink (|CC| = 1) via hops <= Rt
 */
#ifndef QTSWE_H
#define QTSWE_H
#include <stdint.h>

typedef struct {
    /* map */
    int    rows, cols;
    int    sink_r, sink_c;          /* 0-indexed; -1 => map centre                */
    /* detection / network model (paper Section 3 / 6) */
    double Tg, Rcon, Rmax, alpha, beta, Rt;
    int    use_conn;                /* 1: enforce connectivity to the sink        */
    /* QTSwE parameters (paper Table 13) */
    int    pop;                     /* population p per generation       (100)    */
    int    gens;                    /* generations per attempt           (500)    */
    int    detect;                  /* stagnation iterations before LS   (10)     */
    int    local_radius;            /* neighbour reach in LS: 1 = paper's 8 nbrs  (2)     */
    int    local_rounds;            /* max local-search sweeps           (10)     */
    double delta;                   /* Q-matrix rotation step            (0.1)    */
    double w1, w2;                  /* fitness weights f=w1*f1+w2*f2     (1,1)    */
    /* ablation switches */
    int    entangle;                /* 1: entanglement-inspired initialisation    */
    int    sa_ref;                  /* A: SA acceptance of the reference solution (drives Q update) */
    int    sa_ls;                   /* B: SA acceptance of worse moves inside local search          */
    double sa_T0, sa_Tend;          /* temperature schedule per attempt (geometric): T0 -> Tend      */
    int    local_search;            /* 1: local search on stagnation              */
    /* search control */
    int    retries;                 /* extra attempts at the first failing k      (2)     */
    int    fixed_k;                 /* >0: only test this sensor count            */
    double time_limit;              /* total wall-clock seconds, 0 = unlimited    */
    uint64_t seed;                  /* 0 => derived from clock                    */
    int    include_pg;              /* include Pg matrix in JSON output           */
    int    threads;                 /* 0 = auto (OpenMP)                          */
    /* extension hooks (future ADQISA etc.).  May be NULL. */
    void  (*on_generation)(int k, int gen, double best_f, void *ud);
    /* return non-zero to stop the whole search after an attempt */
    int   (*on_attempt)(int k, int feasible, int gens_used, void *ud);
    void  *hook_ud;
} qp_t;

typedef struct {
    int    k, feasible, gens_used;
    int    f1, f2;
    double best_f, sec;
} qattempt_t;

typedef struct {
    int     k;                      /* sensors in the returned layout             */
    int     feasible;               /* 1 if all constraints verified in double    */
    int    *pos;                    /* k cell indices (r*cols+c, 0-indexed)       */
    int     f1, f2;                 /* covered cells, sensors connected to sink   */
    double  pg_min, pg_avg;
    double *pg;                     /* rows*cols, only if include_pg              */
    qattempt_t *att; int natt;      /* every (k) tried                            */
    float  *curve; int ncurve;      /* best-fitness per generation, last attempt  */
    int     k_start;                /* initial guess for k                        */
    double  s_prop;                 /* proper lattice spacing derived from Tg     */
    double  sec;                    /* total seconds                              */
    uint64_t evals;                 /* fitness evaluations                        */
    int     timed_out;
    int     rows, cols, sink_r, sink_c;
    uint64_t seed_used;
} qr_t;

void qtswe_defaults(qp_t *p);
/* returns 0 ok, <0 unknown key */
int  qtswe_set_param(qp_t *p, const char *key, const char *val);
/* returns NULL if ok else error text */
const char *qtswe_validate(const qp_t *p);
int  qtswe_optimize(const qp_t *p, qr_t *out);
void qtswe_free(qr_t *r);
int  qtswe_threads(void);

#endif
