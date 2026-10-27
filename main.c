/* main.c -- CLI + tiny HTTP server for the QTSwE backend
 *
 *   qtswe serve [--port 8080] [--host 127.0.0.1]
 *   qtswe run   --rows 25 --cols 25 --tg 0.3 [--json] [--seed N] ...
 *   qtswe bench --rows 25 --cols 25 --tg 0.3 --runs 10   (ablation table, CSV)
 */
#define _GNU_SOURCE
#include "qtswe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/stat.h>
#ifdef _WIN32
  #include <direct.h>
  #define MKDIR(d) _mkdir(d)
#else
  #define MKDIR(d) mkdir((d), 0755)
#endif
#ifdef _WIN32
  #include <winsock2.h>
  #include <ws2tcpip.h>
  typedef SOCKET sock_t;
  #define sock_close closesocket
  #define strtok_r strtok_s
  static long sock_read(sock_t s, char *b, size_t n)        { return recv(s, b, (int)n, 0); }
  static long sock_write(sock_t s, const char *b, size_t n) { return send(s, b, (int)n, 0); }
#else
  #include <unistd.h>
  #include <signal.h>
  #include <sys/socket.h>
  #include <netinet/in.h>
  #include <arpa/inet.h>
  typedef int sock_t;
  #define sock_close close
  static long sock_read(sock_t s, char *b, size_t n)        { return read(s, b, n); }
  static long sock_write(sock_t s, const char *b, size_t n) { return write(s, b, n); }
#endif

/* ---------------------------------------------------------- string buffer */
typedef struct { char *s; size_t n, cap; } sb_t;
static void sb_put(sb_t *b, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
#include <stdarg.h>
static void sb_put(sb_t *b, const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    for (;;) {
        va_list cp; va_copy(cp, ap);
        int n = vsnprintf(b->s + b->n, b->cap - b->n, fmt, cp); va_end(cp);
        if (n >= 0 && (size_t)n < b->cap - b->n) { b->n += n; break; }
        b->cap = b->cap ? b->cap * 2 : 4096;
        if (n > 0 && b->n + n + 1 > b->cap) b->cap = b->n + n + 1;
        b->s = realloc(b->s, b->cap);
    }
    va_end(ap);
}

static void result_json(sb_t *b, const qp_t *P, const qr_t *r) {
    sb_put(b, "{\"ok\":true,\"rows\":%d,\"cols\":%d,\"sink\":[%d,%d],", r->rows, r->cols, r->sink_r + 1, r->sink_c + 1);
    sb_put(b, "\"params\":{\"Tg\":%g,\"Rcon\":%g,\"Rmax\":%g,\"alpha\":%g,\"beta\":%g,\"Rt\":%g,\"conn\":%d,"
              "\"pop\":%d,\"gens\":%d,\"detect\":%d,\"local_rounds\":%d,\"delta\":%g,\"entangle\":%d,\"local\":%d,\"sa_ref\":%d,\"sa_ls\":%d,\"sa_T0\":%g,\"sa_Tend\":%g,\"seed\":%llu},",
           P->Tg, P->Rcon, P->Rmax, P->alpha, P->beta, P->Rt, P->use_conn, P->pop, P->gens, P->detect,
           P->local_rounds, P->delta, P->entangle, P->local_search, P->sa_ref, P->sa_ls, P->sa_T0, P->sa_Tend, (unsigned long long)r->seed_used);
    sb_put(b, "\"feasible\":%s,\"timed_out\":%s,\"k\":%d,\"covered\":%d,\"cells\":%d,\"connected\":%d,"
              "\"pg_min\":%.6f,\"pg_avg\":%.6f,\"k_start\":%d,\"s_prop\":%.3f,\"seconds\":%.3f,\"evals\":%llu,\"threads\":%d,",
           r->feasible ? "true" : "false", r->timed_out ? "true" : "false", r->k, r->f1, r->rows * r->cols, r->f2,
           r->pg_min, r->pg_avg, r->k_start, r->s_prop, r->sec, (unsigned long long)r->evals, qtswe_threads());
    sb_put(b, "\"sensors\":[");
    for (int i = 0; i < r->k; i++) sb_put(b, "%s[%d,%d]", i ? "," : "", r->pos[i] / r->cols + 1, r->pos[i] % r->cols + 1);
    sb_put(b, "],\"attempts\":[");
    for (int i = 0; i < r->natt; i++) {
        const qattempt_t *a = &r->att[i];
        sb_put(b, "%s{\"k\":%d,\"feasible\":%s,\"gens\":%d,\"covered\":%d,\"connected\":%d,\"fitness\":%g,\"seconds\":%.3f}",
               i ? "," : "", a->k, a->feasible ? "true" : "false", a->gens_used, a->f1, a->f2, a->best_f, a->sec);
    }
    sb_put(b, "],\"curve\":[");
    int stride = r->ncurve > 600 ? (r->ncurve + 599) / 600 : 1, first = 1;
    for (int i = 0; i < r->ncurve; i += stride) { sb_put(b, "%s%g", first ? "" : ",", r->curve[i]); first = 0; }
    sb_put(b, "]");
    if (r->pg) {
        sb_put(b, ",\"pg\":[");
        for (int i = 0; i < r->rows * r->cols; i++) sb_put(b, "%s%.4f", i ? "," : "", r->pg[i]);
        sb_put(b, "]");
    }
    sb_put(b, "}");
}

static void ascii_map(const qp_t *P, const qr_t *r) {
    (void)P;
    if (r->cols > 120 || r->rows > 80) { printf("(map too large for ASCII view)\n"); return; }
    char *g = malloc((size_t)r->rows * r->cols);
    memset(g, '.', (size_t)r->rows * r->cols);
    for (int i = 0; i < r->k; i++) g[r->pos[i]] = 'o';
    g[r->sink_r * r->cols + r->sink_c] = 'S';
    for (int y = 0; y < r->rows; y++) { fwrite(g + (size_t)y * r->cols, 1, r->cols, stdout); putchar('\n'); }
    free(g);
}

static int parse_args(qp_t *P, int argc, char **argv, int *json) {
    for (int i = 0; i < argc; i++) {
        if (strncmp(argv[i], "--", 2)) { fprintf(stderr, "bad arg %s\n", argv[i]); return -1; }
        const char *k = argv[i] + 2;
        if (!strcmp(k, "json")) { *json = 1; continue; }
        if (i + 1 >= argc) { fprintf(stderr, "missing value for %s\n", argv[i]); return -1; }
        if (qtswe_set_param(P, k, argv[++i])) { fprintf(stderr, "unknown option --%s\n", k); return -1; }
    }
    return 0;
}

static int cmd_run(int argc, char **argv) {
    qp_t P; qtswe_defaults(&P); int json = 0;
    if (parse_args(&P, argc, argv, &json)) return 2;
    const char *e = qtswe_validate(&P); if (e) { fprintf(stderr, "error: %s\n", e); return 2; }
    qr_t r; if (qtswe_optimize(&P, &r)) return 1;
    if (json) { sb_t b = {0}; result_json(&b, &P, &r); puts(b.s); free(b.s); }
    else {
        printf("map %dx%d  Tg=%.2f  sink=(%d,%d)  seed=%llu  threads=%d\n", r.rows, r.cols, P.Tg, r.sink_r + 1, r.sink_c + 1,
               (unsigned long long)r.seed_used, qtswe_threads());
        for (int i = 0; i < r.natt; i++)
            printf("  try k=%-5d %-10s gens=%-4d covered=%d/%d conn=%d  %.2fs\n", r.att[i].k, r.att[i].feasible ? "FEASIBLE" : "infeasible",
                   r.att[i].gens_used, r.att[i].f1, r.rows * r.cols, r.att[i].f2, r.att[i].sec);
        printf("RESULT: %s  sensors=%d  covered=%d/%d  connected=%d/%d  Pg_min=%.4f  %.2fs\n",
               r.feasible ? "feasible" : "NOT feasible", r.k, r.f1, r.rows * r.cols, r.f2, r.k, r.pg_min, r.sec);
        ascii_map(&P, &r);
    }
    qtswe_free(&r); return 0;
}

/* ------------------------------------------------------------- ablation */
static int cmd_bench(int argc, char **argv) {
    qp_t P; qtswe_defaults(&P); int json = 0, runs = 5;
    /* extract --runs */
    char **a2 = malloc(sizeof(char *) * (argc + 1)); int n2 = 0;
    for (int i = 0; i < argc; i++) { if (!strcmp(argv[i], "--runs") && i + 1 < argc) runs = atoi(argv[++i]); else a2[n2++] = argv[i]; }
    uint64_t base = 1;
    if (parse_args(&P, n2, a2, &json)) return 2;
    if (P.seed) base = P.seed;
    const char *e = qtswe_validate(&P); if (e) { fprintf(stderr, "error: %s\n", e); return 2; }
    static const struct { const char *name; int ent, loc, sa, sb; } cf[] = {
        {"full", 1, 1, 0, 0}, {"no_entangle", 0, 1, 0, 0}, {"no_local", 1, 0, 0, 0}, {"plain_QTS", 0, 0, 0, 0},
        {"full+SA_A", 1, 1, 1, 0}, {"full+SA_B", 1, 1, 0, 1}, {"full+SA_AB", 1, 1, 1, 1} };
    const int NC = (int)(sizeof cf / sizeof cf[0]);
    printf("config,seed,rows,cols,Tg,k,feasible,seconds,evals,attempts,gens_last,covered_last\n");
    double sumk[16] = {0}, sumt[16] = {0}, sumg[16] = {0}, sumc[16] = {0}; int okc[16] = {0};
    for (int c = 0; c < NC; c++) for (int s = 0; s < runs; s++) {
        qp_t Q = P; Q.entangle = cf[c].ent; Q.local_search = cf[c].loc; Q.sa_ref = cf[c].sa; Q.sa_ls = cf[c].sb;
        Q.seed = base + s; Q.include_pg = 0;
        qr_t r; if (qtswe_optimize(&Q, &r)) return 1;
        const qattempt_t *la = &r.att[r.natt - 1];
        printf("%s,%llu,%d,%d,%g,%d,%d,%.3f,%llu,%d,%d,%d\n", cf[c].name, (unsigned long long)Q.seed, r.rows, r.cols, Q.Tg, r.k, r.feasible,
               r.sec, (unsigned long long)r.evals, r.natt, la->gens_used, la->f1);
        fflush(stdout);
        sumk[c] += r.k; sumt[c] += r.sec; okc[c] += r.feasible; sumg[c] += la->gens_used; sumc[c] += la->f1; qtswe_free(&r);
    }
    if (P.fixed_k > 0) {
        fprintf(stderr, "\nfixed k=%d  (success = all cells met within %d generations)\n%-12s %10s %10s %12s\n", P.fixed_k, P.gens, "config", "success", "mean_gens", "mean_covered");
        for (int c = 0; c < NC; c++) fprintf(stderr, "%-12s %7d/%-2d %10.0f %12.1f\n", cf[c].name, okc[c], runs, sumg[c] / runs, sumc[c] / runs);
    } else {
        fprintf(stderr, "\n%-12s %8s %8s %6s\n", "config", "mean_k", "mean_s", "feas");
        for (int c = 0; c < NC; c++) fprintf(stderr, "%-12s %8.2f %8.2f %3d/%d\n", cf[c].name, sumk[c] / runs, sumt[c] / runs, okc[c], runs);
    }
    free(a2); return 0;
}

/* ----------------------------------------------------------------- HTTP */
static void url_decode(char *s) {
    char *o = s;
    for (; *s; s++) {
        if (*s == '+') *o++ = ' ';
        else if (*s == '%' && isxdigit((unsigned char)s[1]) && isxdigit((unsigned char)s[2])) {
            char h[3] = { s[1], s[2], 0 }; *o++ = (char)strtol(h, NULL, 16); s += 2;
        } else *o++ = *s;
    }
    *o = 0;
}
static void send_all(sock_t fd, const char *p, size_t n) {
    while (n) { long w = sock_write(fd, p, n); if (w <= 0) return; p += w; n -= (size_t)w; }
}
static void respond(sock_t fd, int code, const char *status, const char *body) {
    char h[512];
    int n = snprintf(h, sizeof h,
        "HTTP/1.1 %d %s\r\nContent-Type: application/json; charset=utf-8\r\nContent-Length: %zu\r\n"
        "Access-Control-Allow-Origin: *\r\nAccess-Control-Allow-Methods: GET, POST, OPTIONS\r\n"
        "Access-Control-Allow-Headers: Content-Type\r\nAccess-Control-Allow-Private-Network: true\r\n"
        "Cache-Control: no-store\r\nConnection: close\r\n\r\n", code, status, strlen(body));
    send_all(fd, h, n); send_all(fd, body, strlen(body));
}
static void json_err(sock_t fd, int code, const char *status, const char *msg) {
    char b[300]; snprintf(b, sizeof b, "{\"ok\":false,\"error\":\"%s\"}", msg); respond(fd, code, status, b);
}


/* ---------------------------------------------------------- result cache
 * Same parameters (incl. seed) always give the same layout, so a finished run is
 * stored as cache/<hash>.json and replayed instantly.  threads/time are not part
 * of the key (results are thread-count independent); runs that hit the time limit
 * and runs with seed=0 (clock-derived) are never cached. */
static char g_cache_dir[512] = "cache";
static int  g_cache_on = 1;

static void cache_key(char *k, size_t n, const qp_t *P) {
    snprintf(k, n, "v1|%d|%d|%d|%d|%.17g|%.17g|%.17g|%.17g|%.17g|%.17g|%d|%d|%d|%d|%d|%d|%d|%.17g|%.17g|%.17g|"
                   "%d|%d|%d|%.17g|%.17g|%d|%d|%d|%llu|%d",
        P->rows, P->cols, P->sink_r, P->sink_c, P->Tg, P->Rcon, P->Rmax, P->alpha, P->beta, P->Rt,
        P->use_conn, P->pop, P->gens, P->detect, P->local_radius, P->local_rounds, 0, P->delta, P->w1, P->w2,
        P->entangle, P->sa_ref, P->sa_ls, P->sa_T0, P->sa_Tend, P->local_search, P->retries, P->fixed_k,
        (unsigned long long)P->seed, P->include_pg);
}
static void cache_path(char *p, size_t n, const char *key) {
    uint64_t h = 1469598103934665603ULL;
    for (const char *c = key; *c; c++) { h ^= (unsigned char)*c; h *= 1099511628211ULL; }
    snprintf(p, n, "%s/%016llx.json", g_cache_dir, (unsigned long long)h);
}
/* file layout: first line = key, rest = JSON body.  returns malloc'd body or NULL */
static char *cache_load(const char *key) {
    char path[640]; cache_path(path, sizeof path, key);
    FILE *f = fopen(path, "rb"); if (!f) return NULL;
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    if (sz <= 0) { fclose(f); return NULL; }
    char *buf = malloc((size_t)sz + 1);
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) { fclose(f); free(buf); return NULL; }
    fclose(f); buf[sz] = 0;
    char *nl = strchr(buf, '\n');
    if (!nl) { free(buf); return NULL; }
    *nl = 0;
    if (strcmp(buf, key)) { free(buf); return NULL; }   /* hash collision */
    char *body = strdup(nl + 1); free(buf); return body;
}
static void cache_store(const char *key, const char *body) {
    MKDIR(g_cache_dir);
    char path[640], tmp[660]; cache_path(path, sizeof path, key);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "wb"); if (!f) return;
    fprintf(f, "%s\n%s", key, body); fclose(f);
    remove(path); rename(tmp, path);
}

static void handle(sock_t fd) {
    char buf[16384]; size_t n = 0;
    while (n < sizeof buf - 1) {
        long r = sock_read(fd, buf + n, sizeof buf - 1 - n);
        if (r <= 0) break;
        n += r; buf[n] = 0;
        if (strstr(buf, "\r\n\r\n")) break;
    }
    buf[n] = 0;
    char method[16] = {0}, path[8192] = {0};
    if (sscanf(buf, "%15s %8191s", method, path) != 2) { json_err(fd, 400, "Bad Request", "bad request"); return; }
    if (!strcmp(method, "OPTIONS")) { respond(fd, 200, "OK", "{}"); return; }
    char *q = strchr(path, '?'); if (q) *q++ = 0;
    if (!strcmp(path, "/health")) {
        char b[160]; snprintf(b, sizeof b, "{\"ok\":true,\"name\":\"qtswe\",\"threads\":%d,\"max_cells\":4000000}", qtswe_threads());
        respond(fd, 200, "OK", b); return;
    }
    if (strcmp(path, "/optimize")) { json_err(fd, 404, "Not Found", "unknown route"); return; }
    qp_t P; qtswe_defaults(&P); P.time_limit = 120;
    int use_cache = g_cache_on;
    if (q) {
        char *save = NULL;
        for (char *kv = strtok_r(q, "&", &save); kv; kv = strtok_r(NULL, "&", &save)) {
            char *eq = strchr(kv, '='); if (!eq) continue; *eq = 0;
            url_decode(kv); url_decode(eq + 1);
            if (!strcmp(kv, "nocache")) { if (atoi(eq + 1)) use_cache = 0; continue; }
            if (qtswe_set_param(&P, kv, eq + 1)) { char m[120]; snprintf(m, sizeof m, "unknown parameter %.40s", kv); json_err(fd, 400, "Bad Request", m); return; }
        }
    }
    if (P.time_limit <= 0 || P.time_limit > 3600) P.time_limit = 3600;
    const char *e = qtswe_validate(&P);
    if (e) { json_err(fd, 400, "Bad Request", e); return; }
    fprintf(stderr, "[optimize] %dx%d Tg=%.2f conn=%d entangle=%d local=%d seed=%llu\n", P.rows, P.cols, P.Tg, P.use_conn,
            P.entangle, P.local_search, (unsigned long long)P.seed);
    char key[1024]; cache_key(key, sizeof key, &P);
    if (P.seed == 0) use_cache = 0;
    if (use_cache) {
        char *hit = cache_load(key);
        if (hit) {
            sb_t o = {0}; sb_put(&o, "{\"cached\":true,%s", hit + 1);
            fprintf(stderr, "[optimize] cache hit\n");
            respond(fd, 200, "OK", o.s); free(o.s); free(hit); return;
        }
    }
    qr_t r; sb_t b = {0};
    if (qtswe_optimize(&P, &r)) { json_err(fd, 500, "Internal Server Error", "optimize failed"); return; }
    result_json(&b, &P, &r);
    fprintf(stderr, "[optimize] -> k=%d feasible=%d %.2fs\n", r.k, r.feasible, r.sec);
    if (use_cache && !r.timed_out) cache_store(key, b.s);
    { sb_t o = {0}; sb_put(&o, "{\"cached\":false,%s", b.s + 1); respond(fd, 200, "OK", o.s); free(o.s); }
    free(b.s); qtswe_free(&r);
}

static int cmd_serve(int argc, char **argv) {
    int port = 8080; const char *host = "127.0.0.1";
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "--port") && i + 1 < argc) port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--host") && i + 1 < argc) host = argv[++i];
        else if (!strcmp(argv[i], "--cache") && i + 1 < argc) {
            ++i; if (!strcmp(argv[i], "off")) g_cache_on = 0; else snprintf(g_cache_dir, sizeof g_cache_dir, "%s", argv[i]);
        }
    }
#ifdef _WIN32
    WSADATA wsa; if (WSAStartup(MAKEWORD(2, 2), &wsa)) { fprintf(stderr, "WSAStartup failed\n"); return 1; }
#else
    signal(SIGPIPE, SIG_IGN);
#endif
    sock_t s = socket(AF_INET, SOCK_STREAM, 0); int on = 1;
    setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char *)&on, sizeof on);
    struct sockaddr_in a = {0}; a.sin_family = AF_INET; a.sin_port = htons(port);
    if (inet_pton(AF_INET, host, &a.sin_addr) != 1) { fprintf(stderr, "bad host\n"); return 2; }
    if (bind(s, (struct sockaddr *)&a, sizeof a) || listen(s, 16)) { perror("bind/listen"); return 1; }
    fprintf(stderr, "QTSwE backend listening on http://%s:%d  (threads=%d)\n  GET /health\n  GET /optimize?rows=25&cols=25&tg=0.3\n  result cache: %s\n", host, port, qtswe_threads(), g_cache_on ? g_cache_dir : "off");
    for (;;) {
        sock_t c = accept(s, NULL, NULL);
#ifdef _WIN32
        if (c == INVALID_SOCKET) continue;
#else
        if (c < 0) continue;
#endif
        handle(c); sock_close(c);
    }
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "serve")) return cmd_serve(argc - 2, argv + 2);
    if (argc >= 2 && !strcmp(argv[1], "run")) return cmd_run(argc - 2, argv + 2);
    if (argc >= 2 && !strcmp(argv[1], "bench")) return cmd_bench(argc - 2, argv + 2);
    fprintf(stderr,
        "usage: %s serve [--port 8080] [--host 127.0.0.1] [--cache DIR|off]\n"
        "       %s run   --rows M --cols N --tg 0.3 [--json] [--seed S] [--conn 0|1] [--entangle 0|1] [--local 0|1] [--sa_ref 0|1] [--sa_ls 0|1]\n"
        "                [--pop 100 --gens 500 --detect 10 --delta 0.1 --rcon 1.414 --rmax 5 --rt 5 --k FIXED_K --time SEC]\n"
        "       %s bench --rows M --cols N --tg 0.3 --runs 10 [--k FIXED_K → success-rate mode]   (ablation: full / no_entangle / no_local / plain_QTS / +SA_A / +SA_B / +SA_AB)\n",
        argv[0], argv[0], argv[0]);
    return 2;
}
