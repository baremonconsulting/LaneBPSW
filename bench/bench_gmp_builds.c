// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// bench_gmp_builds.c: the same GMP functions from two shared libraries (dlopen), rotated order, thread CPU time.
// usage: bench_gmp_builds libA libB  (prints time ratios A/B)
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
typedef void (*powm_t)(mpz_ptr, mpz_srcptr, mpz_srcptr, mpz_srcptr);
typedef int (*mr_t)(mpz_srcptr, int);
typedef void (*np_t)(mpz_ptr, mpz_srcptr);
static double cpu(void){struct timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec;}
int main(int argc, char **argv){
  if (argc < 3) { fprintf(stderr, "usage: bench_gmp_builds libA libB\n"); return 2; }
  void *h[2]; const char *nm[2] = {argv[1], argv[2]};
  powm_t P[2]; mr_t M[2]; np_t N[2];
  for (int i=0;i<2;i++){ h[i]=dlopen(nm[i], RTLD_NOW|RTLD_LOCAL); if(!h[i]){printf("dlopen %s\n", dlerror()); return 1;}
    P[i]=(powm_t)dlsym(h[i],"__gmpz_powm"); M[i]=(mr_t)dlsym(h[i],"__gmpz_millerrabin"); N[i]=(np_t)dlsym(h[i],"__gmpz_nextprime"); }
  /* use library 0's mpz functions for setup: each library has its own allocator, so build numbers per library */
  typedef void (*init_t)(mpz_ptr); typedef void (*setstr_t)(mpz_ptr,const char*,int); typedef void (*setui_t)(mpz_ptr, unsigned long);
  typedef void (*subui_t)(mpz_ptr, mpz_srcptr, unsigned long);
  int sizes[] = {256, 1024, 2048, 4096}; uint64_t s = 12345;
  for (int si=0; si<4; si++){ int k = sizes[si]; int reps = k<=1024? 60 : (k<=2048? 20 : 6);
    double tp[2]={0,0}, tm[2]={0,0}, tn[2]={0,0};
    for (int r=0;r<reps;r++){
      char hex[1100]; int nd = k/4; for (int j=0;j<nd;j++){ s = s*6364136223846793005ULL+1442695040888963407ULL; hex[j]="0123456789abcdef"[(s>>59)&15]; } hex[0]="89abcdef"[(s>>40)&7]; hex[nd-1]='1'+2*((s>>50)&3); hex[nd]=0;
      for (int i=0;i<2;i++){ int o = (r+i)&1;  /* rotate order */
        init_t in=(init_t)dlsym(h[o],"__gmpz_init"); setstr_t ss=(setstr_t)dlsym(h[o],"__gmpz_set_str"); setui_t su=(setui_t)dlsym(h[o],"__gmpz_set_ui"); subui_t sb=(subui_t)dlsym(h[o],"__gmpz_sub_ui");
        mpz_t n, p, two, e, y; in(n); in(p); in(two); in(e); in(y); ss(n, hex, 16); su(two, 2);
        double t0=cpu(); N[o](p, n); double t1=cpu(); M[o](p, 25); double t2=cpu(); sb(e, p, 1); P[o](y, two, e, p); double t3=cpu();
        tn[o]+=t1-t0; tm[o]+=t2-t1; tp[o]+=t3-t2; }
    }
    printf("k=%4d  nextprime ms: sys %.3f fat %.3f (sys/fat %.2f) | millerrabin(p,25) us: sys %.1f fat %.1f (%.2f) | powm us: sys %.1f fat %.1f (%.2f)\n",
      k, 1e3*tn[0]/reps, 1e3*tn[1]/reps, tn[0]/tn[1], 1e6*tm[0]/reps, 1e6*tm[1]/reps, tm[0]/tm[1], 1e6*tp[0]/reps, 1e6*tp[1]/reps, tp[0]/tp[1]);
  }
  return 0;
}
