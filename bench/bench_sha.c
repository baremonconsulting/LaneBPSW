// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jan Suchánek, Baremon s.r.o.
// bench_sha.c: OpenSSL one-shot SHA256() against SHA256_Init/Update/Final on 32-byte messages.
#define OPENSSL_SUPPRESS_DEPRECATED
#include <openssl/sha.h>
#include <stdio.h>
#include <time.h>
static double cpu(void){struct timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec;}
int main(void){ unsigned char in[33]={0}, h[32]; int N=1000000; volatile unsigned x=0;
  for (int rep=0; rep<3; rep++){
  double t0=cpu(); for(int i=0;i<N;i++){ in[31]=i; SHA256(in,32,h); x+=h[0]; }
  double t1=cpu(); for(int i=0;i<N;i++){ in[31]=i; SHA256_CTX c; SHA256_Init(&c); SHA256_Update(&c,in,32); SHA256_Final(h,&c); x+=h[0]; }
  double t2=cpu(); printf("one-shot SHA256(): %.1f ns   SHA256_Init/Update/Final: %.1f ns   ratio %.2f\n", 1e9*(t1-t0)/N, 1e9*(t2-t1)/N, (t1-t0)/(t2-t1)); }
  return 0; }
