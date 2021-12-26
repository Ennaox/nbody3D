//
#include <stdio.h>
#if __AVX512F__

#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

#define ALIGN 64

//
typedef struct particle_s {

  f32 * restrict x,* restrict y,* restrict z;
  f32 * restrict vx,* restrict vy,* restrict vz;
  
} particle_t;

//
void init(particle_t *p, u64 n)
{
  for (u64 i = 0; i < n; i++)
    {
      //
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->y[i] = (f32)rand() / (f32)RAND_MAX;
      p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

      //
      p->vx[i] = (f32)rand() / (f32)RAND_MAX;
      p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t *p, const f32 dt, u64 n)
{
  //
  const f32 softening = 1e-20;

  __m512 vsoftening = _mm512_load_ps(&softening);

  __m512 pxi, pyi, pzi, pxj, pyj, pzj, vxi, vyi, vzi;

  __m512 vdt = _mm512_load_ps(&dt);

  //
  for (u64 i = 0; i < n; i++)
    {
      pxi = _mm512_loadu_ps(&p->x[i]);
      pyi = _mm512_loadu_ps(&p->y[i]);
      pzi = _mm512_loadu_ps(&p->z[i]);

      vxi = _mm512_loadu_ps(&p->vx[i]);
      vyi = _mm512_loadu_ps(&p->vy[i]);
      vzi = _mm512_loadu_ps(&p->vz[i]);

      //
      __m512 vfx = _mm512_setzero_ps();
      __m512 vfy = _mm512_setzero_ps();
      __m512 vfz = _mm512_setzero_ps();
      //23 floating-point operations
    for (u64 j = 0; j < n; j+= 16)
    {
      pxj = _mm512_loadu_ps(&p->x[j]);
      pyj = _mm512_loadu_ps(&p->y[j]);
      pzj = _mm512_loadu_ps(&p->z[j]);
      //Newton's law
      const __m512 dx = _mm512_sub_ps(pxj,pxi); //1
      const __m512 dy = _mm512_sub_ps(pyj,pyi); //2
      const __m512 dz = _mm512_sub_ps(pzj,pzi); //3
      const __m512 d_2 = _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(dx,dx),_mm512_mul_ps(dy,dy)),_mm512_add_ps(_mm512_mul_ps(dz,dz),vsoftening)); //9
      #if __AVX512ER__
          __m512 tmp = _mm512_rsqrt28_ps(d_2);
        #else
          __m512 tmp = _mm512_rsqrt14_ps(d_2);
        #endif  //10
      const __m512 d_3_over_2 = _mm512_mul_ps(_mm512_mul_ps(tmp,tmp),tmp); //12

      //Net force
      vfx = _mm512_fmadd_ps(dx,d_3_over_2,vfx); //14
      vfy = _mm512_fmadd_ps(dy,d_3_over_2,vfy); //16
      vfz = _mm512_fmadd_ps(dz,d_3_over_2,vfz); //18
    }
      //
      vxi = _mm512_fmadd_ps(vdt,vfx,vxi); //20
      vyi = _mm512_fmadd_ps(vdt,vfy,vyi); //22
      vzi = _mm512_fmadd_ps(vdt,vfz,vzi); //24
       
      _mm512_storeu_ps(&p->vx[i],vxi);
      _mm512_storeu_ps(&p->vy[i],vyi);
      _mm512_storeu_ps(&p->vz[i],vzi);
  }
  //3 floating-point operations
  for (u64 i = 0; i < n; i+=8)
    {
      pxi = _mm512_loadu_ps(&p->x[i]);
      pyi = _mm512_loadu_ps(&p->y[i]);
      pzi = _mm512_loadu_ps(&p->z[i]);

      vxi = _mm512_loadu_ps(&p->vx[i]);
      vyi = _mm512_loadu_ps(&p->vy[i]);
      vzi = _mm512_loadu_ps(&p->vz[i]);


      vxi = _mm512_fmadd_ps(vdt,vxi,pxi);
      vyi = _mm512_fmadd_ps(vdt,vyi,pyi);
      vzi = _mm512_fmadd_ps(vdt,vzi,pzi);
      _mm512_storeu_ps(&p->x[i],pxi);
      _mm512_storeu_ps(&p->y[i],pyi);
      _mm512_storeu_ps(&p->z[i],pzi); 
    }
}

//
int main(int argc, char **argv)
{
  //
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 steps= 10;
  const f32 dt = 0.01;

  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t *p = malloc(sizeof(particle_t));
  p->x = aligned_alloc(ALIGN,sizeof(f32) * n);
  p->y = aligned_alloc(ALIGN,sizeof(f32) * n);
  p->z = aligned_alloc(ALIGN,sizeof(f32) * n);
  p->vx = aligned_alloc(ALIGN,sizeof(f32) * n);
  p->vy = aligned_alloc(ALIGN,sizeof(f32) * n);
  p->vz = aligned_alloc(ALIGN,sizeof(f32) * n);

  //
  init(p, n);

  const u64 s = sizeof(particle_t) * n;
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  //
  for (u64 i = 0; i < steps; i++)
    {
      //Measure
      const f64 start = omp_get_wtime();

      move_particles(p, dt, n);

      const f64 end = omp_get_wtime();

      //Number of interactions/iterations
      const f32 h1 = (f32)(n) * (f32)(n - 1);

      //GFLOPS
      const f32 h2 = (24.0 * h1 + 3.0 * (f32)n) * 1e-9;
      
      if (i >= warmup)
  {
    rate += h2 / (end - start);
    drate += (h2 * h2) / ((end - start) * (end - start));
  }

      //
      printf("%5llu %10.3e %10.3e %8.1f %s\n",
       i,
       (end - start),
       h1 / (end - start),
       h2 / (end - start),
       (i < warmup) ? "*" : "");
      
      fflush(stdout);
    }

  //
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
   "Average performance:", "", rate, drate);
  printf("-----------------------------------------------------\n");
  
  //
  free(p);

  //
  return 0;
}

#endif

int main()
{
  printf("No AVX512 on this computer\n");
}