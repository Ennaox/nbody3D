//PAS FINI

//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

  f32* restrict fx = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict fy = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict fz = aligned_alloc(ALIGN,sizeof(f32) * 8);

  f32* restrict pxi = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict pyi = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict pzi = aligned_alloc(ALIGN,sizeof(f32) * 8);  

  f32* restrict dx = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict dy = aligned_alloc(ALIGN,sizeof(f32) * 8);
  f32* restrict dz = aligned_alloc(ALIGN,sizeof(f32) * 8);
  //
  for (u64 i = 0; i < n; i+=8)
    {
      //

      fx[0] = 0.0; fy[0] = 0.0; fz[0] = 0.0;
      fx[1] = 0.0; fy[1] = 0.0; fz[1] = 0.0;
      fx[2] = 0.0; fy[2] = 0.0; fz[2] = 0.0;
      fx[3] = 0.0; fy[3] = 0.0; fz[3] = 0.0;
      fx[4] = 0.0; fy[4] = 0.0; fz[4] = 0.0;
      fx[5] = 0.0; fy[5] = 0.0; fz[5] = 0.0;
      fx[6] = 0.0; fy[6] = 0.0; fz[6] = 0.0;
      fx[7] = 0.0; fy[7] = 0.0; fz[7] = 0.0;

      //23 floating-point operations
      pxi[0] = p->x[i]; pyi[0] = p->y[i]; pzi[0] = p->z[i];
      pxi[1] = p->x[i+1]; pyi[1] = p->y[i+1]; pzi[1] = p->z[i+1];
      pxi[2] = p->x[i+2]; pyi[2] = p->y[i+2]; pzi[2] = p->z[i+2];
      pxi[3] = p->x[i+3]; pyi[3] = p->y[i+3]; pzi[3] = p->z[i+3];
      pxi[4] = p->x[i+4]; pyi[4] = p->y[i+4]; pzi[4] = p->z[i+4];
      pxi[5] = p->x[i+5]; pyi[5] = p->y[i+5]; pzi[5] = p->z[i+5];
      pxi[6] = p->x[i+6]; pyi[6] = p->y[i+6]; pzi[6] = p->z[i+6];
      pxi[7] = p->x[i+7]; pyi[7] = p->y[i+7]; pzi[7] = p->z[i+7];

      for (u64 j = 0; j < n; j++)
	{
	  //Newton's law
	  dx = p->x[j] - pxi; dy = p->y[j] - pyi; dz = p->z[j] - pzi;
	  
    const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9
    f32 tmp = sqrt(d_2);  //10
    const f32 d_3_over_2 = tmp * tmp * tmp; //12

	  //Net force
	  fx += dx / d_3_over_2; //14
	  fy += dy / d_3_over_2; //16
	  fz += dz / d_3_over_2; //18
	}

      //
      p->vx[i] += dt * fx; //20
      p->vy[i] += dt * fy; //22
      p->vz[i] += dt * fz; //24
    }

  //3 floating-point operations
  for (u64 i = 0; i < n; i++)
    {
      p->x[i] += dt * p->vx[i];
      p->y[i] += dt * p->vy[i];
      p->z[i] += dt * p->vz[i];
    }
    printf("%f\n",p->x[2]);
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