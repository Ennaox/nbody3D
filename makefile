all: basenbody.g v1nbody.g v2nbody.g v3nbody256.g v3nbody512.g

compile_all: basenbody.g v1nbody.g v2nbody.g v3nbody256.g v3nbody512.g basenbody.i v1nbody.i v2nbody.i v3nbody256.i v3nbody512.i 

bench: base v1 v2 v3_256 v3_512

bench_all: base_all v1_all v2_all v3_256_all v3_512_all

base: basenbody.g 
	taskset -c 3 ./basenbody.g

v1: v1nbody.g
	taskset -c 3 ./v1nbody.g

v2: v2nbody.g
	taskset -c 3 ./v2nbody.g

v3_256: v3nbody256.g
	taskset -c 3 ./v3nbody256.g

v3_512: v3nbody512.g
	taskset -c 3 ./v3nbody512.g

base_all: basenbody.g basenbody.i 
	taskset -c 3 ./basenbody.g
	taskset -c 3 ./basenbody.i

v1_all: v1nbody.g v1nbody.i
	taskset -c 3 ./v1nbody.g
	taskset -c 3 ./v1nbody.i

v2_all: v2nbody.g v2nbody.i
	taskset -c 3 ./v2nbody.g
	taskset -c 3 ./v2nbody.i

v3_256_all: v3nbody256.g v3nbody256.i
	taskset -c 3 ./v3nbody256.g
	taskset -c 3 ./v3nbody256.i

v3_512_all: v3nbody512.g v3nbody512.i
	taskset -c 3 ./v3nbody512.g
	taskset -c 3 ./v3nbody512.i

basenbody.g: basenbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

v1nbody.g: v1nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

v2nbody.g: v2nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

v3nbody256.g: v3nbody256.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp -g

v3nbody512.g: v3nbody512.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp -g

basenbody.i: basenbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

v1nbody.i: v1nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

v2nbody.i: v2nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

v3nbody256.i: v3nbody256.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

v3nbody512.i: v3nbody512.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	rm -f *~ *.g *.i *.optrpt

