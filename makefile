all: basenbody.g v1nbody.g v2nbody.g v3nbody256.g v3nbody512.g

bench: base v1 v2 v3_256 v3_512

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

