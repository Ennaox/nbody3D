all: basenbody.g v1enbody.g

bench:
	taskset -c 3 ./basenbody.g
	taskset -c 3 ./v1enbody.g

basenbody.g: basenbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

v1enbody.g: v1nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

basenbody.i: basenbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

v1nbody.i: v1nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	rm -Rf *~ nbody.g nbody.i *.optrpt *.g

