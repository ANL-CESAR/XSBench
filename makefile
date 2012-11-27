all: do_flops.o
	gcc -fopenmp -Wall -std=c99 -O3 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c do_flops.o -o XSBench -lm
papi: do_flops.o
	gcc -fopenmp -Wall -std=c99 -O3 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c do_flops.o /usr/local/lib/libpapi.a -o XSBench -lm
clean:
	rm -rf XSBench XSBench.dSYM do_flops.o
run:
	./XSBench
edit: 
	vim -p Main.c GridInit.c XSutils.c CalculateXS.c Materials.c papi.c XSbench_header.h do_flops.c
profile: do_flops.o
	gcc -fopenmp -Wall -pg -std=c99 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c do_flops.o -o XSBench -lm
debug: do_flops.o
	gcc -fopenmp -Wall -std=c99 -g Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c do_flops.o -o XSBench -lm -g
do_flops.o: XSbench_header.h
	gcc -c -Wall -std=c99 -O0 do_flops.c
