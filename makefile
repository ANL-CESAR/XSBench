all:
	gcc -fopenmp -Wall -std=c99 -O3 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c /usr/local/lib/libpapi.a -o XSBench -lm
clean:
	rm -rf XSBench XSBench.dSYM
run:
	./XSBench
edit:
	vim -p Main.c GridInit.c XSutils.c CalculateXS.c Materials.c papi.c XSbench_header.h
profile:
	gcc -fopenmp -Wall -pg -std=c99 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c -o XSBench -lm
debug:
	gcc -fopenmp -Wall -std=c99 -g Materials.c Main.c GridInit.c XSutils.c CalculateXS.c papi.c -o XSBench -lm -g
