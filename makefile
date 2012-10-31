all:
	gcc -fopenmp -Wall -std=c99 -O3 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c -o XSBench -lm
clean:
	rm -f XSBench
run:
	./XSBench
edit:
	vim -p Main.c GridInit.c XSutils.c CalculateXS.c Materials.c XSBench_header.h
