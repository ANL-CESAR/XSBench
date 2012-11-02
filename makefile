all:
	gcc -fopenmp -Wall -std=c99 -O3 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c -o XSBench -lm
clean:
	rm -f XSBench
run:
	./XSBench
edit:
	vim -p Main.c GridInit.c XSutils.c CalculateXS.c Materials.c XSbench_header.h
profile:
	gcc -fopenmp -Wall -pg -std=c99 Materials.c Main.c GridInit.c XSutils.c CalculateXS.c -o XSBench -lm
	
