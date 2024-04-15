CC = gcc
CFLAGS = -Wall -O3 -std=c99 -g

omp: main.c ga.c io.c
	$(CC) $(CFLAGS) main.c ga.c io.c -lm -fopenmp -o omp -DTIMEGenPobIni

run_omp:
	./omp 20000 0 6 1 1 >> ../output/out
	./omp 20000 0 6 2 1 >> ../output/out
	./omp 20000 0 6 4 2 >> ../output/out	

clean:
	rm -f sec
