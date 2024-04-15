CC = gcc
CFLAGS = -fopenmp -O3 -Wall -std=c99 -g

N = 
N_GEN = 
T_POB = 
M_RATE = 

EXEC = sec
OUTFILE = ../output/out.txt
C_FILES = main.c io.c ga-7.c

sec: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG

test_sec:
	./$(EXEC) $(N) $(N_GEN) $(T_POB) $(M_RATE) >> $(OUTFILE)

clean:
	rm -f $(EXEC)
