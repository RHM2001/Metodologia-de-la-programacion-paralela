CC = mpicc
CFLAGS = -O3 -Wall -std=c99 -g

N = 2500
N_GEN = 50
T_POB = 800
M_RATE = 0.15

EXEC = sec
OUTFILE = ../output/out.txt
C_FILES = main.c io.c ga-5.c

all: $(EXEC)

$(EXEC): $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG

run: $(EXEC)
	mpirun -np 6 $(EXEC) $(N) $(N_GEN) $(T_POB) $(M_RATE) >> $(OUTFILE)

clean:
	rm -f $(EXEC)
