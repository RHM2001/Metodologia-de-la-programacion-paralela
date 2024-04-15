#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

#include "../include/io.h"

extern double aplicar_ga(const double *, int, int, int, int *, double, int, int);

static double mseconds() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec*1000 + t.tv_usec/1000;
}

int main(int argc, char **argv)
{
//	Check Number of Input Args
	if(argc < 4) {
		fprintf(stderr,"Ayuda:\n"); 
		fprintf(stderr,"  ./programa n nGen tamPob mRate\n");
		return(EXIT_FAILURE);
	}
	
	int n;
	int n_gen;
	int tam_pob;
	double m_rate;
	
//	Generate matrix D with distance values among elements
	double *d;

	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// El cero lee los argumentos
	if (rank == 0){
		n = atoi(argv[1]);
		n_gen = atoi(argv[2]);
		tam_pob = atoi(argv[3]);
		m_rate = atof(argv[4]);
	}
	
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tam_pob, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if( rank == 0){
		d = generar_matriz_distancias(n);
	}
	else {
		d = (double *) malloc(n * n *sizeof(double));
	}

	MPI_Bcast(d, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	#ifdef DEBUG
//		print_matrix(d, n);
	#endif
	
//	Allocate memory for output data
	int *sol = (int *) malloc(n * sizeof(int));
	
	#ifdef TIME
		MPI_Barrier(MPI_COMM_WORLD);
		double ti = 0;
		if(rank == 0){
			ti = mseconds();
		}
	#endif
	
//	Call Metaheuristic
	double value = aplicar_ga(d, n, n_gen, tam_pob, sol, m_rate, rank, size);
	
	#ifdef TIME
		MPI_Barrier(MPI_COMM_WORLD);
		double tf = 0;
		if(rank == 0){
			tf = mseconds();
			printf("Execution Time: %.2lf sec\n", (tf - ti)/1000);
		}
	#endif
	
    
	#ifdef DEBUG
		if(rank == 0){
			print_solution(n, sol, value);
		}
	#endif

	    // Free Allocated Memory	
	free(sol);
	free(d);

	MPI_Finalize();
	return(EXIT_SUCCESS);
}
