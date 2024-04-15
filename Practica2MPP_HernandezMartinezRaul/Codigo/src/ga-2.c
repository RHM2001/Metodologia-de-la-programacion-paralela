#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "../include/ga.h"

#define PRINT 1

#define NGM 10
#define NEM 5

int aleatorio(int n) {
	return (rand() % n);  // genera un numero aleatorio entre 0 y n-1
}

int search_element(int *array, int end, int element)
{
	int i=0;
	int found=0;
	
	// comprueba que un elemento no está incluido en el individuo (el cual no admite enteros repetidos)
	while((i < end) && ! found) {
		if(array[i] == element) {
			found = 1;
		}
		i++;
	}
        return found;
}

int find_element(int *array, int end, int element)
{
        int pos = 0;
	for(int i = 0; i < end; i++) {
             if(array[i] == element) {
                 pos = i;
                 break;
             }
        }
        return pos; // Posición del elemento encontrado
}

int *crear_individuo(int n)
{
        // El primer elemento del individuo siempre será el 0, por ejemplo.
	int i=1, value;
	int *individuo = (int *) malloc(n * sizeof(int));
	
	// inicializa array de elementos
	memset(individuo, 0, n * sizeof(int));
	
	while(i < n) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!search_element(individuo, i, value)) {
			individuo[i] = value;  // lo incluimos
			i++;
		}
	}
	return individuo;
}

int comp_fitness(const void *a, const void *b) {
	/* qsort pasa un puntero al elemento que está ordenando */
	return (*(Individuo *)a).fitness - (*(Individuo *)b).fitness;
}

/**
 *  Crea un nuevo tipo de datos derivado en MPI
 *  Necesario para el envio/recepcion de mensajes con datos de tipo Individuo
 **/
void crear_tipo_datos(int n, MPI_Datatype *individuo_type)
{
	int blocklen[2] = {n, 1};
	MPI_Datatype dtype[2] = {MPI_INT, MPI_DOUBLE};

	MPI_Aint disp[2];
	disp[0] = offsetof(Individuo, array_int);
	disp[1] = offsetof(Individuo, fitness);

	MPI_Type_create_struct(2, blocklen, disp, dtype, individuo_type);
	MPI_Type_commit(individuo_type);
}

double aplicar_ga(const double *d, int n, int n_gen, int tam_pob, int *sol, double m_rate, int rank, int size)
{
	int i, g, mutation_start;

	// Crear tipo de datos derivado MPI para Individuo
	MPI_Datatype individuo_type;
	crear_tipo_datos(n, &individuo_type);

	Individuo *isla_poblacion;

	// Dividimos la población global entre los procesos
	int isla_poblacion_size = tam_pob / size;
	// Nos guardamos el resto
	int resto = tam_pob % size;
	
	if(rank == 0){

		// crea poblacion inicial (array de individuos)
		isla_poblacion = (Individuo *)malloc(tam_pob * sizeof(Individuo));
		assert(isla_poblacion);

		// crea cada individuo (array de enteros aleatorios)
		for (i = 0; i < tam_pob; i++)
		{
			//QUITAR MALLOC
			//sub_poblacion[i] = (Individuo *)malloc(sizeof(Individuo));
			int *array = crear_individuo(n);
			for (int j = 0; j < n; j++) {
				isla_poblacion[i].array_int[j] = array[j];
			}
			free(array);

			// calcula el fitness del individuo
			fitness(d, &isla_poblacion[i], n);
		}

		// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
		qsort(isla_poblacion, tam_pob, sizeof(Individuo), comp_fitness);

		// Le enviamos a las islas la parte del array correspondiente a partir de la isla 0
		for(i=1; i < size; i++){
			MPI_Send(isla_poblacion+(isla_poblacion_size*i)+resto, isla_poblacion_size, individuo_type, i, 0, MPI_COMM_WORLD);
		}
		isla_poblacion_size = isla_poblacion_size + resto;

	}
	else{
		isla_poblacion = (Individuo *)malloc(isla_poblacion_size * sizeof(Individuo));
		assert(isla_poblacion);
		// Recibimos la sub-poblacion
		MPI_Recv(isla_poblacion, isla_poblacion_size, individuo_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	//double fitness_anterior = sub_poblacion[0]->fitness;

	int contador_NGM = 0;

	// evoluciona la poblacion durante un numero de generaciones
	for (g = 0; g < n_gen; g++)
	{

		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		for (i = 0; i < (isla_poblacion_size / 2) - 1; i += 2)
		{
			cruzar(&isla_poblacion[i], &isla_poblacion[i + 1], &isla_poblacion[isla_poblacion_size / 2 + i], &isla_poblacion[isla_poblacion_size / 2 + i + 1], n);
		}

		// por ejemplo, inicia la mutacion a partir de 1/4 de la poblacion.
		// puede haber otras opciones pero dejar el primer individuo sin modificar siempre
		mutation_start = isla_poblacion_size / 4;

		// muta 3/4 partes de la poblacion
		for (i = mutation_start; i < isla_poblacion_size; i++)
		{
			mutar(&isla_poblacion[i], n, m_rate);
		}

		// recalcula el fitness del individuo
		for (i = 0; i < isla_poblacion_size; i++)
		{
			fitness(d, &isla_poblacion[i], n);
		}

		// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
		qsort(isla_poblacion, isla_poblacion_size, sizeof(Individuo), comp_fitness);

		if (PRINT && rank == 0)
		{
			printf("Generacion %d - ", g);
			printf("Fitness = %.0lf\n", (isla_poblacion[0].fitness));
		}

		if (NGM == contador_NGM){
			if (rank == 0){
				isla_poblacion_size = isla_poblacion_size - resto;

				for(i=1 ; i < size ; i++){
					// Obtenemos los mejores individuos de cada poblacion de cada isla
					MPI_Recv(isla_poblacion + (isla_poblacion_size*i) + resto, NEM, individuo_type, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}

				qsort(isla_poblacion, tam_pob, sizeof(Individuo), comp_fitness);

				for(i=1 ; i < size ; i++){
					// Se envian a las islas los mejores individuos de todos los que hay (de las demas islas)
					MPI_Send(isla_poblacion, NEM, individuo_type, i, 0, MPI_COMM_WORLD);
				}
				isla_poblacion_size = isla_poblacion_size + resto;
			}
			else{
				// Se envian a la isla 0 los mejores individuos
				MPI_Send(isla_poblacion, NEM, individuo_type, 0, 0, MPI_COMM_WORLD);
				// Se reciben los mejores individuos que hay
				MPI_Recv(isla_poblacion, NEM, individuo_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			contador_NGM = 0;
		}
		else{
			contador_NGM = contador_NGM + 1;
		}
	}

	if (rank == 0){
		isla_poblacion_size = isla_poblacion_size - resto;

		for(i=1 ; i < size ; i++){
			// Obtenemos los mejores individuos de cada poblacion de cada isla
			MPI_Recv(isla_poblacion + (isla_poblacion_size*i) + resto, NEM, individuo_type, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		isla_poblacion_size = isla_poblacion_size + resto;

		qsort(isla_poblacion, tam_pob, sizeof(Individuo), comp_fitness);
		
	}
	else{
		// Se envian a la isla 0 los mejores individuos
		MPI_Send(isla_poblacion, NEM, individuo_type, 0, 0, MPI_COMM_WORLD);
	}

	memmove(sol, isla_poblacion[0].array_int, n*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (isla_poblacion[0].fitness);
	
	free(isla_poblacion);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n)
{
	// Elegir un punto (o puntos) de corte aleatorio a partir del que se realiza el intercambio de los genes. 
	int puntoCorteAleatorio = aleatorio(n-1);
	
	// Entonces, por ejemplo, los primeros genes del padre1 van al hijo1, y los primeros del padre2 al hijo2.
        // Se debe evitar en cada paso la introduccion de duplicados en los hijos
	// Y los restantes genes de cada hijo son del otro padre, respectivamente.

	for(int i = 0 ; i<puntoCorteAleatorio ; i++){
		hijo1->array_int[i] = padre1->array_int[i];
		hijo2->array_int[i] = padre2->array_int[i];
	}

	for(int i = puntoCorteAleatorio ; i<n ; i++){
		hijo1->array_int[i] = padre2->array_int[i];
		hijo2->array_int[i] = padre1->array_int[i];
	}

        // Otra opción podría ser factibilizar a posteriori, despues de generar los descendientes: eliminar posibles 
        // repetidos de ambos hijos. Si encuentro algún elemento repetido en el hijo, lo cambio por otro que no este el array

}

void invertir(int *a, int k)
{
        //int t;
	// Uno por uno invierte los elementos de a[0..k-1]
}

void mutar(Individuo *actual, int n, double m_rate) {
    int i, j;

    srand(time(NULL));

    // Seleccionar dos posiciones aleatorias (i, j) tal que i < j
    i = rand() % (n - 1) + 1; // i debe ser > 0 y < n
    j = rand() % (n - i) + i + 1;

    // Realizar la mutación: invertir el orden de los elementos entre i y j
    while (i < j) {

	// Determinar si mutar o no basado en la tasa de mutación
    double random_value = ((double)rand() / RAND_MAX);
    if (random_value > m_rate) {
        // Intercambiar los elementos en las posiciones i y j
        int temp = actual->array_int[i];
        actual->array_int[i] = actual->array_int[j];
        actual->array_int[j] = temp;
    }
        // Mover los índices hacia el centro
        i++;
        j--;
    }
}

double distancia_ij(const double *d, int i, int j, int n)
{
	double dist = 0.0;

	if (i >= 0 && i < n && j >= 0 && j < n && i != j)  
    {
        dist = d[i*n + j];
    }
    return dist; 
	
	// Devuelve la distancia entre dos elementos de la matriz 'd'
	return dist;
}

void fitness(const double *d, Individuo *individuo, int n)
{
	double temp_fitness = 0.0;

	for(int i=0 ; i < n ; i++){
		double distancia = distancia_ij(d, individuo->array_int[i], individuo->array_int[i + 1], n);
		temp_fitness += distancia;
	}

	// Finalmente, se cierra el circulo sumando la distancia de la última ciudad con la ciudad inicial
    double  distancia_ult = distancia_ij(d, individuo->array_int[n-1], individuo->array_int[0], n);
    temp_fitness += distancia_ult;

    // Se establece el valor fitness calculado al individuo
    individuo->fitness = temp_fitness;
}
