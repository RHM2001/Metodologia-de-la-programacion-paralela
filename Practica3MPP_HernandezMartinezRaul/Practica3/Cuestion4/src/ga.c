#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#include "../include/ga.h"

#define PRINT 1
#define NUM_HILOS 16

unsigned int seed;
#pragma omp threadprivate(seed)
double t_total_time;

static double mseconds() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec*1000 + t.tv_usec/1000;
}


void inicializar_seed() {

	#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        unsigned int thread_seed = (unsigned int)(time(NULL) + thread_id);
        seed = thread_seed;
    }
}

int aleatorio(int n) {
	return (rand_r(&seed) % n);  // genera un numero aleatorio entre 0 y n-1
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
	return (*(Individuo **)a)->fitness - (*(Individuo **)b)->fitness;
}

double aplicar_ga(const double *d, int n, int n_gen, int tam_pob, int *sol, double m_rate, int n_hilos_ini, int n_hilos_fit)
{
	int i, g, mutation_start;
	
	// crea poblacion inicial (array de individuos)
	Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
	assert(poblacion);

	//omp_set_num_threads(NUM_HILOS);

	inicializar_seed();
	
	// crea cada individuo (array de enteros aleatorios)

	//Cuestion 3.1.c)
	double t_ini = mseconds();

	//Cuestion 3.1.a)
	//#pragma omp parallel for private(i)
	//Cuestion 3.1.b
	omp_set_nested(1);
	#pragma omp parallel for private(i) num_threads(n_hilos_ini)
	for(i = 0; i < tam_pob; i++) {
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		poblacion[i]->array_int = crear_individuo(n);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n, n_hilos_fit);
	}

	//Cuestion 3.1.c)
	double t_fin = mseconds();
    t_total_time += (t_fin - t_ini) / 1000;
	printf("Tiempo de FOR_INI = %.4lf\n", t_total_time);
	

	
	// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);

	//double fitness_anterior = poblacion[0]->fitness;

	double temp_total = 0.0;
	
	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen; g++)
	{

		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		double t_ini = mseconds();
		
		
		for(i = 0; i < (tam_pob/2) - 1; i += 2) {
			cruzar(poblacion[i], poblacion[i+1], poblacion[tam_pob/2 + i], poblacion[tam_pob/2 + i + 1], n);
		}
		double t_fin = mseconds();
   		t_total_time += (t_fin - t_ini) / 1000;
		temp_total += t_total_time;
		
		// por ejemplo, inicia la mutacion a partir de 1/4 de la poblacion.
        // puede haber otras opciones pero dejar el primer individuo sin modificar siempre
		mutation_start = tam_pob/4;

		for(i = mutation_start; i < tam_pob; i++) {
			mutar(poblacion[i], n, m_rate);
		}

		for(i = 0; i < tam_pob; i++) {
			fitness(d, poblacion[i], n, n_hilos_fit);
		}
		
		// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
		qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
		
		if (PRINT) {
			printf("Generacion %d - ", g);
			printf("Fitness = %.0lf\n", (poblacion[0]->fitness));

			//printf("Diferencia con Fitness Anterior = %.2e%c\n", mejora_porcentual, 37);
		}

		/*if (mejora_porcentual < m_rate) {
            break;
        }*/

		// se guarda el anterior fitness para la siguiente iteración.

        //fitness_anterior = poblacion[0]->fitness;

	}

	//printf("Tiempo total for cruzar = %.4lf\n", temp_total);
	
	//printf("Tiempo total fitness=%.4lf\n", t_total_time);
	memmove(sol, poblacion[0]->array_int, n*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0]->fitness);
	
	// se libera la memoria reservada
	for(i = 0; i < tam_pob; i++) {
		free(poblacion[i]->array_int);
		free(poblacion[i]);
	}
	
	free(poblacion);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n)
{
	// Elegir un punto (o puntos) de corte aleato   rio a partir del que se realiza el intercambio de los genes. 
	int puntoCorteAleatorio = aleatorio(n-1);
	
	// Entonces, por ejemplo, los primeros genes del padre1 van al hijo1, y los primeros del padre2 al hijo2.
        // Se debe evitar en cada paso la introduccion de duplicados en los hijos
	// Y los restantes genes de cada hijo son del otro padre, respectivamente.
	
	// -------------------------------------------- EJERCICIO 6 ------------------------------------------------
	
	#pragma omp parallel sections shared(padre1, padre2, hijo1, hijo2, puntoCorteAleatorio, n)
	{
		#pragma omp section
		{
			for(int i = 0 ; i<puntoCorteAleatorio ; i++){
				hijo1->array_int[i] = padre1->array_int[i];
				hijo2->array_int[i] = padre2->array_int[i];
			}
		}

		#pragma omp section
		{
			for(int i = puntoCorteAleatorio ; i<n ; i++){
				hijo1->array_int[i] = padre2->array_int[i];
				hijo2->array_int[i] = padre1->array_int[i];
			}
		}
		
	}
        // Otra opción podría ser factibilizar a posteriori, despues de generar los descendientes: eliminar posibles 
        // repetidos de ambos hijos. Si encuentro algún elemento repetido en el hijo, lo cambio por otro que no este el array

}

void invertir(int *a, int k)
{
        int t;
	// Uno por uno invierte los elementos de a[0..k-1]
}

void mutar(Individuo *actual, int n, double m_rate) {
    int i, j;

    //srand(time(NULL));

    // Seleccionar dos posiciones aleatorias (i, j) tal que i < j
    i = rand_r(&seed) % (n - 1) + 1; // i debe ser > 0 y < n
    j = rand_r(&seed) % (n - i) + i + 1;

    // Realizar la mutación: invertir el orden de los elementos entre i y j
    while (i < j) {

	// Determinar si mutar o no basado en la tasa de mutación
    double random_value = ((double)rand_r(&seed) / RAND_MAX);
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

//REDUCTION

void fitness(const double *d, Individuo *individuo, int n, int n_hilos_fit)
{
	double t_ini = mseconds();

    double temp_fitness = 0.0;
	
    
	// Se suma la distancia entre todas las ciudades para calcular su valor fitness

	#pragma omp parallel for reduction(+:temp_fitness) num_threads(n_hilos_fit)
    for(int i = 0; i < n-1; i++) {
		double distancia = distancia_ij(d, individuo->array_int[i], individuo->array_int[i + 1], n);
        temp_fitness += distancia;
    }

	// Se realiza a parte la ultima iteracion para unir la ultima ciudad con la primera.
	double distancia_ult = distancia_ij(d, individuo->array_int[n-1], individuo->array_int[0], n);
	temp_fitness += distancia_ult;


	// Se guarda el fitness en el fitness del individuo
	individuo->fitness = temp_fitness;

	double t_fin = mseconds();
    t_total_time += (t_fin - t_ini) / 1000;
}
