#ifndef _GA
#define _GA
	
	typedef struct {
		int *array_int;
		double fitness;
	} Individuo;
	
	void cruzar(Individuo *, Individuo *, Individuo *, Individuo *, int);
	void fitness(const double *, Individuo *, int);
	void mutar(Individuo *, int, double);
#endif
