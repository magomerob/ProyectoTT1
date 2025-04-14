// $Source$
/**
 * @file eccAnom.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación del cálculo de la anomalía excéntrica.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/eccAnom.h"

double eccAnom(double M, double e)
{
    int maxit = 15;
    int i = 1;

    M = fmod(M, 2.0*M_PI);

    double E = 0;

	if (e<0.8){
		E = M;
	}else{
		E = M_PI;
	}

	double f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );
	
	double eps=1e-10;
	
	while (abs(f) > 1e2*eps){    
		f = E - e*sin(E) - M;
		E = E - f / ( 1.0 - e*cos(E) );
		i = i+1;
		if (i==maxit){
			cout<<" convergence problems in EccAnom";
			exit(EXIT_FAILURE);
		}  
	}
	
	return E;
}