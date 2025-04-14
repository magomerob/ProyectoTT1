// $Source$
/**
 * @file cheb3d.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la aproximación de chebyshev.
 * @version 0.1
 * @date 2025-04-11
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/cheb3d.h"

Matrix Cheb3D (double t,int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz)
{
    //Comprueba validez
	if ( (t<Ta) || (Tb<t) ) {
		cout<<"Cheb3d: Time out of range::Value\n";
		exit(EXIT_FAILURE);
	}
		
	// Algoritmo de Clenshaw
	double tau = (2*t-Ta-Tb)/(Tb-Ta);  

	Matrix f1(3);
	Matrix f2(3);
    Matrix aux(3);
    Matrix old_f1(3);
	for(int i=N;i>=2;i--){
		old_f1 = f1;
		aux(1)=Cx(i);
        aux(2)=Cy(i);
        aux(3)=Cz(i);
		f1 = f1*2*tau-f2+aux;
		f2 = old_f1;
	}
	Matrix aux2(3);
	aux2(1)=Cx(1);aux2(2)=Cy(1);aux2(3)=Cz(1);
	Matrix ChebApp(3);
	ChebApp = f1*tau-f2+aux2;
	return ChebApp;
}