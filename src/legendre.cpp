// $Source$
/**
 * @file legendre.cpp
 * @author Marcosf Gómez Robres
 * @brief Implementación de la función legendre.
 * @version 0.1
 * @date 2025-04-15
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/legendre.h"

tuple<Matrix, Matrix> legendre(int n, int m, double fi)
{
    // Inicialización de las matrices
    Matrix pnm = zeros(n+1, m+1);
    Matrix dpnm = zeros(n+1, m+1);

    // Inicialización de los valores iniciales
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2,2)=sqrt(3.0)*cosf(fi);
    dpnm(2,2)=-sqrt(3.0)*sinf(fi);
    
    for(size_t i = 2; i<= n; i++){
        pnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*cosf(fi)*pnm(i,i);
        dpnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*((cosf(fi)*dpnm(i,i))-(sinf(fi)*pnm(i,i)));
    }
    // horizontal first step coefficients
    for(size_t i = 1; i<= n; i++){ 
        pnm(i+1,i)= sqrt(2.0*i+1)*sinf(fi)*pnm(i,i);
    
        dpnm(i+1,i)= sqrt(2.0*i+1)*((cosf(fi)*pnm(i,i))+(sinf(fi)*dpnm(i,i)));
    }
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(1){
        for(size_t i = k; i<= n; i++){         
            pnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)*1.0))*((sqrt(2.0*i-1)*sinf(fi)*pnm(i,j+1))
            -(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*pnm(i-1,j+1)));
            dpnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)*1.0))*((sqrt(2.0*i-1.0)*sinf(fi)*dpnm(i,j+1))
                +(sqrt(2.0*i-1.0)*cosf(fi)*pnm(i,j+1))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
    return make_tuple(pnm, dpnm);
}