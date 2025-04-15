// $Source$
/**
 * @file legendre.cpp
 * @author Marcos Gómez Robres
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
    Matrix pnm(n+1, m+1);
    Matrix dpnm(n+1, m+1);

    // Inicialización de los valores iniciales
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2,2)=sqrt(3)*cos(fi);
    dpnm(2,2)=-sqrt(3)*sin(fi);
    

    for(size_t i = 2; i<= n; i++){
        pnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*cos(fi)*pnm(i,i);
        dpnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*((cos(fi)*dpnm(i,i))-(sin(fi)*pnm(i,i)));
    }
    // horizontal first step coefficients
    for(size_t i = 1; i<= n; i++){ 
        pnm(i+1,i)= sqrt(2*i+1)*sin(fi)*pnm(i,i);
    
        dpnm(i+1,i)= sqrt(2*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
    }
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(1){
        for(size_t i = k; i<= n; i++){         
            pnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*pnm(i,j+1))
                -(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-1,j+1)));
            
            dpnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*dpnm(i,j+1))
                +(sqrt(2*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
    
    //cout<< pnm << endl;
    cout<<"---"<<endl;
    return {pnm, dpnm};
}