/**
 * @file main.cpp
 * @author Marcos Gómez Robres
 * @brief Archivo principal
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "../include/matrix.h"
#include "../include/cheb3d.h"
#include "../include/eccAnom.h"
#include "../include/accelPointMass.h"
#include "../include/legendre.h"
#include <iostream>
#include <tuple>

using namespace std;

int main(){
    cout<<"auto"<<endl;
    auto [p10, p20] = legendre(2, 2, 0);
    cout << p10 << endl;

    cout<<"tie"<<endl;
    Matrix p1(3, 3);
    Matrix p2(3, 3);
    tie(p1, p2) = legendre(2, 2, 0);
    cout<<p1<<endl;
    cout<<"get"<<endl;
    auto res = legendre(2, 2, 0);
    Matrix p11(get<0>(res));
    Matrix p22(get<1>(res));
    cout<<p11<<endl;
    return 0;
}