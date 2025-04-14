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
#include <iostream>

using namespace std;

int main(){
    double M = 2;
    double e = 0.6;

    double r = eccAnom(M, e);
    cout<<r;
    return 0;
}