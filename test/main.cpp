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
#include "../include/azelpa.h"
#include "../include/eccAnom.h"
#include "../include/accelPointMass.h"
#include "../include/legendre.h"
#include "../include/nutAngles.h"
#include "../include/accelHarmonic.h"
#include "../include/JPL_Eph.h"
#include "../include/global.h"
#include <iostream>
#include <tuple>

using namespace std;

int main(){

    DE430Coeff();
   
    auto [a, b, c, d, e, f, g, h, i, j, k] = JPL_Eph(50000);
    
    return 0;
}
