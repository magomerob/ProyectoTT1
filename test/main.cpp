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
#include <iostream>
#include <tuple>

using namespace std;

int main(){

    Matrix v(3);

    v(1) = 1;
    v(2) = 2;
    v(3) = 3;

    //auto [a,b,c,d] = azelpa(v);

    //cout << a <<endl<< b<<endl << c << d << endl;

    auto [A,B] = legendre(2, 2, 1);
    cout << A <<endl<< B << endl;

    /*
    resultado esperado
    1.00000000000000000000 0.00000000000000000000 0.00000000000000000000 
    1.45747045027529753547 0.93583099453595375294 0.00000000000000000000
    1.25691629764617052167 1.76084674147066455596 0.56531333344858780698 

    0.00000000000000000000 0.00000000000000000000 0.00000000000000000000
    0.93583099453595375294 -1.45747045027529753547 0.00000000000000000000
    3.04987602056929052452 -1.61172970742831900282 -1.76084674147066455596
    */

    return 0;
}
