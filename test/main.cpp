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
#include <iostream>
#include <tuple>

using namespace std;

int main(){

    Matrix v(3,3);

    v(1,1) = 1;
    v(2,1) = 4;
    v(3,1) = 7;
    v(1,2) = 2;
    v(2,2) = 5;
    v(3,2) = 8;
    v(1,3) = 3;
    v(2,3) = 6;
    v(3,3) = 9;
    cout<<v<<endl;

    Matrix E(3);
    E(1) = 10;
    E(2) = 10;
    E(3) = 10;

    Matrix E2 = transpose(E);
    //auto [a,b,c,d] = azelpa(v);

    //cout << a <<endl<< b<<endl << c << d << endl;

    /*auto [A,B] = legendre(2, 2, 1);
    cout << A <<endl<< B << endl;

    /*
    resultado esperado
    1.00000000000000000000 0.00000000000000000000 0.00000000000000000000 
    1.45747045027529753547 0.93583099453595375294 0.00000000000000000000
    1.25691629764617052167 1.76084674147066455596 0.56531333344858780698 

    0.00000000000000000000 0.00000000000000000000 0.00000000000000000000
    0.93583099453595375294 -1.45747045027529753547 0.00000000000000000000
    3.04987602056929052452 -1.61172970742831900282 -1.76084674147066455596
    

    auto [f,g] = NutAngles(2003);
    cout << f <<endl<< g << endl;
    */
    Matrix m = accelHarmonic(E2,v, 1,10);

    /*
    -38518152281.9955
    -45925489259.3023
    -53332826236.6091
    */
    cout<<m<<endl;
    return 0;
}
