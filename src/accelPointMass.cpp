// $Source$
/**
 * @file accelPointMass.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la aceleración respecto a un punto de masa.
 * @version 0.1
 * @date 2025-04-11
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/accelPointMass.h"

Matrix AccelPointMass (Matrix &r, Matrix &s, double GM)
{
    //Posición relativa
    Matrix d(3);
    d = r-s;

    //Aceleración
    Matrix a(3);
    a = ((d/(pow(norm(d),3))) + (s/(pow(norm(s),3))))*(-GM);
    return a;
}