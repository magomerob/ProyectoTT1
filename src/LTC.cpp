// $Source$
/**
 * @file nutMatrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de LTC.
 * @version 0.1
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/LTC.h"

Matrix LTC(double lon, double lat)
{
    Matrix M = r_y(-1.0*lat) * r_z(lon);

    double aux;
    for(size_t i = 1; i <= 3; i++){
        aux=M(1,i);
        M(1,i)=M(2,i);
        M(2,i)=M(3,i);
        M(3,i)=aux;
    }

    return M;  
}