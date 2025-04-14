// $Header$
/**
 * @file position.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase mjday, esta implementa un el cálculo de la posición en coordenadas geodésicas.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _position_
#define _position_
#include <cmath>
#include "matrix.h"
using namespace std;

/**
 * @brief Calcula la posición en coordenadas geodésicas.
 * 
 * @param lon Longitud en radianes.
 * @param lat Latitud en radianes.
 * @param h Altura en metros.
 * @return Matrix Vector con las coordenadas geodésicas.
 */
Matrix Position(double lon, double lat, double h);

#endif
