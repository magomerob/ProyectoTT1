// $Header$
/**
 * @file g_accelHarmonic.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase g_accelHarmonic, esta calcula el gradiente del campo gravitacional armónico de la Tierra
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _G_ACCELHARMONIC_
#define _G_ACCELHARMONIC_
#include "matrix.h"
#include "accelHarmonic.h"

/**
 * @brief Esta es la cabecera de la clase g_accelHarmonic, esta calcula el gradiente del campo gravitacional armónico de la Tierra
 * 
 * @param r posición del satelite
 * @param u matriz de transformación
 * @param n grado del modelo de gravedad
 * @param m orden del modelo de gravedad
 * @return Matrix matriz de rotación
 */

Matrix g_accelHarmonic(Matrix r, Matrix u, int n, int m);

#endif