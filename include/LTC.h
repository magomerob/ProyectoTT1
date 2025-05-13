// $Header$
/**
 * @file LTC.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase poleMatrix, esta calcula la Transformación del sistema del meridiano de Greenwich a coordenadas tangentes locales
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _LTC_
#define _LTC_
#include "matrix.h"
#include "r_y.h"
#include "r_z.h"

/**
 * @brief calcula la transformación del sistema del meridiano de Greenwich a coordenadas tangentes locales
 * 
 * @param lon longitud
 * @param lat latitud
 * @return Matrix matriz de rotación
 */
Matrix LTC(double lon, double lat);

#endif