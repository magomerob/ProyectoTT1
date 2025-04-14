// $Header$
/**
 * @file mjday_tdb.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase mjday, esta implementa un el cálculo de la fecha juliana para tiempo dinámico baricentrico.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _mjdayTdb_
#define _mjdayTdb_
#include <cmath>

using namespace std;

/**
 * @brief Calcula la fecha juliana para tiempo dinámico baricentrico.
 * 
 * @param Mjd_TT Fecha Juliana modificada (Tiempo terrestre).
 * @return double Fecha Juliana modificada (Tiempo dinámico baricéntrico).
 */
double Mjd_TDB(double Mjd_TT);

#endif
