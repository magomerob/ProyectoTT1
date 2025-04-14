// $Header$
/**
 * @file mjday.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase mjday, esta implementa un el cálculo de la fecha juliana.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _mjday_
#define _mjday_
#include <cmath>

using namespace std;

/**
 * @brief Calcula la fecha juliana.
 * 
 * @param yr Año.
 * @param mon Mes.
 * @param day Día.
 * @param hr Hora.
 * @param min Minuto.
 * @param sec Seguntos
 * @return double Fecha Juliana modificada.
 */
double Mjday(int yr, int mon, int day, int hr, int min, double sec);

#endif
