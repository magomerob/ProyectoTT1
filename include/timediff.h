// $Header$
/**
 * @file timediff.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase timediff, esta implementa una lista de diferencias horarias.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _timediff_
#define _timediff_

/**
 * @brief Lista con los diferentes usos horarios.
 * 
 */
struct TimeDifferences {
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
};

/**
 * @brief Función que devuelve las diferencias horarias.
 * 
 * @param UT1_UTC 
 * @param TAI_UTC 
 * @return TimeDifferences Diccionario con los resultados.
 */
TimeDifferences timediff(double UT1_UTC, double TAI_UTC);
#endif
