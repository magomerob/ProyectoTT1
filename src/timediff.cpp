// $Source$
/**
 * @file timediff.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de las diferencias horarias.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/timediff.h"


tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC) {
    const double TT_TAI  =  32.184;   // TT - TAI [s]
    const double GPS_TAI = -19.0;     // GPS - TAI [s]

    double TT_GPS  = TT_TAI - GPS_TAI;  // TT - GPS [s]
    double TAI_GPS = -GPS_TAI;          // TAI - GPS [s]

    double UT1_TAI = UT1_UTC - TAI_UTC;       // UT1 - TAI [s]
    double UTC_TAI = -TAI_UTC;                // UTC - TAI [s]
    double UTC_GPS = UTC_TAI - GPS_TAI;       // UTC - GPS [s]
    double UT1_GPS = UT1_TAI - GPS_TAI;       // UT1 - GPS [s]
    double TT_UTC  = TT_TAI - UTC_TAI;        // TT - UTC [s]
    double GPS_UTC = GPS_TAI - UTC_TAI;       // GPS - UTC [s]

    return { UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC };
}