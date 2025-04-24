// $Source$
/**
 * @file iers.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la función iers.
 * @version 0.1
 * @date 2025-04-15
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/iers.h"

int find(Matrix m, double v)
{
    for (size_t   i = 1; i <= m.col; i++)
    {
        if (m(i) == v)
            return i;
    }
    return -1;
}

tuple<double, double, double, double, double, double, double, double, double> iers(Matrix& eop, double Mjd_UTC, char interp)
{
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, mjd;

    if (interp =='l'){
        // linear interpolation
        mjd = (floor(Mjd_UTC));

        int i = find(eop.extract_column(4),mjd);

        Matrix preeop(eop.extract_row(i));

        Matrix nexteop(eop.extract_row(i+1));

        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
        y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
        UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
        LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
        dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
        deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
        dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
        dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
        TAI_UTC = preeop(13);
        
        x_pole  = x_pole/SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole  = y_pole/SAT_Const::Arcs;  // Pole coordinate [rad]
        dpsi    = dpsi/SAT_Const::Arcs;
        deps    = deps/SAT_Const::Arcs;
        dx_pole = dx_pole/SAT_Const::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole/SAT_Const::Arcs; // Pole coordinate [rad]
    }
    else if (interp =='n'){ 
        mjd = (floor(Mjd_UTC));
        int i = find(eop.extract_column(4),mjd);
        Matrix eop2 (eop.extract_row(i));
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = eop2(5)/SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole  = eop2(6)/SAT_Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = eop2(7);                   // UT1-UTC time difference [s]
        LOD     = eop2(8);                   // Length of day [s]
        dpsi    = eop2(9)/SAT_Const::Arcs;
        deps    = eop2(10)/SAT_Const::Arcs;
        dx_pole = eop2(11)/SAT_Const::Arcs; // Pole coordinate [rad]
        dy_pole = eop2(12)/SAT_Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = eop2(13);                  // TAI-UTC time difference [s]
    }
    else{
        cout << "IERS: wrong interpolation method\n";
        exit(EXIT_FAILURE);
    }
    return {x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC}; 
}
