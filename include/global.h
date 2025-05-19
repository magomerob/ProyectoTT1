// $Header$
/**
 * @file global.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase global.
 * @date 2025-04-23
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #ifndef _GLOBAL_
 #define _GLOBAL_
 #include <cmath>
 #include "matrix.h"
 #include "mjday.h"
 #include "sat_const.h"
 #include <string.h>
typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Matrix eopdata;
extern Matrix PC;
extern Param AuxParam;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix obs;

void eop19620101(int c);
void DE430Coeff();
void auxparam();
void GGM03S();
void GEOS3();
 #endif
 