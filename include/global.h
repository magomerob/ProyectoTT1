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
 

extern Matrix eopdata;
extern Matrix PC;

void eop19620101(int c);
void DE430Coeff();
 #endif
 