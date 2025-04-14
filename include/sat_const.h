// $Header$
/**
 * @file sat_const.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase SAT_Const, esta recoge las constantes usadas en el resto de programas.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _SAT_Const_
#define _SAT_Const_
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * @brief Lista de constantes.
 * 
 */
class SAT_Const {
public:
    double pi2,Rad,Deg,Arcs,MJD_J2000,T_B1950,c_light,AU,R_Earth,f_Earth,R_Sun,R_Moon,omega_Earth,GM_Earth,GM_Sun,GM_Moon,GM_Mercury,GM_Venus,GM_Mars,GM_Jupiter,GM_Saturn,GM_Uranus,GM_Neptune,GM_Pluto,P_Sol;

    /**
     * @brief Constructor del objeto SAT_const
     * 
     */
    SAT_Const();
};
#endif