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

#ifndef _SAT_CONST_
#define _SAT_CONST_
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <cmath>

namespace SAT_Const {
    // Mathematical constants
    /**
     * @brief Constante 2*pi
     */
    constexpr double pi2       = M_PI * 2;                 // 2π
    /**
     * @brief Constante radianes por grado
     */
    constexpr double Rad       = M_PI / 180;               // Radians per degree
    /**
     * @brief Constante grados por radian
     */
    constexpr double Deg       = 180 / M_PI;               // Degrees per radian
    /**
     * @brief Constante arcosegundos por radian
     */
    constexpr double Arcs      = 3600 * 180 / M_PI;        // Arcseconds per radian

    // General
    /**
     * @brief Constante fecha juliana de J2000
     */
    constexpr double MJD_J2000 = 51544.5;                  // Modified Julian Date of J2000
    /**
     * @brief Constante época B1950
     */
    constexpr double T_B1950   = -0.500002108;             // Epoch B1950
    /**
     * @brief Constante velocidad de la luz
     */
    constexpr double c_light   = 299792458.0;              // Speed of light [m/s]
    /**
     * @brief Constante unidad astronómica
     */
    constexpr double AU        = 149597870700.0;           // Astronomical unit [m]

    // Equatorial radius and flattening
    /**
     * @brief Constante radio de la Tierra
     */
    constexpr double R_Earth   = 6378.1363e3;              // Earth's radius [m]
    /**
     * @brief Constante aplanamiento de la Tierra
     */
    constexpr double f_Earth   = 1.0 / 298.257223563;      // Flattening
    /**
     * @brief Constante radio del Sol
     */
    constexpr double R_Sun     = 696000e3;                 // Sun's radius [m]
    /**
     * @brief Constante radio de la Luna
     */
    constexpr double R_Moon    = 1738e3;                   // Moon's radius [m]

    // Earth rotation
    /**
     * @brief Constante velocidad angular de la Tierra
     */
    constexpr double omega_Earth = 15.04106717866910 / 3600 * Rad; // [rad/s]

    // Gravitational constants
    /**
     * @brief Constante gravitacional de la Tierra
     */
    constexpr double GM_Earth    = 398600.435436e9;                           // [m^3/s^2]
    /**
     * @brief Constante gravitacional del Sol
     */
    constexpr double GM_Sun      = 132712440041.939400e9;                     // [m^3/s^2]
    /**
     * @brief Constante gravitacional de la Luna
     */
    constexpr double GM_Moon     = GM_Earth / 81.30056907419062;             // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Mercurio
     */
    constexpr double GM_Mercury  = 22031.780000e9;                            // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Venus
     */
    constexpr double GM_Venus    = 324858.592000e9;                           // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Marte
     */
    constexpr double GM_Mars     = 42828.375214e9;                            // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Júpiter
     */
    constexpr double GM_Jupiter  = 126712764.800000e9;                        // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Saturno
     */
    constexpr double GM_Saturn   = 37940585.200000e9;                         // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Urano
     */
    constexpr double GM_Uranus   = 5794548.600000e9;                          // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Neptuno
     */
    constexpr double GM_Neptune  = 6836527.100580e9;                          // [m^3/s^2]
    /**
     * @brief Constante gravitacional de Plutón
     */
    constexpr double GM_Pluto    = 977.0000000000009e9;                       // [m^3/s^2]

    // Solar radiation pressure
    /**
     * @brief Constante presión de radiación solar
     */
    constexpr double P_Sol       = 1367 / c_light;           // [N/m^2]
}

#endif