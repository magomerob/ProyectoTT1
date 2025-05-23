// $Source$
/**
 * @file azelpa.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la función azelpa.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/azelpa.h"

using namespace std;

tuple<double, double, Matrix, Matrix> azelpa(Matrix& s)
{
    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    double Az = atan2(s(1),s(2));

    double pi2 = M_PI*2;

    if (Az<0.0){ 
        Az = Az+pi2;
    }

    double El = atan ( s(3) / rho );

    // Partials
    Matrix dAds = zeros(3);
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;
	
    Matrix dEds = zeros(3);
    dEds(1) = -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    dEds(3) = rho;

    if(s.fil>1){
        s = transpose(s);
    }
    dEds = dEds / dot(s,s);

    return {Az, El, dAds, dEds};
}