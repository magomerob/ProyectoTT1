// $Source$
/**
 * @file measUpdate.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de measUpdate.
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/measUpdate.h"

using namespace std;

tuple<Matrix, Matrix, Matrix> measUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n)
{
    int m = 1;
    Matrix Inv_W = zeros(m,m);

    Inv_W(1,1) = s*s;// Inverse weight (measurement covariance)
    

    // Kalman gain
    Matrix K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));

    // State update
    Matrix newx = x + K*(z-g);

    // Covariance update
    Matrix newP = (eye(n,n)-K*G)*P;

    return tie(K, newx, newP);
}
