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
    double m = 1;
    Matrix Inv_W = zeros(m,m);

    Inv_W(1,1) = s*s;    // Inverse weight (measurement covariance)
    //cout<<"P: "<<endl<<P<<endl;
    //cout<<"G: "<<endl<<G<<endl;
    //cout<<"Inv_W: "<<endl<<Inv_W<<endl;
    // Kalman gain
    Matrix K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));
	//cout<<"K: "<<endl<<K<<endl;
    // State update
    Matrix x2 = x + K*(z-g);
    //cout<<"x2: "<<endl<<K<<endl;
    //Covariance update
    Matrix P2 = (eye(n,n)-K*G)*P;
    //cout<<"p2: "<<endl<<K<<endl;
    return tie(K, x2, P2);
}
