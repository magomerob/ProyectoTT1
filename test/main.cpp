#include "../include/matrix.h"
#include <iostream>

using namespace std;

int main(){
    Matrix M1(3,2);
    Matrix M2(3,2);

    M1(1,2) = 5.0;
    M2(1,1) = 2.0;
    Matrix M3 = (M1+M2);
    cout << M3;

    return 0;
}