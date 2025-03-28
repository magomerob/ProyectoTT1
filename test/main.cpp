#include "../include/matrix.hpp"
#include <iostream>

using namespace std;

int main(){
    Matrix M1(3,2);

    M1(1,2) = 5.0;

    cout << M1;

    return 0;
}