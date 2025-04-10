#include "../include/matrix.h"
#include <iostream>

using namespace std;

int main(){
    Matrix M1(3,2);
    Matrix M2(3,2);

    M1(1,2) = 5.0;
    M2(1,1) = 5.0;
    Matrix M3 = (M1+M2);
    cout << M3;

    cout<<"---------------"<<endl;

    Matrix M4(2,2);
    M4(1,1)=2;
    M4(1,2)=7;
    M4(2,1)=2;
    M4(2,2)=8;
    cout<<M4;
    M4 = inv(M4);
    cout<<"---------------"<<endl;
    cout << M4;
    return 0;
}