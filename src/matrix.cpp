#include "../include/matrix.hpp"

using namespace std;


Matrix::Matrix(int n_row, int n_column){
    if(n_row <= 0 || n_column <= 0){
        cout<<"Matrix: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }

    this->n_row = n_row;
    this->n_column = n_column;
    this->data = (double**)malloc(n_row*sizeof(double*));

    if(this->data == NULL){
        cout<<"Matrix: error in data\n";
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i<n_row; i++){
        this->data[i] = (double*)malloc(n_column*sizeof(double));
    }
}

double& Matrix::operator () (const int n_row, const int n_column){
        if(n_row <= 0 || n_row <= 0 || n_row >= this->n_row ||n_column > this->n_column > n_column){
            cout<<"Matrix get: error in n_row/n_column\n";
            exit(EXIT_FAILURE);
        }
        return this->data[n_row-1][n_column-1];
}

double& Matrix::operator + (Matrix &m){
        if(m.n_row != this->n_row ||this->n_column != m.n_column){
            cout<<"Matrix add: error in n_row/n_column\n";
            exit(EXIT_FAILURE);
        }

        Matrix *m_aux = new Matrix(this->n_row, this->n_column);

        for(int i = 0; i < m.n_row; i++){
            for(int j = 0; j < m.n_column; j++){
                (*m_aux)(i,j) = (*m)(t,j) + (*this)(i,j)
            }

        }
}


ostream& Matrix::operator << (ostream &o, Matrix &m){
    for(int i = 0; i < m.n_row; i++){
        for(int j = 0; j < m.n_column; j++){
            printf("%5.20lf ",m(i,j));
        }
        o << "\n"
    }

    return 0;
}