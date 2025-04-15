// $Source$
/**
 * @file matrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de una matriz en C++
 * @version 0.1
 * @date 2025-04-10
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.fil; i++) {
        for (int j = 1; j <= m.col; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
} 

void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(this->fil, this->col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}
 
 
double& Matrix::operator()(const int i, const int j) const
{
    if (i <= 0 || i > this->fil || j <= 0 || j > this->col) {
		cout << "Matrix get: position out of bounds\n";
        exit(EXIT_FAILURE);
	}
    return matrix[i-1][j-1];
}

double& Matrix::operator()(const int n) const
{
    if (n <= 0 || n > this->col) {
		cout << "Vector get: position out of bounds\n";
        exit(EXIT_FAILURE);
	}
    return matrix[0][n-1];
}

Matrix Matrix::operator/(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    Matrix m2inv = inv(matrix2);

    result = (*this)*m2inv;
 
    return result;
}

Matrix inv(const Matrix& matrix)
{
    if (matrix.col != matrix.fil) {
		cout << "Matrix inv: not square matrix\n";
        exit(EXIT_FAILURE);
	}

    //Método de Gauss-Jordan
    Matrix extendida(matrix.fil, matrix.col*2);

    //LLenamos el lado izquierdo con la matriz original
    for (size_t i = 1; i <= matrix.fil; ++i) {
        for (size_t j = 1; j <= matrix.col; ++j) {
            extendida(i, j) = matrix(i, j);
        }
    }
    
    //LLenamos el lado izquierdo con la matriz identidad
    for (size_t i = 1; i <= matrix.fil; ++i) {
            extendida(i, i + matrix.col) = 1.0;
    }
    //Aplicamos el método
    for (size_t i = 1; i <= matrix.fil; ++i) {
        // Busca pivote
        size_t fila_pivote = i;
        double max_val = abs(extendida(i, i));
        
        // Encuentra la fila con el mayor valor en la columna
        for (size_t j = i + 1; j <= matrix.fil; ++j) {
            if (abs(extendida(j, i)) > max_val) {
                max_val = abs(extendida(j, i));
                fila_pivote = j;
            }
        }
        // Comprobar si la matriz es singular
        if (max_val < 1e-10) {
            cout << "Matrix inv: singular matrix\n";
            exit(EXIT_FAILURE);
        }
        // Cambiar filas si es necesario
        if (fila_pivote != i) {
            Matrix _f1 = extendida.extract_row(fila_pivote);
            Matrix _f2 = extendida.extract_row(i);
            extendida.assign_row(_f1, i);
            extendida.assign_row(_f2, fila_pivote);
        }

        // Escalar la fila
        double pivote = extendida(i, i);
        for (size_t j = 1; j <= 2 * matrix.col; ++j) {
            extendida(i, j) /= pivote;
        }

        // Eliminar otras filas
        for (size_t j = 1; j <= matrix.fil; ++j) {
            if (j != i) {
                double factor = extendida(j, i);
                for (size_t k = 1; k <= 2 * matrix.col; ++k) {
                    extendida(j, k) -= factor * extendida(i, k);
                }
            }
        }
    }

    // El resultado es la parte derecha de la matriz extendida
    Matrix result(matrix.fil, matrix.col);
    for (size_t i = 1; i <= matrix.fil; ++i) {
        for (size_t j = 1; j <= matrix.col; ++j) {
            result(i, j) = extendida(i, j + matrix.col);
        }
    }
    
    return result;
}

Matrix Matrix::extract_row(const int row)
{
    if (row <= 0 || row > this->fil) {
		cout << "Row get: position out of bounds\n";
        exit(EXIT_FAILURE);
	}

    Matrix res(1,this->col);
    for(size_t i=1; i<=this->col; i++){
        res(i)=(*this)(row, i);
    }
    return res;
}

Matrix Matrix::extract_column(const int column)
{
    if (column <= 0 || column > this->fil) {
		cout << "Column get: position out of bounds\n";
        exit(EXIT_FAILURE);
	}

    Matrix res(1,this->fil);
    for(size_t i=1; i<=this->fil; i++){
        res(i)=(*this)(i, column);
    }
    return res;
}

void Matrix::assign_row(const Matrix& fil, const int ind)
{
    if (ind <= 0 || ind > this->fil) {
		cout << "Row set: position out of bounds\n";
        exit(EXIT_FAILURE);
	}
    if(fil.col != this->col || fil.fil != 1){
        cout << "Row set: wrong size vector\n";
        exit(EXIT_FAILURE);
    }

    for(size_t i=1; i<=this->col; i++){
        (*this)(ind, i) = fil(i);
    }
}

void Matrix::assign_column(const Matrix& col, const int ind)
{
    if (ind <= 0 || ind > this->col) {
		cout << "Column set: position out of bounds\n";
        exit(EXIT_FAILURE);
	}
    if(col.col != this->fil || col.fil != 1){
        cout << "Column set: wrong size vector\n";
        exit(EXIT_FAILURE);
    }

    for(size_t i=1; i<=this->fil; i++){
        (*this)(i, ind) = col(i);
    }
}

Matrix& eye(const int fil, const int col) {
	Matrix *m_aux = new Matrix(fil, col);
	
	for(size_t i = 1; i <= fil; i++) {
		for(size_t j = 1; j <= col; j++) {
			(*m_aux)(i,j) = 1;
		}
	}
	
	return (*m_aux);
} 

Matrix& transpose(const Matrix& matrix)
{
    Matrix *m_aux = new Matrix(matrix.col, matrix.fil);
	
	for(size_t i = 1; i <= matrix.fil; i++) {
		for(size_t j = 1; j <= matrix.col; j++) {
			(*m_aux)(j,i) = matrix(i,j);
		}
	}
	
	return (*m_aux);
}

Matrix Matrix::operator+(const double n)
{
    Matrix result(fil, col);
    
    for (size_t i = 0; i < fil; i++)
        for (size_t j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + n;
 
    return result;
}
 
Matrix Matrix::operator-(const double n)
{
    Matrix result(fil, col);
    
    for (size_t i = 0; i < fil; i++)
        for (size_t j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - n;
 
    return result;
}

Matrix Matrix::operator*(const double n)
{
    Matrix result(fil, col);
    
    for (size_t i = 0; i < fil; i++)
        for (size_t j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] * n;
 
    return result;
}
 
Matrix Matrix::operator/(const double n)
{
    if(abs(n)<1e-10){
        cout << "Div: can't divide by 0\n";
        exit(EXIT_FAILURE);
    }

    Matrix result(fil, col);
    
    for (size_t i = 0; i < fil; i++)
        for (size_t j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] / n;
 
    return result;
}

Matrix::Matrix(int n) : fil(1), col(n)
{
    initMatrix();
}

double norm(Matrix& vec)
{
    if (vec.fil != 1){
		cout << "Vector norm: not a vector\n";
        exit(EXIT_FAILURE);
	}

    double sum = 0.0;
	for(size_t i = 1; i <=vec.col; i++) {
		sum+=pow(vec(i), 2);
	}
	return sqrt(sum);
}

double dot(Matrix &vec1, Matrix &vec2) {
	if (vec1.fil!=1 || vec2.fil!=1) {
		cout << "Vector dot: not a vector" << endl;
		exit(EXIT_FAILURE);
	}
    if(vec1.col!=vec2.col){
        cout << "Vector dot: different size vectors" << endl;
		exit(EXIT_FAILURE);
    }
	double sum = 0.0;
	for(size_t i = 1; i <=vec1.col; i++) {
		sum+=vec1(i)*vec2(i);
	}

	return sum;
}

Matrix& cross(Matrix &vec1,Matrix &vec2){
	if(vec1.col!=3 || vec2.col!=3 || vec1.fil!=1 || vec2.fil!=1){
		cout<<"Vector cross: not size 3 vectors\n";
		exit(EXIT_FAILURE);
	}
	Matrix *m_aux=new Matrix(3);
	(*m_aux)(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2);
	(*m_aux)(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3);
	(*m_aux)(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1);
	
	return *m_aux;
}

Matrix extract_vector (Matrix &v,int i, int j)
{
    if(v.fil!=1){
		cout<<"Vector extract: not a vector\n";
		exit(EXIT_FAILURE);
	}
    if(v.col < j || i<1 || i>j){
        cout<<"Vector extract: wrong indexes\n";
		exit(EXIT_FAILURE);
    }
    Matrix res(j-i+1);
    for(size_t k = i; k<=j; k++){
        res(k-i+1) = v(k);
    }

    return res;
}

Matrix union_vector(Matrix &v1, Matrix &v2)
{
    if(v1.fil!=1 || v2.fil!=1){
		cout<<"Vector union: not a vector\n";
		exit(EXIT_FAILURE);
	}
    Matrix res(v1.col+v2.col);
    for(size_t i = 1; i<=v1.col; i++){
        res(i) = v1(i);
    }
    for(size_t i = 1; i<=v2.col; i++){
        res(i+v1.col) = v2(i);
    }
    return res;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {

    if (this == &other) return *this;

    // Free old memory
    for (int i = 0; i < fil; ++i)
        delete[] matrix[i];
    delete[] matrix;

    // Steal data
    matrix = other.matrix;
    fil = other.fil;
    col = other.col;

    other.matrix = nullptr;
    other.fil = 0;
    other.col = 0;

    return *this;
}