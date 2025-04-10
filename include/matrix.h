/**
 * @file matrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase Matriz, esta implementa este objeto en C++.
 * @version 0.1
 * @date 2025-04-10
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix
{
    public:
    /**
     * @brief Construye un nuevo objeto Matrix vacío.
     * 
     * @param [in] fil Número de filas.
     * @param [in] col Número de filas.
     */
        Matrix(int fil, int col);

    /**
     * @brief Construye un nuevo objeto Matrix a partir de un vector inicial.
     * 
     * @param [in] fil Número de filas.
     * @param [in] col Número de filas.
     * @param [in] v Vector con los datos de la matriz.
     * @param [in] n Longitud del vector.
     */
        Matrix(int fil, int col, double v[], int n);
    /**
     * @brief Construye un nuevo objeto Matrix a partir de otra matriz.
     * 
     * @param [in] m La otra matriz.
     */
        Matrix(const Matrix& m);
    /**
     * @brief Destroy the Matrix object
     * 
     */
        ~Matrix();
    
    /**
     * @brief Sobreescribe el operador =
     * 
     * @param [in] matrix2 Otra matriz para igualar.
     * @return Matrix& La matriz igual a matrix2.
     */
        Matrix& operator=(const Matrix& matrix2);
    /**
     * @brief Sobreescribe el operador + para la suma.
     * 
     * @param [in] matrix2 La matriz a la que sumar.
     * @return Matrix La matriz resultado de la suma.
     */
        Matrix  operator+(const Matrix& matrix2);
    /**
     * @brief Sobreescribe el operador - para la resta.
     * 
     * @param [in] matrix2 La matriz a la que restar.
     * @return Matrix La matriz resultado de la resta.
     */
        Matrix  operator-(const Matrix& matrix2);
    /**
     * @brief Sobreescribe el operador + para el producto.
     * 
     * @param [in] matrix2 La matriz a la que multiplicar.
     * @return Matrix La matriz resultado del producto.
     */
        Matrix  operator*(const Matrix& matrix2);
    
    /**
     * @brief Sobreescribe el operador() para obtener un elemento de la matriz
     * 
     * @param [in] i fila del elemento.
     * @param [in] j columna del elemento.
     * @return double El elemento en la posición (i,j)
     */
        double& operator()(const int i, const int j) const;
 
 
    private:
        void initMatrix();
 
    public:
        int fil;
        int col;
        double **matrix;
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int fil, const int col);

#endif 
