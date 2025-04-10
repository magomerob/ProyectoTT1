// $Header$
/**
 * @file matrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase Matriz, esta implementa este objeto en C++.
 * @version 0.1
 * @date 2025-04-10
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
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
        Matrix operator+(const Matrix& matrix2);
    
    /**
     * @brief Sobreescribe el operador - para la resta.
     * 
     * @param [in] matrix2 La matriz a la que restar.
     * @return Matrix La matriz resultado de la resta.
     */
        Matrix operator-(const Matrix& matrix2);
    
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
    
    /**
     * @brief Sobreescribe el operador() para obtener un elemento de un vector
     * 
     * @param [in] n posición del elemento.
     * @return double El elemento en la posición n
     */
        double& operator()(const int n) const;
    
    /**
     * @brief Sobreescribe el operador / para la división.
     * 
     * @param [in] matrix2 La matriz entre la que dividir.
     * @return Matrix La matriz resultado de la división.
     */
        Matrix  operator/(const Matrix& matrix2);
    
    /**
     * @brief función para extraer una fila.
     * 
     * @param [in] row el índice de la fila a extraer.
     * @return Matrix Fila extraida en forma de vector.
     */
        Matrix extract_row(const int row);
    
    /**
     * @brief función para extraer una columna.
     * 
     * @param [in] column el índice de la columna a extraer.
     * @return Matrix Columna extraida en forma de vector.
     */
        Matrix extract_column(const int column);
    
    /**
     * @brief función para sustituir una columna.
     * 
     * @param [in] col La columna a insertar.
     * @param [in] ind el índice de la columna a sustituir.
     */
        void assign_column(const Matrix& col, const int ind);
    
    /**
     * @brief función para sustituir una fila.
     * 
     * @param [in] fil La fila a insertar.
     * @param [in] ind el índice de la fila a sustituir.
     */
        void assign_row(const Matrix& fil, const int ind);
    private:
        void initMatrix();
 
    public:
        int fil;
        int col;
        double **matrix;
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

/**
 * @brief Genera una matriz de ceros con tamaño fil*col.
 * 
 * @param fil Numero de filas de la matriz.
 * @param col Numero de columnas de la matriz.
 * @return Matrix& matriz de ceros.
 */
Matrix& zeros(const int fil, const int col);

/**
 * @brief función para calcular la matriz inversa usando el método de Gauss-Jordan
 * 
 * @param [in] matrix La matriz a invertir.
 * @return Matrix La matriz inversa.
 */
Matrix inv(const Matrix& matrix);

#endif 
