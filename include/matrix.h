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
     * @brief Construye un nuevo objeto Matrix del tipo vector vacío.
     * 
     * @param [in] n Número de elementos.
     */
        Matrix(int n);

    /**
     * @brief Construye un nuevo objeto Matrix vacío.
     * 
     */
        Matrix();

    /**
     * @brief Construye un nuevo objeto Matrix a partir de un vector inicial.
     * 
     * @param [in] fil Número de filas.
     * @param [in] col Número de filas.
     * @param [in] v Vector con los datos de la matriz.
     * @param [in] n Longitud del vector.
     */
        Matrix(int fil, int col, double v[], int n); //skipcq: CXX-W2066
    
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
    
    /**
     * @brief Sobreescribe el operador + para la suma.
     * 
     * @param [in] n El real que sumar.
     * @return Matrix La matriz resultado de la suma.
     */
        Matrix operator+(const double n);
    
    /**
     * @brief Sobreescribe el operador - para la resta.
     * 
     * @param [in] n El real que restar.
     * @return Matrix La matriz resultado de la resta.
     */
        Matrix operator-(const double n);
    
    /**
     * @brief Sobreescribe el operador + para la multiplicación.
     * 
     * @param [in] n El real que multiplicar.
     * @return Matrix La matriz resultado del producto.
     */
        Matrix operator*(const double n);
    
    /**
     * @brief Sobreescribe el operador / para la división.
     * 
     * @param [in] n El real que divide.
     * @return Matrix La matriz resultado de la división.
     */
        Matrix operator/(const double n);
    
        Matrix& operator=(Matrix&& other) noexcept;

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
 * @brief Genera una vector de ceros con tamaño n.
 * 
 * @param n Numero de filas de la matriz.
 * @return Matrix& matriz de ceros.
 */
Matrix& zeros(const int n);

/**
 * @brief función para calcular la matriz inversa usando el método de Gauss-Jordan
 * 
 * @param [in] matrix La matriz a invertir.
 * @return Matrix La matriz inversa.
 */
Matrix inv(const Matrix& matrix);

/**
 * @brief Genera una matriz de unos con tamaño fil*col.
 * 
 * @param fil Numero de filas de la matriz.
 * @param col Numero de columnas de la matriz.
 * @return Matrix& matriz de ceros.
 */
Matrix& eye(const int fil, const int col);

/**
 * @brief Función para transponer una matriz.
 * 
 * @param [in] matrix La matriz a transponer.
 * @return Matrix La matriz transpuesta.
 */
Matrix& transpose(const Matrix& matrix);

/**
 * @brief Función para calcualr la norma de un vector.
 * 
 * @param vec Vector al que se le calcula la norma.
 * @return double Norma del vector.
 */
double norm(Matrix& vec);

/**
 * @brief Función para calcualr el producto escalar de dos vectores.
 * 
 * @param vec1 Vector 1.
 * @param vec2 Vector 2.
 * @return double Norma del vector.
 */
double dot(Matrix& vec1, Matrix& vec2);

/**
 * @brief Función para calcualr el producto vectorial de dos vectores.
 * 
 * @param vec1 Vector 1.
 * @param vec2 Vector 2.
 * @return Matrix& Vector resultado.
 */
Matrix& cross(Matrix& vec1, Matrix& vec2);

/**
 * @brief Extrae un subvector de un vector.
 * 
 * @param v Vector del cual se extrae.
 * @param i Índice a partir del cual se extrae.
 * @param j Índice hasta el que se extrae.
 * @return Matrix& Subvector resultado.
 */
Matrix extract_vector (Matrix &v,int i, int j);

/**
 * @brief Une dos vectores.
 * 
 * @param v1 Primer vector.
 * @param v2 Segundo vector
 * @return Matrix Vector resultado de la concatenación.
 */
Matrix union_vector(Matrix &v1, Matrix &v2);
#endif 
