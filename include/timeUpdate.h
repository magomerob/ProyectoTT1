// $Header$
/**
 * @file timeUpdate.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la función timeUpdate, esta actualiza el tiempo.
 * @version 0.1
 * @date 2025-04-24
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include <cmath>
#include "matrix.h"


/**
 * @brief Actualiza el tiempo.
 * 
 *
 * @param [in] P Matriz.
 * @param [in] Phi Matriz.  
 * @param [in] Qdt Double.  
 * @return Matrix tiempo actualizado.
 */
Matrix& TimeUpdate (Matrix &P, Matrix &Phi, double Qdt);	

/**
 * @brief Actualiza el tiempo, da el valoR 0 a Qdt.
 * 
 *
 * @param [in] P Matriz.
 * @param [in] Phi Matriz.  
 * @return Matrix tiempo actualizado. 
 */
Matrix& TimeUpdate (Matrix &P, Matrix &Phi);	
#endif