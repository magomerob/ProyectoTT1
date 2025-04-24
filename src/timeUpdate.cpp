// $Source$
/**
 * @file timeUpdate.cpp
 * @author Marcos Gómez Robres
 * @brief Esta es la implementación de la función TimeUpdate.
 * @version 0.1
 * @date 2025-04-24
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "..\include\timeUpdate.h"

Matrix& TimeUpdate (Matrix &P, Matrix &Phi, double Qdt){
	return Phi*P*transpose(Phi) + Qdt;
	
}
Matrix& TimeUpdate (Matrix &P, Matrix &Phi){
	return TimeUpdate(P,Phi,0.0);
}