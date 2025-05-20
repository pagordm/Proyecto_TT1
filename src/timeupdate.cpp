/**
 * @file timeupdate.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function timeupdate.
 * @date 2025-05-20
 * 
 * 
 */
#include "../include/timeupdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt) {
    
    Matrix& nP = Phi * P * Phi.transpose() + Qdt;
    return nP;
}

Matrix& TimeUpdate(Matrix& P, Matrix Phi) {
    return TimeUpdate(P, Phi, 0);
}
