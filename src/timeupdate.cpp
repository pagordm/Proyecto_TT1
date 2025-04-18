#include "../include/timeupdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt) {
    
    P = Phi * P * Phi.transpose() + Qdt;
    return P;
}

Matrix& TimeUpdate(Matrix& P, Matrix Phi) {
    return TimeUpdate(P, Phi, 0.0);
}