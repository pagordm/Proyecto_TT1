#include "../include/timeupdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt) {
    
    Matrix& nP = Phi * P * Phi.transpose() + Qdt;
    return nP;
}

Matrix& TimeUpdate(Matrix& P, Matrix Phi) {
    return TimeUpdate(P, Phi, 0);
}
