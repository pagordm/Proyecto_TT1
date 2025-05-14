#include "..\include\unit.hpp"

/**
 * @brief Calculates a unit vector given the original vector
 * 
 * If a zero vector is input, the vector is set to zero.
 * 
 * @param vec Vector to normalize
 * @return Matrix& Unit vector output
 */

Matrix& unit(Matrix vec) {
    double small = 0.000001;
    double magv = norm(vec.transpose());
    Matrix* outvec = new Matrix(3);
    if ( magv > small ) {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= vec(i)/magv;
        }
    } else {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= 0.0;
        }
    }

    return *outvec;

}