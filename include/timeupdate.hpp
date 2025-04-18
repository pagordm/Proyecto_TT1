#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "matrix.h"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt);

Matrix& TimeUpdate(Matrix& P, Matrix Phi);

#endif
