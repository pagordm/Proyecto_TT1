/**
 * @file global.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains global variables used throughout the project, and functions to initialize them.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _GLOBAL_
#define _GLOBAL_

#include "matrix.h"
#include <cmath>
#include <string.h>
#include "mjday.hpp"
#include "const.hpp"

typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Matrix eopdata;

extern Matrix Cnm;
extern Matrix Snm;

extern Matrix PC;
extern Param AuxParam;
extern Matrix obs;

void eop19620101(int c);

void GGM03S();

void DE430Coeff();

void auxparam();

void GEOS3(int f);

#endif