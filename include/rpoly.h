/**
 * @file rpoly.h
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the code for calculating roots of polynomials.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef BEGIN_C_DECLS
#define BEGIN_C_DECLS

// Opaque state handle
struct RPoly_State;

// Allocate state for finding roots of a real polynomial.
// Degrees up to max_degree are supported by returned state.
struct RPoly_State *real_poly_alloc(int max_degree);

// Release state
void real_poly_release(struct RPoly_State *s);


// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Caller must have already allocated 'state'.
// Returns number of roots stored in zeror[], zeroi[].
int real_poly_roots_compute(const double p[], int degree,
                            struct RPoly_State *state,
                            double zeror[], double zeroi[]);

/**
 * @brief Convenience function to find roots of a polynomial with real coefficients.
 *
 * This function computes the roots of a polynomial with real coefficients.
 * The coefficients are provided in the array \p p, where:
 *   p(x) = p[0]*x^degree + ... + p[degree]
 *
 * @param[in]  p       Array of polynomial coefficients.
 * @param[in]  degree  Degree of the polynomial.
 * @param[out] zeror   Array to store the real parts of the roots.
 * @param[out] zeroi   Array to store the imaginary parts of the roots.
 * @return Number of roots stored in \p zeror and \p zeroi.
 */
static inline int real_poly_roots(const double p[], int degree,
                                  double zeror[], double zeroi[])
{
    struct RPoly_State *state = real_poly_alloc(degree);
    int nr = real_poly_roots_compute(p, degree, state, zeror, zeroi);
    real_poly_release(state);

    return nr;
}

#endif
