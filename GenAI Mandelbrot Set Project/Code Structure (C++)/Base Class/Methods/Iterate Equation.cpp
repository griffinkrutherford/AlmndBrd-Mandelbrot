#include "../Properties/Complex Number Handling (z, c).h"

/**
 * Performs a single iteration of the Mandelbrot equation: z = z² + c
 * 
 * @param z Current z value, will be updated with the new value
 * @param c The constant complex value for this point
 */
void iterateEquation(Complex& z, const Complex& c) {
    // z = z² + c
    z = z * z + c;
}

/**
 * Performs multiple iterations of the Mandelbrot equation until
 * either divergence is detected or maximum iterations are reached
 * 
 * @param c The constant complex value for this point
 * @param maxIterations Maximum number of iterations to perform
 * @param threshold Divergence threshold
 * @return Number of iterations performed before divergence or maxIterations
 */
int computeFull(const Complex& c, int maxIterations, double threshold) {
    // Start with z = 0
    Complex z(0.0, 0.0);
    
    // Iterate until divergence or max iterations
    for (int i = 0; i < maxIterations; i++) {
        // Perform one iteration
        iterateEquation(z, c);
        
        // Check for divergence
        if (z.magnitudeSquared() > threshold * threshold) {
            return i + 1;
        }
    }
    
    // If we reach here, the point did not diverge within maxIterations
    return maxIterations;
}