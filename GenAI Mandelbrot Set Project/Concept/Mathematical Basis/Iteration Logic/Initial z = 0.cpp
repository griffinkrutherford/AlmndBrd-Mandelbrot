#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Explanation and demonstration of why z₀ = 0 is used for the Mandelbrot set
 * 
 * This file explains the mathematical reasoning behind starting with z₀ = 0
 * for the Mandelbrot set calculation, and demonstrates the impact of alternative
 * starting values.
 */

/**
 * @brief Demonstrate why z₀ = 0 is used for the Mandelbrot set
 */
void explainInitialZero() {
    std::cout << "==== Why z₀ = 0 is Used for the Mandelbrot Set ====" << std::endl;
    std::cout << std::endl;
    std::cout << "Definition:" << std::endl;
    std::cout << "  The Mandelbrot set is defined as the set of complex numbers c for which" << std::endl;
    std::cout << "  the function f_c(z) = z² + c does not diverge when iterated from z = 0." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Mathematical Reasoning:" << std::endl;
    std::cout << "  1. The Mandelbrot set is a catalog of Julia sets." << std::endl;
    std::cout << "  2. For each c, there is a corresponding Julia set for f_c(z) = z² + c." << std::endl;
    std::cout << "  3. A key property: c is in the Mandelbrot set if and only if" << std::endl;
    std::cout << "     the corresponding Julia set is connected." << std::endl;
    std::cout << "  4. This property is determined by the behavior of the critical point." << std::endl;
    std::cout << "  5. For quadratic functions, the critical point is z = 0." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Critical Points:" << std::endl;
    std::cout << "  - Critical points are where the derivative f'(z) = 0." << std::endl;
    std::cout << "  - For f_c(z) = z² + c, the derivative is f'_c(z) = 2z." << std::endl;
    std::cout << "  - Setting 2z = 0 gives us z = 0 as the only critical point." << std::endl;
    std::cout << std::endl;
    
    std::cout << "The Importance of Critical Points:" << std::endl;
    std::cout << "  - The behavior of critical orbits determines the topology of Julia sets." << std::endl;
    std::cout << "  - If the orbit of z = 0 under f_c(z) = z² + c remains bounded," << std::endl;
    std::cout << "    then the Julia set for parameter c is connected." << std::endl;
    std::cout << "  - If the orbit of z = 0 escapes to infinity, the Julia set is disconnected" << std::endl;
    std::cout << "    (a Cantor dust)." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Summary:" << std::endl;
    std::cout << "  Starting with z₀ = 0 allows us to determine whether the" << std::endl;
    std::cout << "  corresponding Julia set is connected or disconnected, which is" << std::endl;
    std::cout << "  the fundamental property that defines the Mandelbrot set." << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Demonstrate the impact of using different initial z values
 * 
 * @param cValue Complex parameter to test
 * @param initialZValues Vector of different initial z values to try
 * @param iterations Number of iterations to perform
 */
void demonstrateDifferentInitialZ(
    const Complex& cValue, 
    const std::vector<Complex>& initialZValues,
    int iterations = 20
) {
    std::cout << "==== Impact of Different Initial z Values ====" << std::endl;
    std::cout << "Using c = " << cValue.getReal() << " + " << cValue.getImag() << "i" << std::endl;
    std::cout << std::endl;
    
    // Process each initial z value
    for (size_t i = 0; i < initialZValues.size(); i++) {
        Complex z = initialZValues[i];
        
        std::cout << "Starting with z₀ = " << z.getReal() << " + " << z.getImag() << "i:" << std::endl;
        
        bool escaped = false;
        int escapeIteration = 0;
        
        // Iterate the function f_c(z) = z² + c
        for (int n = 0; n < iterations; n++) {
            // Calculate next iteration
            z = z * z + cValue;
            
            // Check if the sequence has escaped
            if (z.magnitude() > 2.0) {
                escaped = true;
                escapeIteration = n + 1;
                break;
            }
        }
        
        // Report the outcome
        if (escaped) {
            std::cout << "  Sequence escapes after " << escapeIteration << " iterations." << std::endl;
        } else {
            std::cout << "  Sequence remains bounded after " << iterations << " iterations." << std::endl;
            std::cout << "  Final |z| = " << z.magnitude() << std::endl;
        }
        
        std::cout << std::endl;
    }
}

/**
 * @brief Demonstrate the relationship between different initial z values
 * 
 * @param cValues Vector of different c parameters to test
 * @param initialZValues Vector of different initial z values to try
 * @param iterations Number of iterations to perform
 */
void compareInitialZValues(
    const std::vector<Complex>& cValues,
    const std::vector<Complex>& initialZValues,
    int iterations = 30
) {
    std::cout << "==== Comparing Results with Different Initial z Values ====" << std::endl;
    std::cout << std::endl;
    
    // Process each c value
    for (size_t ci = 0; ci < cValues.size(); ci++) {
        Complex c = cValues[ci];
        
        std::cout << "For c = " << c.getReal() << " + " << c.getImag() << "i:" << std::endl;
        std::cout << std::string(50, '-') << std::endl;
        std::cout << std::setw(20) << "Initial z" << " | " 
                  << std::setw(15) << "Final |z|" << " | "
                  << std::setw(12) << "Iterations" << " | "
                  << std::setw(10) << "Outcome" << std::endl;
        std::cout << std::string(65, '-') << std::endl;
        
        // Test with each initial z value
        for (size_t zi = 0; zi < initialZValues.size(); zi++) {
            Complex z = initialZValues[zi];
            std::string initialZStr = std::to_string(z.getReal()) + " + " + 
                                     std::to_string(z.getImag()) + "i";
            
            bool escaped = false;
            int escapeIteration = 0;
            
            // Iterate the function f_c(z) = z² + c
            for (int n = 0; n < iterations; n++) {
                z = z * z + c;
                
                if (z.magnitude() > 2.0) {
                    escaped = true;
                    escapeIteration = n + 1;
                    break;
                }
            }
            
            // Report outcome in table format
            std::cout << std::setw(20) << initialZStr << " | " 
                      << std::setw(15) << z.magnitude() << " | ";
                      
            if (escaped) {
                std::cout << std::setw(12) << escapeIteration << " | "
                          << std::setw(10) << "Escapes" << std::endl;
            } else {
                std::cout << std::setw(12) << iterations << " | "
                          << std::setw(10) << "Bounded" << std::endl;
            }
        }
        
        std::cout << std::endl;
    }
}

/**
 * @brief Demonstrate the orbit of critical point z = 0
 * 
 * @param c Complex parameter
 * @param iterations Number of iterations to trace
 */
void traceCriticalOrbit(const Complex& c, int iterations = 20) {
    std::cout << "==== Tracing the Orbit of Critical Point z = 0 ====" << std::endl;
    std::cout << "For c = " << c.getReal() << " + " << c.getImag() << "i" << std::endl;
    std::cout << std::endl;
    
    // Start with critical point z = 0
    Complex z(0.0, 0.0);
    
    std::cout << std::setw(5) << "n" << " | " 
              << std::setw(30) << "z_n" << " | "
              << std::setw(15) << "|z_n|" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    // Print initial value
    std::cout << std::setw(5) << 0 << " | " 
              << std::setw(30) << "0 + 0i" << " | "
              << std::setw(15) << "0" << std::endl;
    
    // Trace the orbit of z = 0
    for (int n = 1; n <= iterations; n++) {
        // Calculate next iteration
        z = z * z + c;
        
        // Format the complex number
        std::string zStr = std::to_string(z.getReal()) + " + " + 
                           std::to_string(z.getImag()) + "i";
        
        // Print current iteration
        std::cout << std::setw(5) << n << " | " 
                  << std::setw(30) << zStr << " | "
                  << std::setw(15) << z.magnitude() << std::endl;
        
        // Exit if diverging
        if (z.magnitude() > 2.0) {
            std::cout << std::endl;
            std::cout << "The critical orbit escapes, so c is not in the Mandelbrot set." << std::endl;
            break;
        }
    }
    
    // Check if the orbit remained bounded
    if (z.magnitude() <= 2.0) {
        std::cout << std::endl;
        std::cout << "The critical orbit remains bounded, so c might be in the Mandelbrot set." << std::endl;
    }
}

/**
 * @brief Main demonstration function
 */
int main() {
    // Explain the mathematical basis
    explainInitialZero();
    std::cout << std::string(70, '=') << std::endl << std::endl;
    
    // Demonstrate the impact of different initial z values
    Complex testC(-0.75, 0.1);
    std::vector<Complex> initialZValues = {
        Complex(0.0, 0.0),     // Critical point
        Complex(1.0, 0.0),     // Simple non-zero value
        Complex(0.0, 1.0),     // Imaginary unit
        Complex(-1.0, 0.0),    // Negative real
        Complex(0.5, 0.5)      // Mixed real and imaginary
    };
    demonstrateDifferentInitialZ(testC, initialZValues, 20);
    std::cout << std::string(70, '=') << std::endl << std::endl;
    
    // Compare results across different c values and initial z values
    std::vector<Complex> cValues = {
        Complex(0.0, 0.0),     // Center of main cardioid (in the set)
        Complex(-1.0, 0.0),    // Period-2 bulb (in the set)
        Complex(-0.75, 0.1),   // Near boundary
        Complex(0.3, 0.5)      // Outside the set
    };
    compareInitialZValues(cValues, initialZValues, 30);
    std::cout << std::string(70, '=') << std::endl << std::endl;
    
    // Trace the orbit of critical point for points in and out of the set
    traceCriticalOrbit(Complex(-0.75, 0.1), 20);  // Near boundary
    std::cout << std::endl;
    traceCriticalOrbit(Complex(0.3, 0.5), 20);    // Outside the set
    
    return 0;
}