#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Demonstration of the Mandelbrot set iteration formula z_{n+1} = z_n² + c
 * 
 * This file demonstrates the recursive formula that defines the Mandelbrot set,
 * showing step-by-step how iterations evolve for various points.
 */

/**
 * @brief Print complex number in a readable format
 * 
 * @param z Complex number to print
 * @param precision Number of decimal places to show
 * @return std::string Formatted complex number string
 */
std::string formatComplex(const Complex& z, int precision = 6) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision);
    
    double real = z.getReal();
    double imag = z.getImag();
    
    // Format the real part
    if (std::abs(real) < 1e-10) {
        oss << "0";
    } else {
        oss << real;
    }
    
    // Format the imaginary part
    if (std::abs(imag) < 1e-10) {
        // If both parts are zero, just return "0"
        if (std::abs(real) < 1e-10) {
            return "0";
        }
        return oss.str();
    }
    
    if (imag > 0) {
        oss << " + " << imag << "i";
    } else {
        oss << " - " << std::abs(imag) << "i";
    }
    
    return oss.str();
}

/**
 * @brief Demonstrate iteration formula for the Mandelbrot set
 * 
 * @param c Complex constant
 * @param iterations Number of iterations to perform
 */
void demonstrateMandelbrotFormula(const Complex& c, int iterations = 10) {
    std::cout << "Demonstrating Mandelbrot iteration formula: z_{n+1} = z_n² + c" << std::endl;
    std::cout << "For c = " << formatComplex(c) << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::setw(5) << "n" << " | " 
              << std::setw(30) << "z_n" << " | "
              << std::setw(30) << "z_n²" << " | "
              << std::setw(30) << "z_n² + c" << " | "
              << std::setw(12) << "|z_n|" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Start with z_0 = 0
    Complex z(0.0, 0.0);
    
    for (int n = 0; n <= iterations; n++) {
        // Calculate z_n²
        Complex z_squared = z * z;
        
        // Calculate z_{n+1} = z_n² + c
        Complex z_next = z_squared + c;
        
        // Calculate magnitude of z_n
        double magnitude = z.magnitude();
        
        // Print current iteration
        std::cout << std::setw(5) << n << " | " 
                  << std::setw(30) << formatComplex(z) << " | "
                  << std::setw(30) << formatComplex(z_squared) << " | "
                  << std::setw(30) << formatComplex(z_next) << " | "
                  << std::setw(12) << magnitude << std::endl;
        
        // Exit if diverging
        if (magnitude > 2.0) {
            std::cout << "Point escapes after " << n << " iterations! (|z| > 2)" << std::endl;
            break;
        }
        
        // Set z for next iteration
        z = z_next;
    }
    
    // Check if point might be in the set
    if (z.magnitude() <= 2.0) {
        std::cout << "Point might be in the Mandelbrot set after " << iterations << " iterations." << std::endl;
    }
}

/**
 * @brief Demonstrate the mathematical relationship between Mandelbrot and Julia sets
 * 
 * @param c Complex constant for the Julia set
 * @param startingPoints Vector of initial z values to trace
 * @param iterations Number of iterations to perform
 */
void demonstrateJuliaRelationship(
    const Complex& c, 
    const std::vector<Complex>& startingPoints,
    int iterations = 10
) {
    std::cout << "Demonstrating Julia set formula z_{n+1} = z_n² + c" << std::endl;
    std::cout << "For fixed c = " << formatComplex(c) << std::endl;
    std::cout << "Multiple starting points are traced to show the set's shape." << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Process each starting point
    for (size_t pointIndex = 0; pointIndex < startingPoints.size(); pointIndex++) {
        Complex z = startingPoints[pointIndex];
        
        std::cout << "Starting point " << pointIndex + 1 << ": z_0 = " << formatComplex(z) << std::endl;
        std::cout << std::setw(5) << "n" << " | " 
                  << std::setw(30) << "z_n" << " | "
                  << std::setw(12) << "|z_n|" << std::endl;
        std::cout << std::string(50, '-') << std::endl;
        
        // Track iterations for this starting point
        for (int n = 0; n <= iterations; n++) {
            // Print current point
            std::cout << std::setw(5) << n << " | " 
                      << std::setw(30) << formatComplex(z) << " | "
                      << std::setw(12) << z.magnitude() << std::endl;
            
            // Check if diverging
            if (z.magnitude() > 2.0) {
                std::cout << "Point escapes after " << n << " iterations!" << std::endl;
                break;
            }
            
            // Calculate next iteration: z_{n+1} = z_n² + c
            z = z * z + c;
        }
        
        // Check if point might be in the Julia set
        if (z.magnitude() <= 2.0) {
            std::cout << "Point might be in the Julia set after " << iterations << " iterations." << std::endl;
        }
        
        std::cout << std::endl;
    }
}

/**
 * @brief Demonstrate how the Mandelbrot set relates to parameter space
 * 
 * For each point c in the complex plane, the Mandelbrot set contains points
 * where the orbit of 0 under z -> z² + c doesn't escape to infinity.
 * 
 * @param testPoints Vector of c values to test
 * @param iterations Number of iterations to perform
 */
void demonstrateParameterSpace(
    const std::vector<Complex>& testPoints,
    int iterations = 20
) {
    std::cout << "Demonstrating Mandelbrot set as a map of parameter space:" << std::endl;
    std::cout << "For each point c, we check if the sequence z_0 = 0, z_{n+1} = z_n² + c escapes." << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Process each test point as a potential c value
    for (size_t pointIndex = 0; pointIndex < testPoints.size(); pointIndex++) {
        Complex c = testPoints[pointIndex];
        
        std::cout << "Test point " << pointIndex + 1 << ": c = " << formatComplex(c) << std::endl;
        
        // Start with z_0 = 0
        Complex z(0.0, 0.0);
        bool escaped = false;
        
        // Track iterations for this c value
        for (int n = 0; n < iterations; n++) {
            // Calculate next iteration: z_{n+1} = z_n² + c
            z = z * z + c;
            
            // Check if diverging
            if (z.magnitude() > 2.0) {
                std::cout << "  Escapes after " << n + 1 << " iterations (|z| = " 
                          << z.magnitude() << " > 2)" << std::endl;
                escaped = true;
                break;
            }
        }
        
        // Check if point might be in the Mandelbrot set
        if (!escaped) {
            std::cout << "  May be in the Mandelbrot set (no escape after " 
                      << iterations << " iterations, |z| = " << z.magnitude() << ")" << std::endl;
        }
        
        std::cout << std::endl;
    }
}

/**
 * @brief Main demonstration function
 */
int main() {
    // Demonstrate Mandelbrot iteration formula for a point outside the set
    Complex c1(-0.5, 0.6);
    demonstrateMandelbrotFormula(c1, 10);
    std::cout << std::endl;
    
    // Demonstrate Mandelbrot iteration formula for a point inside the set
    Complex c2(-1.0, 0.0);
    demonstrateMandelbrotFormula(c2, 20);
    std::cout << std::endl;
    
    // Demonstrate Julia set with multiple starting points
    Complex juliaC(-0.8, 0.156);
    std::vector<Complex> startingPoints = {
        Complex(0.0, 0.0),
        Complex(0.5, 0.0),
        Complex(0.0, 0.5),
        Complex(-0.5, 0.0),
        Complex(0.0, -0.5)
    };
    demonstrateJuliaRelationship(juliaC, startingPoints, 15);
    std::cout << std::endl;
    
    // Demonstrate parameter space exploration
    std::vector<Complex> testPoints = {
        Complex(0.0, 0.0),      // Center of main cardioid (in set)
        Complex(-1.0, 0.0),     // Center of period-2 bulb (in set)
        Complex(-2.0, 0.0),     // Just outside the set
        Complex(0.3, 0.6),      // Outside the set
        Complex(-0.75, 0.1),    // Near the boundary
        Complex(-1.77, 0.0)     // Mini Mandelbrot
    };
    demonstrateParameterSpace(testPoints, 30);
    
    return 0;
}