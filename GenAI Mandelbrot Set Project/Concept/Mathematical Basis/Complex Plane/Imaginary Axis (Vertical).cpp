#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Exploration of the imaginary axis (vertical) in the Mandelbrot set
 * 
 * This file demonstrates the properties of the Mandelbrot set along the imaginary axis
 * and explains the mathematical significance of key points.
 */

/**
 * @brief Calculate the iteration count for a point in the Mandelbrot set
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum number of iterations
 * @param threshold Divergence threshold
 * @return int Number of iterations before escape, or maxIterations if bounded
 */
int calculateIterations(const Complex& c, int maxIterations, double threshold) {
    Complex z(0.0, 0.0);
    double thresholdSquared = threshold * threshold;
    
    for (int i = 0; i < maxIterations; i++) {
        z = z * z + c;
        
        if (z.magnitudeSquared() > thresholdSquared) {
            return i + 1;
        }
    }
    
    return maxIterations;
}

/**
 * @brief Explore the imaginary axis of the Mandelbrot set
 * 
 * @param start Starting imaginary value
 * @param end Ending imaginary value
 * @param steps Number of steps to sample
 * @param maxIterations Maximum iterations for each point
 */
void exploreImaginaryAxis(double start, double end, int steps, int maxIterations = 100) {
    std::cout << "=== Exploring the Imaginary Axis of the Mandelbrot Set ===" << std::endl;
    std::cout << "Range: [" << start << "i, " << end << "i]" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(15) << "c (imaginary)" << " | "
              << std::setw(12) << "Iterations" << " | "
              << std::setw(10) << "Status" << " | "
              << "Visualization" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    double step = (end - start) / steps;
    
    for (int i = 0; i <= steps; i++) {
        double imag = start + i * step;
        Complex c(0.0, imag);
        
        int iterations = calculateIterations(c, maxIterations, 2.0);
        
        std::string status = (iterations == maxIterations) ? "In Set" : "Escapes";
        
        // Create a simple visualization of iteration count
        std::string viz;
        if (iterations == maxIterations) {
            viz = std::string(20, '#');  // Likely in the set
        } else {
            int barLength = static_cast<int>(20.0 * iterations / maxIterations);
            viz = std::string(barLength, '=');
        }
        
        std::cout << std::setw(15) << imag << "i | "
                  << std::setw(12) << iterations << " | "
                  << std::setw(10) << status << " | "
                  << viz << std::endl;
    }
}

/**
 * @brief Explain the significance of key points on the imaginary axis
 */
void explainImaginaryAxisPoints() {
    std::cout << std::endl;
    std::cout << "=== Key Points on the Imaginary Axis ===" << std::endl;
    std::cout << std::endl;
    
    std::cout << "The imaginary axis reveals the vertical symmetry of the Mandelbrot set:" << std::endl;
    std::cout << std::endl;
    
    std::cout << "1. c = 0i" << std::endl;
    std::cout << "   This is the origin of the complex plane." << std::endl;
    std::cout << "   The point is within the main cardioid of the Mandelbrot set." << std::endl;
    std::cout << std::endl;
    
    std::cout << "2. c = ±i" << std::endl;
    std::cout << "   These points lie at the boundary of the Mandelbrot set." << std::endl;
    std::cout << "   They mark the transition between bounded and unbounded behavior." << std::endl;
    std::cout << std::endl;
    
    std::cout << "3. c = ±1.1i" << std::endl;
    std::cout << "   These points are outside the Mandelbrot set." << std::endl;
    std::cout << "   The orbit escapes to infinity relatively quickly." << std::endl;
    std::cout << std::endl;
    
    std::cout << "4. Symmetry along the real axis" << std::endl;
    std::cout << "   The Mandelbrot set is symmetric about the real axis." << std::endl;
    std::cout << "   This means if c = a + bi is in the set, then c = a - bi is also in the set." << std::endl;
    std::cout << "   This is because complex conjugation preserves the behavior of z² + c." << std::endl;
    std::cout << std::endl;
    
    std::cout << "5. The boundary at i" << std::endl;
    std::cout << "   The boundary at c = i marks an interesting region of the Mandelbrot set." << std::endl;
    std::cout << "   The dynamics of orbits change significantly around this value." << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Demonstrate the symmetry around the real axis
 * 
 * @param imagValues Vector of imaginary values to test
 * @param maxIterations Maximum iterations for each point
 */
void demonstrateSymmetry(const std::vector<double>& imagValues, int maxIterations = 100) {
    std::cout << "=== Demonstrating Symmetry Around the Real Axis ===" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(15) << "c (positive)" << " | "
              << std::setw(12) << "Iterations" << " | "
              << std::setw(15) << "c (negative)" << " | "
              << std::setw(12) << "Iterations" << " | "
              << std::setw(10) << "Symmetric?" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (double imag : imagValues) {
        Complex cPositive(0.0, imag);
        Complex cNegative(0.0, -imag);
        
        int iterPositive = calculateIterations(cPositive, maxIterations, 2.0);
        int iterNegative = calculateIterations(cNegative, maxIterations, 2.0);
        
        std::string symmetric = (iterPositive == iterNegative) ? "Yes" : "No";
        
        std::cout << std::setw(15) << imag << "i | "
                  << std::setw(12) << iterPositive << " | "
                  << std::setw(15) << -imag << "i | "
                  << std::setw(12) << iterNegative << " | "
                  << std::setw(10) << symmetric << std::endl;
    }
}

/**
 * @brief Demonstrate orbits for key points on the imaginary axis
 * 
 * @param iterations Number of iterations to trace
 */
void demonstrateKeyPointOrbits(int iterations = 20) {
    std::vector<double> keyPoints = {0.0, 1.0, -1.0, 1.1, -1.1};
    
    for (double imag : keyPoints) {
        Complex c(0.0, imag);
        Complex z(0.0, 0.0);
        
        std::cout << "=== Orbit for c = " << imag << "i ===" << std::endl;
        std::cout << std::setw(5) << "n" << " | "
                  << std::setw(30) << "z_n" << " | "
                  << std::setw(15) << "|z_n|" << std::endl;
        std::cout << std::string(55, '-') << std::endl;
        
        // Print initial value
        std::cout << std::setw(5) << 0 << " | "
                  << std::setw(30) << "0" << " | "
                  << std::setw(15) << "0" << std::endl;
        
        bool escaped = false;
        
        // Trace the orbit
        for (int n = 1; n <= iterations; n++) {
            // Calculate next iteration
            z = z * z + c;
            
            // Format the complex number
            std::string zStr;
            if (std::abs(z.getReal()) < 1e-10) {
                zStr = std::to_string(z.getImag()) + "i";  // Essentially imaginary
            } else if (std::abs(z.getImag()) < 1e-10) {
                zStr = std::to_string(z.getReal());  // Essentially real
            } else {
                zStr = std::to_string(z.getReal()) + " + " + 
                       std::to_string(z.getImag()) + "i";
            }
            
            // Print current iteration
            std::cout << std::setw(5) << n << " | "
                      << std::setw(30) << zStr << " | "
                      << std::setw(15) << z.magnitude() << std::endl;
            
            // Check for escape
            if (z.magnitude() > 2.0) {
                std::cout << "Orbit escapes after " << n << " iterations." << std::endl;
                escaped = true;
                break;
            }
        }
        
        if (!escaped) {
            std::cout << "Orbit remains bounded after " << iterations << " iterations." << std::endl;
        }
        
        std::cout << std::endl;
    }
}

/**
 * @brief Compare the behavior of pure imaginary vs. pure real inputs
 * 
 * @param values Vector of magnitude values to test
 * @param maxIterations Maximum iterations for each point
 */
void compareRealAndImaginary(const std::vector<double>& values, int maxIterations = 100) {
    std::cout << "=== Comparing Real vs. Imaginary Inputs ===" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(10) << "Value" << " | "
              << std::setw(15) << "Real Input" << " | "
              << std::setw(15) << "Iterations" << " | "
              << std::setw(15) << "Imag Input" << " | "
              << std::setw(15) << "Iterations" << std::endl;
    std::cout << std::string(75, '-') << std::endl;
    
    for (double value : values) {
        Complex cReal(value, 0.0);
        Complex cImag(0.0, value);
        
        int iterReal = calculateIterations(cReal, maxIterations, 2.0);
        int iterImag = calculateIterations(cImag, maxIterations, 2.0);
        
        std::cout << std::setw(10) << value << " | "
                  << std::setw(15) << value << " | "
                  << std::setw(15) << iterReal << " | "
                  << std::setw(15) << value << "i | "
                  << std::setw(15) << iterImag << std::endl;
    }
}

/**
 * @brief Main function to demonstrate imaginary axis properties
 */
int main() {
    // Explore the interesting region of the imaginary axis
    exploreImaginaryAxis(-2.0, 2.0, 40, 100);
    std::cout << std::endl;
    
    // Explain key points
    explainImaginaryAxisPoints();
    std::cout << std::endl;
    
    // Demonstrate symmetry around the real axis
    std::vector<double> symmetryValues = {0.0, 0.5, 1.0, 1.5, 2.0};
    demonstrateSymmetry(symmetryValues, 100);
    std::cout << std::endl;
    
    // Demonstrate orbits for specific points
    demonstrateKeyPointOrbits(20);
    std::cout << std::endl;
    
    // Compare real vs. imaginary inputs
    std::vector<double> compareValues = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5};
    compareRealAndImaginary(compareValues, 100);
    
    return 0;
}