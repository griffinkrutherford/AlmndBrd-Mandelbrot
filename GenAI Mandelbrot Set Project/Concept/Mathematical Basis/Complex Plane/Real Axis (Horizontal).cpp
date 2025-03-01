#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Exploration of the real axis (horizontal) in the Mandelbrot set
 * 
 * This file demonstrates the properties of the Mandelbrot set along the real axis
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
 * @brief Explore the real axis of the Mandelbrot set
 * 
 * @param start Starting real value
 * @param end Ending real value
 * @param steps Number of steps to sample
 * @param maxIterations Maximum iterations for each point
 */
void exploreRealAxis(double start, double end, int steps, int maxIterations = 100) {
    std::cout << "=== Exploring the Real Axis of the Mandelbrot Set ===" << std::endl;
    std::cout << "Range: [" << start << ", " << end << "]" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(15) << "c (real)" << " | "
              << std::setw(12) << "Iterations" << " | "
              << std::setw(10) << "Status" << " | "
              << "Visualization" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    double step = (end - start) / steps;
    
    for (int i = 0; i <= steps; i++) {
        double real = start + i * step;
        Complex c(real, 0.0);
        
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
        
        std::cout << std::setw(15) << real << " | "
                  << std::setw(12) << iterations << " | "
                  << std::setw(10) << status << " | "
                  << viz << std::endl;
    }
}

/**
 * @brief Explain the significance of key points on the real axis
 */
void explainRealAxisPoints() {
    std::cout << std::endl;
    std::cout << "=== Key Points on the Real Axis ===" << std::endl;
    std::cout << std::endl;
    
    std::cout << "The real axis contains several mathematically significant points:" << std::endl;
    std::cout << std::endl;
    
    std::cout << "1. c = -2.0" << std::endl;
    std::cout << "   This is the leftmost point of the Mandelbrot set on the real axis." << std::endl;
    std::cout << "   For values c < -2.0, the sequence always escapes to infinity." << std::endl;
    std::cout << std::endl;
    
    std::cout << "2. c = -1.75" << std::endl;
    std::cout << "   This is approximately the center of the period-3 bulb." << std::endl;
    std::cout << "   Points in this region produce orbits that repeat every 3 iterations." << std::endl;
    std::cout << std::endl;
    
    std::cout << "3. c = -1.5" << std::endl;
    std::cout << "   This point is within the filaments between the main cardioid and period-2 bulb." << std::endl;
    std::cout << std::endl;
    
    std::cout << "4. c = -1.25" << std::endl;
    std::cout << "   This is approximately the center of the period-2 bulb." << std::endl;
    std::cout << "   Points in this region produce orbits that repeat every 2 iterations." << std::endl;
    std::cout << std::endl;
    
    std::cout << "5. c = -0.75" << std::endl;
    std::cout << "   This point is within the main cardioid of the Mandelbrot set." << std::endl;
    std::cout << "   The orbit approaches a single fixed point." << std::endl;
    std::cout << std::endl;
    
    std::cout << "6. c = 0.25" << std::endl;
    std::cout << "   This is the rightmost point of the Mandelbrot set on the real axis." << std::endl;
    std::cout << "   For values c > 0.25, the sequence always escapes to infinity." << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Demonstrate orbits for key points on the real axis
 * 
 * @param iterations Number of iterations to trace
 */
void demonstrateKeyPointOrbits(int iterations = 20) {
    std::vector<double> keyPoints = {-2.0, -1.75, -1.25, -0.75, 0.0, 0.25};
    
    for (double real : keyPoints) {
        Complex c(real, 0.0);
        Complex z(0.0, 0.0);
        
        std::cout << "=== Orbit for c = " << real << " ===" << std::endl;
        std::cout << std::setw(5) << "n" << " | "
                  << std::setw(20) << "z_n" << " | "
                  << std::setw(15) << "|z_n|" << std::endl;
        std::cout << std::string(45, '-') << std::endl;
        
        // Print initial value
        std::cout << std::setw(5) << 0 << " | "
                  << std::setw(20) << "0" << " | "
                  << std::setw(15) << "0" << std::endl;
        
        bool escaped = false;
        
        // Trace the orbit
        for (int n = 1; n <= iterations; n++) {
            // Calculate next iteration
            z = z * z + c;
            
            // Format the complex number
            std::string zStr;
            if (std::abs(z.getImag()) < 1e-10) {
                zStr = std::to_string(z.getReal());  // Essentially real
            } else {
                zStr = std::to_string(z.getReal()) + " + " + 
                       std::to_string(z.getImag()) + "i";
            }
            
            // Print current iteration
            std::cout << std::setw(5) << n << " | "
                      << std::setw(20) << zStr << " | "
                      << std::setw(15) << z.magnitude() << std::endl;
            
            // Check for escape
            if (z.magnitude() > 2.0) {
                std::cout << "Orbit escapes after " << n << " iterations." << std::endl;
                escaped = true;
                break;
            }
            
            // Check for cycle detection (simplified)
            if (n > 5 && z.magnitudeSquared() < 1e-10) {
                std::cout << "Orbit converges to 0 (fixed point)." << std::endl;
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
 * @brief Main function to demonstrate real axis properties
 */
int main() {
    // Explore the entire interesting region of the real axis
    exploreRealAxis(-2.25, 0.5, 30, 100);
    std::cout << std::endl;
    
    // Explain key points
    explainRealAxisPoints();
    std::cout << std::endl;
    
    // Demonstrate orbits for specific points
    demonstrateKeyPointOrbits(20);
    
    return 0;
}