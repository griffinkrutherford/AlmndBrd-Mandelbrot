#include <iostream>
#include <iomanip>
#include "../Base Class/Properties/Complex Number Handling (z, c).h"
#include "../Base Class/Properties/Iteration Parameters.h"
#include "../Base Class/Methods/Iterate Equation.cpp"
#include "../Base Class/Methods/Check Divergence.cpp"

/**
 * Example program to compute Mandelbrot set points within a region
 */
int main() {
    // Iteration parameters
    const int MAX_ITERATIONS = 1000;
    const double THRESHOLD = 2.0;
    
    // Region of complex plane to analyze
    const double X_MIN = -2.0;
    const double X_MAX = 1.0;
    const double Y_MIN = -1.5;
    const double Y_MAX = 1.5;
    
    // Resolution for the grid
    const int WIDTH = 80;
    const int HEIGHT = 40;
    
    // Create iteration parameters
    IterationParameters params(MAX_ITERATIONS, THRESHOLD);
    
    // Characters to display based on iteration count
    const char* CHARSET = " .:-=+*#%@";
    const int CHARSET_LENGTH = 10;
    
    // Step sizes for iterations
    double xStep = (X_MAX - X_MIN) / WIDTH;
    double yStep = (Y_MAX - Y_MIN) / HEIGHT;
    
    // Output header
    std::cout << "Simple Mandelbrot Set Visualization" << std::endl;
    std::cout << "Resolution: " << WIDTH << "x" << HEIGHT << std::endl;
    std::cout << "Max Iterations: " << MAX_ITERATIONS << std::endl;
    std::cout << "Coordinates: [" << X_MIN << "," << X_MAX << "] x ["
              << Y_MIN << "," << Y_MAX << "]" << std::endl << std::endl;
    
    // Generate and display Mandelbrot set
    for (int y = 0; y < HEIGHT; y++) {
        double imag = Y_MAX - y * yStep;
        
        for (int x = 0; x < WIDTH; x++) {
            double real = X_MIN + x * xStep;
            
            // Create point in complex plane
            Complex c(real, imag);
            
            // Perform iteration until divergence or max iterations
            int iterations = computeFull(c, params.maxIterations, params.threshold);
            
            // Map iteration count to character
            char displayChar;
            if (iterations == MAX_ITERATIONS) {
                // Point is likely in the Mandelbrot set
                displayChar = CHARSET[CHARSET_LENGTH - 1];
            } else {
                // Point diverged, map iterations to character
                int charIndex = static_cast<int>(iterations * (CHARSET_LENGTH - 1) / MAX_ITERATIONS);
                displayChar = CHARSET[charIndex];
            }
            
            // Output the character
            std::cout << displayChar;
        }
        std::cout << std::endl;
    }
    
    return 0;
}