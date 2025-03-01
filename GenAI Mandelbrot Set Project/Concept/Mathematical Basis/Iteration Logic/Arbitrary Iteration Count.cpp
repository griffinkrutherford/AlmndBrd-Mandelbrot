#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Analysis of arbitrary iteration counts in Mandelbrot computation
 * 
 * This file demonstrates the impact of different iteration limits on
 * the accuracy and performance of Mandelbrot set computations.
 */

/**
 * @brief Calculate iterations for a point with the given parameters
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum iterations
 * @param threshold Divergence threshold
 * @return int Iteration count
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
 * @brief Compare the effect of different iteration limits on the same point
 * 
 * @param c Complex parameter to test
 * @param iterations Vector of iteration limits to compare
 */
void compareIterationLimits(const Complex& c, const std::vector<int>& iterations) {
    std::cout << "=== Effect of Different Iteration Limits ===" << std::endl;
    std::cout << "Point: c = " << c.getReal() << " + " << c.getImag() << "i" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    std::cout << std::setw(15) << "Max Iterations" << " | "
              << std::setw(15) << "Result" << " | "
              << std::setw(20) << "Execution Time (μs)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (int maxIter : iterations) {
        // Measure execution time
        auto start = std::chrono::high_resolution_clock::now();
        
        int result = calculateIterations(c, maxIter, 2.0);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::string status;
        if (result == maxIter) {
            status = "In Set";
        } else {
            status = "Escapes (" + std::to_string(result) + ")";
        }
        
        std::cout << std::setw(15) << maxIter << " | "
                  << std::setw(15) << status << " | "
                  << std::setw(20) << duration.count() << std::endl;
    }
}

/**
 * @brief Demonstrate how iteration count affects the detail level of boundary regions
 * 
 * @param center Center point of the region
 * @param size Size of the region
 * @param resolution Resolution of the grid
 * @param iterations Vector of iteration limits to compare
 */
void demonstrateBoundaryDetail(
    const Complex& center, 
    double size, 
    int resolution, 
    const std::vector<int>& iterations
) {
    std::cout << std::endl;
    std::cout << "=== Boundary Detail with Different Iteration Limits ===" << std::endl;
    std::cout << "Region centered at c = " << center.getReal() << " + " << center.getImag() << "i" << std::endl;
    std::cout << "Size: " << size << " x " << size << ", Resolution: " << resolution << " x " << resolution << std::endl;
    std::cout << std::endl;
    
    // Constants for display
    const char* CHARSET = " .:+*#%@";
    const int CHARSET_LENGTH = 8;
    
    // Calculate the step size
    double step = size / resolution;
    
    // For each iteration limit
    for (int maxIter : iterations) {
        std::cout << "-- Max Iterations: " << maxIter << " --" << std::endl;
        
        // Measure execution time
        auto start = std::chrono::high_resolution_clock::now();
        
        // Generate a grid of results
        std::vector<std::vector<int>> grid(resolution, std::vector<int>(resolution));
        
        for (int y = 0; y < resolution; y++) {
            double imag = center.getImag() + (y - resolution / 2) * step;
            
            for (int x = 0; x < resolution; x++) {
                double real = center.getReal() + (x - resolution / 2) * step;
                Complex c(real, imag);
                
                grid[y][x] = calculateIterations(c, maxIter, 2.0);
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // Display the grid
        for (int y = 0; y < resolution; y++) {
            for (int x = 0; x < resolution; x++) {
                int iterations = grid[y][x];
                
                // Map iterations to character
                char displayChar;
                if (iterations == maxIter) {
                    // Point is likely in the Mandelbrot set
                    displayChar = CHARSET[CHARSET_LENGTH - 1];
                } else {
                    // Point escaped, map iterations to character
                    int charIndex = static_cast<int>(static_cast<double>(iterations) * (CHARSET_LENGTH - 1) / maxIter);
                    charIndex = std::min(std::max(charIndex, 0), CHARSET_LENGTH - 1);
                    displayChar = CHARSET[charIndex];
                }
                
                std::cout << displayChar;
            }
            std::cout << std::endl;
        }
        
        std::cout << "Computation time: " << duration.count() << " ms" << std::endl;
        std::cout << std::endl;
    }
}

/**
 * @brief Analyze the relationship between iteration limit and accuracy
 * 
 * @param testPoints Vector of points to test
 * @param maxIterations Maximum iteration limit for "ground truth"
 */
void analyzeAccuracy(const std::vector<Complex>& testPoints, int maxIterations) {
    std::cout << "=== Iteration Limit vs. Accuracy Analysis ===" << std::endl;
    std::cout << "Using " << maxIterations << " iterations as 'ground truth'" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(15) << "Iterations" << " | "
              << std::setw(15) << "Avg. Error (%)" << " | "
              << std::setw(15) << "Classification" << " | "
              << std::setw(15) << "Performance" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // Calculate "ground truth" results
    std::vector<int> groundTruth;
    std::vector<bool> inSet;
    
    for (const auto& point : testPoints) {
        int result = calculateIterations(point, maxIterations, 2.0);
        groundTruth.push_back(result);
        inSet.push_back(result == maxIterations);
    }
    
    // Test different iteration limits
    std::vector<int> iterationLimits = {10, 20, 50, 100, 200, 500, 1000, 2000};
    
    for (int iter : iterationLimits) {
        // Skip the ground truth comparison with itself
        if (iter == maxIterations) continue;
        
        // Measure performance
        auto start = std::chrono::high_resolution_clock::now();
        
        std::vector<int> results;
        std::vector<bool> classifications;
        int errorCount = 0;
        double totalError = 0.0;
        
        for (const auto& point : testPoints) {
            int result = calculateIterations(point, iter, 2.0);
            results.push_back(result);
            
            bool isInSet = (result == iter);
            classifications.push_back(isInSet);
            
            // Calculate error for non-escaping points
            if (result < iter) {
                double error = 1.0 - static_cast<double>(result) / groundTruth[results.size() - 1];
                totalError += std::abs(error) * 100.0;  // Error as percentage
            }
            
            // Check if classification is different from ground truth
            if (isInSet != inSet[results.size() - 1]) {
                errorCount++;
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        // Calculate average error percentage for escaping points
        double avgError = totalError / testPoints.size();
        
        // Calculate classification error percentage
        double classError = 100.0 * errorCount / testPoints.size();
        
        // Calculate performance ratio relative to ground truth
        double performanceRatio = 100.0 * iter / maxIterations;
        
        std::cout << std::setw(15) << iter << " | "
                  << std::setw(15) << std::fixed << std::setprecision(2) << avgError << " | "
                  << std::setw(15) << classError << "% | "
                  << std::setw(15) << performanceRatio << "%" << std::endl;
    }
}

/**
 * @brief Recommend an appropriate iteration limit based on error tolerance
 * 
 * @param testPoints Vector of points to test
 * @param errorTolerance Acceptable error percentage
 * @param maxLimit Upper limit for iterations
 * @return int Recommended iteration limit
 */
int recommendIterationLimit(
    const std::vector<Complex>& testPoints, 
    double errorTolerance, 
    int maxLimit = 5000
) {
    std::cout << std::endl;
    std::cout << "=== Iteration Limit Recommendation ===" << std::endl;
    std::cout << "Error tolerance: " << errorTolerance << "%" << std::endl;
    std::cout << std::endl;
    
    // Test different iteration limits
    std::vector<int> iterationLimits = {
        10, 20, 50, 100, 200, 500, 1000, 2000, 5000
    };
    
    // Iterate through each limit and check if it meets the error tolerance
    for (int iter : iterationLimits) {
        if (iter > maxLimit) break;
        
        // Calculate ground truth with a much higher iteration limit
        int groundTruthLimit = std::min(iter * 10, maxLimit);
        
        int errorCount = 0;
        
        for (const auto& point : testPoints) {
            int result = calculateIterations(point, iter, 2.0);
            int groundTruth = calculateIterations(point, groundTruthLimit, 2.0);
            
            bool isInSet = (result == iter);
            bool groundTruthInSet = (groundTruth == groundTruthLimit);
            
            // Check if classification is different from ground truth
            if (isInSet != groundTruthInSet) {
                errorCount++;
            }
        }
        
        // Calculate classification error percentage
        double classError = 100.0 * errorCount / testPoints.size();
        
        std::cout << "Iterations: " << iter << ", Error: " << std::fixed << std::setprecision(2) << classError << "%" << std::endl;
        
        // Check if this meets the error tolerance
        if (classError <= errorTolerance) {
            std::cout << "Recommendation: " << iter << " iterations" << std::endl;
            std::cout << "This provides error below the " << errorTolerance << "% tolerance threshold." << std::endl;
            return iter;
        }
    }
    
    // If we get here, no iteration limit met the tolerance
    std::cout << "Recommendation: Use maximum available iterations (" << maxLimit << ")" << std::endl;
    std::cout << "Even this may not meet the " << errorTolerance << "% error tolerance for all points." << std::endl;
    return maxLimit;
}

/**
 * @brief Main function to demonstrate arbitrary iteration count properties
 */
int main() {
    // Compare different iteration limits for a boundary point
    Complex boundaryPoint(-0.75, 0.1);
    std::vector<int> iterationLimits = {10, 50, 100, 500, 1000, 5000, 10000};
    compareIterationLimits(boundaryPoint, iterationLimits);
    
    // Demonstrate boundary detail for different iteration limits
    Complex boundaryRegion(-0.75, 0.1);
    std::vector<int> detailIterations = {20, 100, 500};
    demonstrateBoundaryDetail(boundaryRegion, 0.1, 40, detailIterations);
    
    // Generate test points for accuracy analysis
    std::vector<Complex> testPoints;
    for (int i = 0; i < 20; i++) {
        double real = -2.0 + i * 0.15;
        for (int j = 0; j < 10; j++) {
            double imag = -1.0 + j * 0.2;
            testPoints.push_back(Complex(real, imag));
        }
    }
    
    // Analyze accuracy with different iteration limits
    analyzeAccuracy(testPoints, 5000);
    
    // Recommend an iteration limit based on error tolerance
    recommendIterationLimit(testPoints, 1.0, 10000);
    
    return 0;
}