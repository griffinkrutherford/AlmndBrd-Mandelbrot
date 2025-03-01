#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Analysis of divergence threshold in Mandelbrot computation
 * 
 * This file explains and demonstrates the mathematical basis for using a
 * pre-defined absolute value (typically 2.0) as the divergence threshold
 * in Mandelbrot set calculations.
 */

/**
 * @brief Calculate iterations for a point with the given threshold
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
 * @brief Explain the mathematical basis for the divergence threshold
 */
void explainDivergenceThreshold() {
    std::cout << "=== The Mathematics of Divergence Thresholds ===" << std::endl << std::endl;
    
    std::cout << "In Mandelbrot set calculations, we iterate z = z² + c starting with z₀ = 0," << std::endl;
    std::cout << "and we need to decide when the sequence {zₙ} has 'escaped' and will" << std::endl;
    std::cout << "definitely diverge to infinity." << std::endl << std::endl;
    
    std::cout << "The standard divergence threshold is |z| > 2, and here's why:" << std::endl << std::endl;
    
    std::cout << "1. Mathematical Proof:" << std::endl;
    std::cout << "   If |z| > 2 and |c| < |z|, then:" << std::endl;
    std::cout << "   |z²+c| ≥ |z²| - |c| > |z²| - |z| = |z|²-|z| = |z|(|z|-1)" << std::endl;
    std::cout << "   Since |z| > 2, we know |z|-1 > 1, so |z|(|z|-1) > |z|" << std::endl;
    std::cout << "   This means |z_(n+1)| > |zₙ| for all future iterations." << std::endl;
    std::cout << "   Therefore, once |z| > 2, the sequence will always diverge to infinity." << std::endl << std::endl;
    
    std::cout << "2. For the standard Mandelbrot set where c is a parameter being tested:" << std::endl;
    std::cout << "   - If we're iterating for a specific c, and |c| ≤ 2, then the threshold |z| > 2 is sufficient." << std::endl;
    std::cout << "   - If |c| > 2, then the sequence escapes immediately (since z₀ = 0, z₁ = c)." << std::endl << std::endl;
    
    std::cout << "3. Computational Efficiency:" << std::endl;
    std::cout << "   - Using |z| > 2 minimizes unnecessary iterations while ensuring accuracy." << std::endl;
    std::cout << "   - A larger threshold would waste computational resources on points that will definitely escape." << std::endl;
    std::cout << "   - A smaller threshold might incorrectly classify some points as escaping when they're in the set." << std::endl << std::endl;
}

/**
 * @brief Compare the effect of different divergence thresholds on the same point
 * 
 * @param c Complex parameter to test
 * @param thresholds Vector of thresholds to compare
 * @param maxIterations Maximum iterations
 */
void compareThresholds(
    const Complex& c, 
    const std::vector<double>& thresholds, 
    int maxIterations = 1000
) {
    std::cout << "=== Effect of Different Divergence Thresholds ===" << std::endl;
    std::cout << "Point: c = " << c.getReal() << " + " << c.getImag() << "i" << std::endl;
    std::cout << "Max Iterations: " << maxIterations << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    std::cout << std::setw(15) << "Threshold" << " | "
              << std::setw(15) << "Iterations" << " | "
              << std::setw(20) << "Execution Time (μs)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (double threshold : thresholds) {
        // Measure execution time
        auto start = std::chrono::high_resolution_clock::now();
        
        int result = calculateIterations(c, maxIterations, threshold);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::string status;
        if (result == maxIterations) {
            status = "In Set";
        } else {
            status = std::to_string(result);
        }
        
        std::cout << std::setw(15) << threshold << " | "
                  << std::setw(15) << status << " | "
                  << std::setw(20) << duration.count() << std::endl;
    }
}

/**
 * @brief Demonstrate the behavior of z = z² + c for |z| > 2
 * 
 * @param c Complex parameter
 * @param initialZ Initial z value with |z| > 2
 * @param iterations Number of iterations to trace
 */
void demonstrateDivergenceBehavior(
    const Complex& c, 
    const Complex& initialZ, 
    int iterations = 10
) {
    std::cout << std::endl;
    std::cout << "=== Demonstrating Divergence for |z| > 2 ===" << std::endl;
    std::cout << "c = " << c.getReal() << " + " << c.getImag() << "i" << std::endl;
    std::cout << "Initial z = " << initialZ.getReal() << " + " << initialZ.getImag() << "i" << std::endl;
    std::cout << "Initial |z| = " << initialZ.magnitude() << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(5) << "n" << " | "
              << std::setw(30) << "zₙ" << " | "
              << std::setw(15) << "|zₙ|" << " | "
              << std::setw(15) << "Growth" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    Complex z = initialZ;
    double previousMagnitude = z.magnitude();
    
    // Print initial value
    std::cout << std::setw(5) << 0 << " | "
              << std::setw(30) << (std::to_string(z.getReal()) + " + " + std::to_string(z.getImag()) + "i") << " | "
              << std::setw(15) << previousMagnitude << " | "
              << std::setw(15) << "-" << std::endl;
    
    // Trace iterations
    for (int n = 1; n <= iterations; n++) {
        // Calculate next iteration
        z = z * z + c;
        
        // Calculate magnitude and growth
        double magnitude = z.magnitude();
        double growth = magnitude / previousMagnitude;
        
        // Format the complex number
        std::string zStr = std::to_string(z.getReal()) + " + " + std::to_string(z.getImag()) + "i";
        
        // Print current iteration
        std::cout << std::setw(5) << n << " | "
                  << std::setw(30) << zStr << " | "
                  << std::setw(15) << magnitude << " | "
                  << std::setw(15) << growth << std::endl;
        
        previousMagnitude = magnitude;
    }
    
    std::cout << std::endl;
    std::cout << "Notice that once |z| > 2, each iteration causes |z| to grow larger," << std::endl;
    std::cout << "demonstrating that the sequence is guaranteed to diverge to infinity." << std::endl;
}

/**
 * @brief Analyze the accuracy of different thresholds
 * 
 * @param testPoints Vector of points to test
 * @param thresholds Vector of thresholds to compare
 * @param referenceThreshold Reference threshold for "ground truth"
 * @param maxIterations Maximum iterations
 */
void analyzeThresholdAccuracy(
    const std::vector<Complex>& testPoints,
    const std::vector<double>& thresholds,
    double referenceThreshold,
    int maxIterations = 1000
) {
    std::cout << std::endl;
    std::cout << "=== Threshold vs. Accuracy Analysis ===" << std::endl;
    std::cout << "Using threshold = " << referenceThreshold << " and " << maxIterations << " iterations as reference" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    std::cout << std::setw(15) << "Threshold" << " | "
              << std::setw(20) << "Classification Errors" << " | "
              << std::setw(15) << "Error Rate (%)" << " | "
              << std::setw(15) << "Performance (%)" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // Calculate reference results
    std::vector<bool> reference;
    auto refStart = std::chrono::high_resolution_clock::now();
    
    for (const auto& point : testPoints) {
        int result = calculateIterations(point, maxIterations, referenceThreshold);
        reference.push_back(result == maxIterations);
    }
    
    auto refEnd = std::chrono::high_resolution_clock::now();
    auto refDuration = std::chrono::duration_cast<std::chrono::microseconds>(refEnd - refStart).count();
    
    // Test different thresholds
    for (double threshold : thresholds) {
        // Skip the reference threshold
        if (std::abs(threshold - referenceThreshold) < 1e-10) continue;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        int errorCount = 0;
        
        for (size_t i = 0; i < testPoints.size(); i++) {
            int result = calculateIterations(testPoints[i], maxIterations, threshold);
            bool inSet = (result == maxIterations);
            
            // Check if classification differs from reference
            if (inSet != reference[i]) {
                errorCount++;
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        double errorRate = 100.0 * errorCount / testPoints.size();
        double performanceRatio = 100.0 * refDuration / duration;
        
        std::cout << std::setw(15) << threshold << " | "
                  << std::setw(20) << errorCount << "/" << testPoints.size() << " | "
                  << std::setw(15) << std::fixed << std::setprecision(2) << errorRate << " | "
                  << std::setw(15) << std::fixed << std::setprecision(2) << performanceRatio << std::endl;
    }
}

/**
 * @brief Main function to demonstrate divergence threshold concepts
 */
int main() {
    // Explain the mathematics behind divergence threshold
    explainDivergenceThreshold();
    std::cout << std::endl;
    
    // Compare different thresholds for a boundary point
    Complex boundaryPoint(-0.75, 0.1);
    std::vector<double> thresholds = {1.0, 1.5, 2.0, 4.0, 10.0, 100.0};
    compareThresholds(boundaryPoint, thresholds, 1000);
    
    // Demonstrate the behavior for |z| > 2
    Complex c(-0.75, 0.1);  // A point in the interesting region
    Complex initialZ(2.1, 0.0);  // |z| > 2
    demonstrateDivergenceBehavior(c, initialZ, 10);
    
    // Generate test points for accuracy analysis
    std::vector<Complex> testPoints;
    for (int i = 0; i < 10; i++) {
        double real = -2.0 + i * 0.25;
        for (int j = 0; j < 8; j++) {
            double imag = -1.0 + j * 0.25;
            testPoints.push_back(Complex(real, imag));
        }
    }
    
    // Analyze accuracy of different thresholds
    std::vector<double> accuracyThresholds = {1.0, 1.5, 2.0, 3.0, 4.0, 10.0};
    analyzeThresholdAccuracy(testPoints, accuracyThresholds, 2.0, 1000);
    
    return 0;
}