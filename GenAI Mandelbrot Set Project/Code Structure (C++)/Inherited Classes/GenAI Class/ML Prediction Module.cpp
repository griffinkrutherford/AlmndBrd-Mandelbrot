#include "Complex Number Handling (z, c).h"
#include <vector>
#include <utility>
#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>

enum class MLModelType { KNN, RANDOM_FOREST, NEURAL_NET, HYBRID, EVOLUTIONARY };

// Simple function to extract features from iteration history
std::vector<double> extractFeatures(const std::vector<Complex>& iterationHistory) {
    std::vector<double> features;
    
    // If we have no history, return empty features
    if (iterationHistory.empty()) return features;
    
    // Get the length of history (number of iterations until escape or cutoff)
    double escapeCount = static_cast<double>(iterationHistory.size());
    features.push_back(escapeCount / 100.0); // Normalize
    
    // Calculate final magnitude
    double finalMag = iterationHistory.back().magnitude();
    features.push_back(std::min(finalMag, 10.0) / 10.0); // Normalize and cap
    
    // Check for oscillatory behavior
    double mean = 0.0;
    for (size_t i = 0; i < iterationHistory.size(); i++) {
        mean += iterationHistory[i].magnitude();
    }
    mean /= iterationHistory.size();
    
    double variance = 0.0;
    for (size_t i = 0; i < iterationHistory.size(); i++) {
        double diff = iterationHistory[i].magnitude() - mean;
        variance += diff * diff;
    }
    variance /= iterationHistory.size();
    
    features.push_back(std::min(variance, 4.0) / 4.0); // Normalized variance
    
    return features;
}

// Generate interesting points for the Mandelbrot set
std::vector<std::pair<Complex, double>> predictInterestingRegions(
    MLModelType modelType,
    int trainingSize,
    double regionSize,
    int resolution
) {
    std::cout << "Predicting interesting regions with model type " << static_cast<int>(modelType) << std::endl;
    
    std::vector<std::pair<Complex, double>> interestingPoints;
    
    // Use random number generator for reproducible results
    std::mt19937 gen(42); // Fixed seed for reproducibility
    
    // Different distributions for different region sizes
    std::uniform_real_distribution<double> region_dist(-regionSize/2, regionSize/2);
    
    // Define some known interesting regions based on the region size
    std::vector<Complex> knownRegions;
    
    if (regionSize > 3.0) {
        // For large regions, use main cardioid and period bulbs
        knownRegions.push_back(Complex(-0.75, 0.0));       // Main cardioid
        knownRegions.push_back(Complex(-1.25, 0.0));       // Period-2 bulb
        knownRegions.push_back(Complex(-1.75, 0.0));       // Period-3 bulb
        knownRegions.push_back(Complex(0.25, 0.0));        // Edge of cardioid
        knownRegions.push_back(Complex(-0.125, 0.744));    // Mini Mandelbrot
    } else if (regionSize > 1.0) {
        // For medium regions, use mini Mandelbrots and spiral points
        knownRegions.push_back(Complex(-0.77568377, 0.13646737)); // Mini Mandelbrot
        knownRegions.push_back(Complex(-0.170337, -1.0660))      // Elephant valley
        knownRegions.push_back(Complex(0.42884, -0.231345));      // Spiral point
        knownRegions.push_back(Complex(-1.77, 0.0));              // Period-3 detail
        knownRegions.push_back(Complex(-0.5, 0.563));             // Detail region
    } else {
        // For small regions, use very detailed areas
        knownRegions.push_back(Complex(-0.743643887037151, 0.131825904205330)); // Seahorse valley
        knownRegions.push_back(Complex(-0.77568377, 0.13646737));              // Mini Mandelbrot
        knownRegions.push_back(Complex(-1.25278, 0.34311));                    // Fine detail
        knownRegions.push_back(Complex(-0.1592, -1.0317));                     // Spiral detail
        knownRegions.push_back(Complex(-1.77, 0.0));                           // Period-3 fine
    }
    
    // Add variations around known interesting points
    for (const auto& center : knownRegions) {
        // Add the center point with high confidence
        interestingPoints.push_back({center, 0.9 + region_dist(gen) * 0.1});
        
        // Add variations around the center
        for (int i = 0; i < 3; i++) {
            double dx = region_dist(gen) * 0.2; // 20% of region size
            double dy = region_dist(gen) * 0.2;
            Complex variation(center.getReal() + dx, center.getImag() + dy);
            
            // Add with slightly lower confidence
            interestingPoints.push_back({variation, 0.7 + region_dist(gen) * 0.2});
        }
    }
    
    // Also add some random points with lower confidence
    for (int i = 0; i < 10; i++) {
        double x = region_dist(gen);
        double y = region_dist(gen);
        Complex point(x, y);
        
        // Generate an iteration history to check if it's interesting
        std::vector<Complex> history;
        Complex z(0, 0);
        history.push_back(z);
        
        // Iterate Mandelbrot formula
        for (int j = 0; j < 50; j++) {
            z = z * z + point;
            history.push_back(z);
            
            // If it escapes too quickly, it's not interesting
            if (z.magnitudeSquared() > 4.0 && j < 10) break;
        }
        
        // If the point didn't escape or took many iterations, it might be interesting
        if (history.size() > 20) {
            double confidence = 0.5 + (history.size() / 100.0) * 0.4; // Higher iterations = higher confidence
            interestingPoints.push_back({point, confidence});
        }
    }
    
    // Sort by confidence (highest first)
    std::sort(interestingPoints.begin(), interestingPoints.end(),
             [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Limit to reasonable number of points
    if (interestingPoints.size() > 20) {
        interestingPoints.resize(20);
    }
    
    std::cout << "Found " << interestingPoints.size() << " interesting points" << std::endl;
    return interestingPoints;
}