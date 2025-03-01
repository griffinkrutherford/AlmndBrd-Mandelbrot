#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <map>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"
#include "../../Code Structure (C++)/Inherited Classes/GenAI Class/ML Prediction Module.cpp"

/**
 * @brief Demonstration of machine learning predictive capabilities for fractal exploration
 * 
 * This file shows how machine learning can be used to predict interesting
 * regions and characteristics of fractals based on their properties.
 */

/**
 * @brief Generate a dataset of points with their fractal properties
 * 
 * @param pointCount Number of points to generate
 * @param iterationLimit Maximum iterations for Mandelbrot calculation
 * @return std::vector<std::pair<Complex, std::vector<double>>> Points and their feature vectors
 */
std::vector<std::pair<Complex, std::vector<double>>> 
generateDataset(int pointCount, int iterationLimit) {
    std::cout << "Generating dataset with " << pointCount << " points..." << std::endl;
    
    std::vector<std::pair<Complex, std::vector<double>>> dataset;
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Generate points in different regions of the complex plane
    std::uniform_real_distribution<double> realDist(-2.5, 1.0);
    std::uniform_real_distribution<double> imagDist(-1.5, 1.5);
    
    for (int i = 0; i < pointCount; i++) {
        // Generate random complex point
        double real = realDist(gen);
        double imag = imagDist(gen);
        Complex c(real, imag);
        
        // Track iteration history
        std::vector<Complex> history;
        Complex z(0.0, 0.0);
        history.push_back(z);
        
        // Iterate Mandelbrot formula
        int iterations = 0;
        while (iterations < iterationLimit && z.magnitudeSquared() <= 4.0) {
            z = z * z + c;
            history.push_back(z);
            iterations++;
        }
        
        // Extract features from iteration history
        std::vector<double> features = extractFeatures(history);
        
        // Add to dataset
        dataset.push_back({c, features});
    }
    
    std::cout << "Dataset generation complete." << std::endl;
    return dataset;
}

/**
 * @brief Train a KNN model to identify points with specific properties
 * 
 * @param dataset Training dataset
 * @param propertyFn Function that determines if a point has the target property
 * @param kNeighbors Number of neighbors for KNN
 * @return KNNModel Trained model
 */
KNNModel trainPropertyModel(
    const std::vector<std::pair<Complex, std::vector<double>>>& dataset,
    const std::function<bool(const Complex&)>& propertyFn,
    int kNeighbors = 5
) {
    std::cout << "Training KNN model with k=" << kNeighbors << "..." << std::endl;
    
    // Extract features and labels
    std::vector<std::vector<double>> features;
    std::vector<int> labels;
    
    for (const auto& [point, feature] : dataset) {
        features.push_back(feature);
        
        // Determine if the point has the target property
        bool hasProperty = propertyFn(point);
        labels.push_back(hasProperty ? 1 : 0);
    }
    
    // Create and train the model
    KNNModel model(kNeighbors);
    model.train(features, labels);
    
    // Calculate training accuracy
    int correctPredictions = 0;
    for (size_t i = 0; i < features.size(); i++) {
        int prediction = model.predict(features[i]);
        if (prediction == labels[i]) {
            correctPredictions++;
        }
    }
    
    double accuracy = 100.0 * correctPredictions / features.size();
    std::cout << "Model trained with " << features.size() << " samples." << std::endl;
    std::cout << "Training accuracy: " << accuracy << "%" << std::endl;
    
    return model;
}

/**
 * @brief Define property functions for different fractal characteristics
 * 
 * @return std::map<std::string, std::function<bool(const Complex&)>> Map of property functions
 */
std::map<std::string, std::function<bool(const Complex&)>> 
definePropertyFunctions() {
    std::map<std::string, std::function<bool(const Complex&)>> properties;
    
    // 1. Points in the main cardioid
    properties["MainCardioid"] = [](const Complex& c) -> bool {
        // p = 1/4 + 1/4 * cos(theta)
        double x = c.getReal();
        double y = c.getImag();
        
        // Check if it's in the cardioid using the quadratic formula
        double q = x - 0.25;
        double q2 = q * q;
        double y2 = y * y;
        
        return q2 + y2 < 0.25 * (q2 / q2 + y2 + 2 * q);
    };
    
    // 2. Points in the period-2 bulb
    properties["Period2Bulb"] = [](const Complex& c) -> bool {
        double x = c.getReal();
        double y = c.getImag();
        
        // Simple check for the period-2 bulb
        double xPlusOne = x + 1.0;
        return xPlusOne * xPlusOne + y * y < 0.0625;  // (x+1)² + y² < 1/16
    };
    
    // 3. Points near the boundary (simplified)
    properties["Boundary"] = [](const Complex& c) -> bool {
        // Approximate boundary check: points that escape after many iterations
        double x = c.getReal();
        double y = c.getImag();
        
        // Quick check to exclude obvious non-boundary points
        if (x < -2.0 || x > 0.5 || y < -1.2 || y > 1.2) {
            return false;
        }
        
        Complex z(0.0, 0.0);
        for (int i = 0; i < 100; i++) {
            z = z * z + c;
            
            if (z.magnitudeSquared() > 4.0) {
                // Escaped, but check if it took many iterations (likely near boundary)
                return i > 20;
            }
        }
        
        // Didn't escape, might be in the set
        return false;
    };
    
    // 4. Filaments (decorations on the main set)
    properties["Filaments"] = [](const Complex& c) -> bool {
        // Check for filaments by identifying areas that are likely in the set
        // but not in the main cardioid or period-2 bulb
        
        // Main cardioid check
        double x = c.getReal();
        double y = c.getImag();
        double q = x - 0.25;
        double q2 = q * q;
        double y2 = y * y;
        bool inCardioid = q2 + y2 < 0.25 * (q2 / q2 + y2 + 2 * q);
        
        // Period-2 bulb check
        double xPlusOne = x + 1.0;
        bool inBulb = xPlusOne * xPlusOne + y * y < 0.0625;
        
        // Check if it might be in the set but not in the main components
        if (inCardioid || inBulb) {
            return false;
        }
        
        // Iterate to see if it's in the set
        Complex z(0.0, 0.0);
        for (int i = 0; i < 300; i++) {
            z = z * z + c;
            
            if (z.magnitudeSquared() > 4.0) {
                return false;  // Not in the set
            }
        }
        
        return true;  // Likely in the set but not in main components, so a filament
    };
    
    return properties;
}

/**
 * @brief Scan a region and identify points with specific properties
 * 
 * @param model KNN model trained for the property
 * @param xMin Minimum x value
 * @param xMax Maximum x value
 * @param yMin Minimum y value
 * @param yMax Maximum y value
 * @param resolution Number of points to scan in each dimension
 * @return std::vector<Complex> Points with the predicted property
 */
std::vector<std::pair<Complex, double>> 
scanRegionForProperty(
    const KNNModel& model,
    double xMin,
    double xMax,
    double yMin,
    double yMax,
    int resolution
) {
    std::cout << "Scanning region [" << xMin << "," << xMax << "] x ["
              << yMin << "," << yMax << "] with resolution " << resolution << "..." << std::endl;
    
    std::vector<std::pair<Complex, double>> interestingPoints;
    
    // Calculate step sizes
    double xStep = (xMax - xMin) / resolution;
    double yStep = (yMax - yMin) / resolution;
    
    // Scan the region
    int totalPoints = 0;
    int interestingCount = 0;
    
    for (int y = 0; y < resolution; y++) {
        double imag = yMin + y * yStep;
        
        for (int x = 0; x < resolution; x++) {
            double real = xMin + x * xStep;
            Complex c(real, imag);
            
            // Generate iteration history
            std::vector<Complex> history;
            Complex z(0.0, 0.0);
            history.push_back(z);
            
            // Iterate Mandelbrot formula with limited iterations for efficiency
            for (int i = 0; i < 50; i++) {
                z = z * z + c;
                history.push_back(z);
                
                // Early escape check
                if (z.magnitudeSquared() > 4.0) {
                    break;
                }
            }
            
            // Extract features
            std::vector<double> features = extractFeatures(history);
            
            // Predict using the model
            int prediction = model.predict(features);
            totalPoints++;
            
            // If the point has the property, add it to the result
            if (prediction == 1) {
                interestingPoints.push_back({c, 1.0});  // Confidence is simply 1.0 for KNN
                interestingCount++;
            }
        }
    }
    
    std::cout << "Scan complete. Found " << interestingCount << " interesting points out of "
              << totalPoints << " (" << (100.0 * interestingCount / totalPoints) << "%)" << std::endl;
    
    return interestingPoints;
}

/**
 * @brief Write points to a CSV file for visualization
 * 
 * @param points Vector of points with their confidence scores
 * @param filename Output filename
 * @param propertyName Name of the property
 */
void writePointsToCSV(
    const std::vector<std::pair<Complex, double>>& points,
    const std::string& filename,
    const std::string& propertyName
) {
    std::cout << "Writing " << points.size() << " points to " << filename << "..." << std::endl;
    
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    // Write header
    outFile << "Real,Imaginary,Confidence,Property" << std::endl;
    
    // Write data
    for (const auto& [point, confidence] : points) {
        outFile << point.getReal() << "," << point.getImag() << ","
                << confidence << "," << propertyName << std::endl;
    }
    
    outFile.close();
    std::cout << "File written successfully." << std::endl;
}

/**
 * @brief Demonstrate property prediction on a validation set
 * 
 * @param model Trained model
 * @param validationSet Validation dataset
 * @param propertyFn Function that determines the ground truth
 * @param propertyName Name of the property
 */
void evaluateModelOnValidation(
    const KNNModel& model,
    const std::vector<std::pair<Complex, std::vector<double>>>& validationSet,
    const std::function<bool(const Complex&)>& propertyFn,
    const std::string& propertyName
) {
    std::cout << std::endl;
    std::cout << "=== Evaluating " << propertyName << " Model on Validation Set ===" << std::endl;
    
    int truePositives = 0;
    int falsePositives = 0;
    int trueNegatives = 0;
    int falseNegatives = 0;
    
    for (const auto& [point, features] : validationSet) {
        int prediction = model.predict(features);
        bool groundTruth = propertyFn(point);
        
        if (prediction == 1 && groundTruth) {
            truePositives++;
        } else if (prediction == 1 && !groundTruth) {
            falsePositives++;
        } else if (prediction == 0 && !groundTruth) {
            trueNegatives++;
        } else if (prediction == 0 && groundTruth) {
            falseNegatives++;
        }
    }
    
    int total = validationSet.size();
    double accuracy = 100.0 * (truePositives + trueNegatives) / total;
    double precision = truePositives > 0 ? 
        100.0 * truePositives / (truePositives + falsePositives) : 0.0;
    double recall = truePositives > 0 ? 
        100.0 * truePositives / (truePositives + falseNegatives) : 0.0;
    
    std::cout << "Validation Results:" << std::endl;
    std::cout << "  Total samples: " << total << std::endl;
    std::cout << "  Accuracy: " << accuracy << "%" << std::endl;
    std::cout << "  Precision: " << precision << "%" << std::endl;
    std::cout << "  Recall: " << recall << "%" << std::endl;
    std::cout << "  True Positives: " << truePositives << std::endl;
    std::cout << "  False Positives: " << falsePositives << std::endl;
    std::cout << "  True Negatives: " << trueNegatives << std::endl;
    std::cout << "  False Negatives: " << falseNegatives << std::endl;
}

/**
 * @brief Recommend regions for exploration based on predicted properties
 * 
 * @param interestingPoints Vector of points with their confidence scores
 * @param regionCount Number of regions to recommend
 * @param regionSize Size of each recommended region
 * @return std::vector<std::pair<Complex, double>> Recommended region centers with confidence scores
 */
std::vector<std::pair<Complex, double>> 
recommendRegions(
    const std::vector<std::pair<Complex, double>>& interestingPoints,
    int regionCount,
    double regionSize
) {
    std::cout << std::endl;
    std::cout << "=== Recommending Regions for Exploration ===" << std::endl;
    
    if (interestingPoints.empty()) {
        std::cout << "No interesting points found. Cannot recommend regions." << std::endl;
        return {};
    }
    
    // Group nearby points by clustering
    std::vector<std::pair<Complex, double>> regionCenters;
    std::vector<bool> assigned(interestingPoints.size(), false);
    
    // Simple greedy clustering
    for (int i = 0; i < regionCount && i < static_cast<int>(interestingPoints.size()); i++) {
        // Find the unassigned point with highest confidence
        int bestIdx = -1;
        double bestConfidence = -1.0;
        
        for (size_t j = 0; j < interestingPoints.size(); j++) {
            if (!assigned[j] && interestingPoints[j].second > bestConfidence) {
                bestIdx = j;
                bestConfidence = interestingPoints[j].second;
            }
        }
        
        if (bestIdx == -1) break;  // No more unassigned points
        
        // Mark this point as assigned
        assigned[bestIdx] = true;
        
        // Find all nearby points
        std::vector<size_t> clusterIndices;
        clusterIndices.push_back(bestIdx);
        
        for (size_t j = 0; j < interestingPoints.size(); j++) {
            if (!assigned[j]) {
                // Calculate distance
                Complex diff = interestingPoints[bestIdx].first + 
                              interestingPoints[j].first * Complex(-1, 0);
                double distance = diff.magnitude();
                
                // If close enough, add to cluster
                if (distance < regionSize) {
                    clusterIndices.push_back(j);
                    assigned[j] = true;
                }
            }
        }
        
        // Calculate cluster center and average confidence
        Complex center(0, 0);
        double totalConfidence = 0.0;
        
        for (size_t idx : clusterIndices) {
            center = center + interestingPoints[idx].first;
            totalConfidence += interestingPoints[idx].second;
        }
        
        center = Complex(
            center.getReal() / clusterIndices.size(),
            center.getImag() / clusterIndices.size()
        );
        
        double avgConfidence = totalConfidence / clusterIndices.size();
        
        // Add to recommended regions
        regionCenters.push_back({center, avgConfidence});
        
        std::cout << "Region " << (i+1) << ": Center = (" << center.getReal() << ", " << center.getImag() << ")" << std::endl;
        std::cout << "  Size: " << regionSize << " x " << regionSize << std::endl;
        std::cout << "  Points in region: " << clusterIndices.size() << std::endl;
        std::cout << "  Confidence: " << avgConfidence << std::endl;
    }
    
    return regionCenters;
}

/**
 * @brief Main function to demonstrate predictive capabilities
 */
int main() {
    // Generate dataset
    std::vector<std::pair<Complex, std::vector<double>>> dataset = generateDataset(5000, 100);
    
    // Split into training and validation sets
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(dataset.begin(), dataset.end(), g);
    
    size_t splitPoint = dataset.size() * 0.8;  // 80% for training, 20% for validation
    std::vector<std::pair<Complex, std::vector<double>>> trainingSet(
        dataset.begin(), dataset.begin() + splitPoint
    );
    
    std::vector<std::pair<Complex, std::vector<double>>> validationSet(
        dataset.begin() + splitPoint, dataset.end()
    );
    
    std::cout << "Dataset split: " << trainingSet.size() << " training samples, "
              << validationSet.size() << " validation samples." << std::endl;
    
    // Define property functions
    auto propertyFunctions = definePropertyFunctions();
    
    // For each property, train a model and evaluate it
    for (const auto& [propertyName, propertyFn] : propertyFunctions) {
        std::cout << std::endl;
        std::cout << "=== Training Model for " << propertyName << " ===" << std::endl;
        
        // Train the model
        KNNModel model = trainPropertyModel(trainingSet, propertyFn, 5);
        
        // Evaluate on validation set
        evaluateModelOnValidation(model, validationSet, propertyFn, propertyName);
        
        // Scan a region for this property
        double xMin, xMax, yMin, yMax;
        
        // Choose appropriate region based on property
        if (propertyName == "MainCardioid") {
            xMin = -0.8; xMax = 0.3; yMin = -0.5; yMax = 0.5;
        } else if (propertyName == "Period2Bulb") {
            xMin = -1.2; xMax = -0.8; yMin = -0.2; yMax = 0.2;
        } else if (propertyName == "Boundary") {
            xMin = -1.8; xMax = -1.7; yMin = -0.05; yMax = 0.05;
        } else if (propertyName == "Filaments") {
            xMin = -1.8; xMax = -1.7; yMin = -0.05; yMax = 0.05;
        } else {
            xMin = -2.0; xMax = 0.5; yMin = -1.2; yMax = 1.2;
        }
        
        auto interestingPoints = scanRegionForProperty(
            model, xMin, xMax, yMin, yMax, 50
        );
        
        // Write points to CSV
        writePointsToCSV(
            interestingPoints, 
            propertyName + "_points.csv", 
            propertyName
        );
        
        // Recommend regions for exploration
        auto recommendedRegions = recommendRegions(interestingPoints, 3, 0.1);
        
        // Write recommended regions to CSV
        writePointsToCSV(
            recommendedRegions,
            propertyName + "_regions.csv",
            propertyName + "_Region"
        );
    }
    
    return 0;
}