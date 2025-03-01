#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include "../../Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Types of machine learning models that can be used for fractal predictions
 */
enum class MLModelType {
    KNN,            // K-Nearest Neighbors
    RANDOM_FOREST,  // Random Forest
    NEURAL_NET,     // Neural Network
    HYBRID,         // Hybrid approach
    EVOLUTIONARY    // Evolutionary algorithm
};

/**
 * @brief Feature extraction from fractal iterations
 * 
 * Extracts numerical features from the iteration history of a point
 * in the Mandelbrot set computation
 * 
 * @param iterationHistory Vector of z values at each iteration
 * @return std::vector<double> Extracted features
 */
std::vector<double> extractFeatures(const std::vector<Complex>& iterationHistory) {
    std::vector<double> features;
    
    // If history is empty, return empty feature vector
    if (iterationHistory.empty()) {
        return features;
    }
    
    // Calculate statistical features from the iteration history
    
    // 1. Trajectory length (sum of distances between consecutive points)
    double trajectoryLength = 0.0;
    for (size_t i = 1; i < iterationHistory.size(); i++) {
        Complex diff = iterationHistory[i] + iterationHistory[i-1] * Complex(-1, 0);
        trajectoryLength += diff.magnitude();
    }
    features.push_back(trajectoryLength);
    
    // 2. Final magnitude
    features.push_back(iterationHistory.back().magnitude());
    
    // 3. Average real component
    double avgReal = 0.0;
    for (const auto& z : iterationHistory) {
        avgReal += z.getReal();
    }
    avgReal /= iterationHistory.size();
    features.push_back(avgReal);
    
    // 4. Average imaginary component
    double avgImag = 0.0;
    for (const auto& z : iterationHistory) {
        avgImag += z.getImag();
    }
    avgImag /= iterationHistory.size();
    features.push_back(avgImag);
    
    // 5. Standard deviation of magnitudes
    double avgMag = 0.0;
    for (const auto& z : iterationHistory) {
        avgMag += z.magnitude();
    }
    avgMag /= iterationHistory.size();
    
    double varMag = 0.0;
    for (const auto& z : iterationHistory) {
        double diff = z.magnitude() - avgMag;
        varMag += diff * diff;
    }
    varMag /= iterationHistory.size();
    features.push_back(std::sqrt(varMag));
    
    // 6. Growth rate (ratio of final to initial magnitude)
    if (iterationHistory.front().magnitude() > 1e-10) {
        features.push_back(iterationHistory.back().magnitude() / iterationHistory.front().magnitude());
    } else {
        features.push_back(iterationHistory.back().magnitude() / 1e-10);
    }
    
    return features;
}

/**
 * @brief Simple KNN model for prediction
 * 
 * This is a simplified simulation of a K-Nearest Neighbors model
 * for educational purposes.
 */
class KNNModel {
private:
    int k;
    std::vector<std::vector<double>> trainingFeatures;
    std::vector<int> trainingLabels;
    
public:
    /**
     * @brief Construct a new KNNModel
     * 
     * @param kNeighbors Number of neighbors to consider
     */
    KNNModel(int kNeighbors = 3) : k(kNeighbors) {}
    
    /**
     * @brief Train the model with features and labels
     * 
     * @param features Vector of feature vectors
     * @param labels Vector of corresponding labels
     */
    void train(const std::vector<std::vector<double>>& features, const std::vector<int>& labels) {
        if (features.size() != labels.size()) {
            throw std::invalid_argument("Features and labels must have the same size");
        }
        
        trainingFeatures = features;
        trainingLabels = labels;
    }
    
    /**
     * @brief Predict the class of a new feature vector
     * 
     * @param features Feature vector to classify
     * @return int Predicted class
     */
    int predict(const std::vector<double>& features) const {
        if (trainingFeatures.empty()) {
            throw std::runtime_error("Model has not been trained");
        }
        
        // Calculate distances to all training points
        std::vector<std::pair<double, int>> distances;
        for (size_t i = 0; i < trainingFeatures.size(); i++) {
            double distance = euclideanDistance(features, trainingFeatures[i]);
            distances.push_back({distance, trainingLabels[i]});
        }
        
        // Sort distances
        std::sort(distances.begin(), distances.end());
        
        // Count votes for the k nearest neighbors
        std::unordered_map<int, int> votes;
        for (int i = 0; i < std::min(k, static_cast<int>(distances.size())); i++) {
            votes[distances[i].second]++;
        }
        
        // Find the class with most votes
        int bestClass = -1;
        int maxVotes = 0;
        for (const auto& [cls, count] : votes) {
            if (count > maxVotes) {
                maxVotes = count;
                bestClass = cls;
            }
        }
        
        return bestClass;
    }
    
private:
    /**
     * @brief Calculate Euclidean distance between two feature vectors
     * 
     * @param a First feature vector
     * @param b Second feature vector
     * @return double Euclidean distance
     */
    double euclideanDistance(const std::vector<double>& a, const std::vector<double>& b) const {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Feature vectors must have the same size");
        }
        
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        
        return std::sqrt(sum);
    }
};

/**
 * @brief Simple neural network for prediction (simulation)
 * 
 * This is a simplified simulation of a neural network model
 * for educational purposes.
 */
class NeuralNetworkModel {
private:
    // Simulated weights for a very simple network
    std::vector<double> weights;
    double bias;
    
    // Random number generator
    std::mt19937 rng;
    
public:
    /**
     * @brief Construct a new NeuralNetworkModel
     * 
     * @param seed Random seed for initialization
     */
    NeuralNetworkModel(unsigned int seed = 42) : bias(0.0) {
        rng = std::mt19937(seed);
    }
    
    /**
     * @brief Initialize the model with random weights
     * 
     * @param inputSize Number of input features
     */
    void initialize(size_t inputSize) {
        // Simulate Xavier initialization
        std::uniform_real_distribution<double> dist(-1.0 / std::sqrt(inputSize), 1.0 / std::sqrt(inputSize));
        
        weights.resize(inputSize);
        for (auto& w : weights) {
            w = dist(rng);
        }
        
        bias = dist(rng);
    }
    
    /**
     * @brief Train the model with features and labels (simulation)
     * 
     * @param features Vector of feature vectors
     * @param labels Vector of corresponding labels
     * @param epochs Number of training epochs
     * @param learningRate Learning rate
     */
    void train(
        const std::vector<std::vector<double>>& features,
        const std::vector<int>& labels,
        int epochs = 100,
        double learningRate = 0.01
    ) {
        if (features.empty()) {
            throw std::invalid_argument("Features cannot be empty");
        }
        
        if (features.size() != labels.size()) {
            throw std::invalid_argument("Features and labels must have the same size");
        }
        
        // Initialize if not done yet
        if (weights.empty()) {
            initialize(features[0].size());
        }
        
        // Simple stochastic gradient descent (simulated)
        std::uniform_int_distribution<size_t> indexDist(0, features.size() - 1);
        
        for (int epoch = 0; epoch < epochs; epoch++) {
            // Randomly select a training example
            size_t idx = indexDist(rng);
            const auto& x = features[idx];
            int y = labels[idx];
            
            // Forward pass
            double output = predict(x, false);
            
            // Gradient (simplistic)
            double error = output - y;
            
            // Update weights and bias
            for (size_t i = 0; i < weights.size(); i++) {
                weights[i] -= learningRate * error * x[i];
            }
            bias -= learningRate * error;
        }
    }
    
    /**
     * @brief Predict the output for a feature vector
     * 
     * @param features Feature vector
     * @param binarize Whether to binarize the output
     * @return double Predicted value
     */
    double predict(const std::vector<double>& features, bool binarize = true) const {
        if (features.size() != weights.size()) {
            throw std::invalid_argument("Feature vector size does not match model");
        }
        
        // Weighted sum
        double sum = bias;
        for (size_t i = 0; i < features.size(); i++) {
            sum += features[i] * weights[i];
        }
        
        // Apply activation function (sigmoid)
        double output = 1.0 / (1.0 + std::exp(-sum));
        
        // Binarize if requested
        if (binarize) {
            return output > 0.5 ? 1.0 : 0.0;
        }
        
        return output;
    }
};

/**
 * @brief Generate training data for ML models using Mandelbrot characteristics
 * 
 * @param numSamples Number of samples to generate
 * @param iterationLimit Maximum iterations for Mandelbrot calculation
 * @return std::pair<std::vector<std::vector<double>>, std::vector<int>> 
 *         Pair of features and labels
 */
std::pair<std::vector<std::vector<double>>, std::vector<int>> 
generateTrainingData(int numSamples, int iterationLimit) {
    std::vector<std::vector<double>> features;
    std::vector<int> labels;
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Generate points in different regions of the complex plane
    std::uniform_real_distribution<double> realDist(-2.5, 1.5);
    std::uniform_real_distribution<double> imagDist(-1.5, 1.5);
    
    for (int i = 0; i < numSamples; i++) {
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
        std::vector<double> featureVector = extractFeatures(history);
        
        // Label: 1 if point is likely in the set, 0 if it escaped
        int label = (iterations == iterationLimit) ? 1 : 0;
        
        features.push_back(featureVector);
        labels.push_back(label);
    }
    
    return {features, labels};
}

/**
 * @brief Predict interesting regions using machine learning
 * 
 * @param modelType Type of ML model to use
 * @param trainingSize Number of training samples
 * @param regionSize Size of region to analyze
 * @param resolution Resolution for analysis
 * @return std::vector<std::pair<Complex, double>> 
 *         Vector of points and their interestingness scores
 */
std::vector<std::pair<Complex, double>> 
predictInterestingRegions(
    MLModelType modelType,
    int trainingSize = 1000,
    double regionSize = 4.0,
    int resolution = 100
) {
    // Generate training data
    auto [trainingFeatures, trainingLabels] = generateTrainingData(trainingSize, 100);
    
    // Create appropriate model based on type
    std::function<double(const std::vector<double>&)> predictFn;
    
    if (modelType == MLModelType::KNN) {
        // Train KNN model
        KNNModel knn(5);
        knn.train(trainingFeatures, trainingLabels);
        
        predictFn = [&knn](const std::vector<double>& features) {
            return static_cast<double>(knn.predict(features));
        };
    }
    else if (modelType == MLModelType::NEURAL_NET) {
        // Train neural network model
        NeuralNetworkModel nn;
        nn.train(trainingFeatures, trainingLabels, 500, 0.01);
        
        predictFn = [&nn](const std::vector<double>& features) {
            return nn.predict(features, false);  // Get raw probability
        };
    }
    else {
        // Default to a random forest-like ensemble (simulated)
        predictFn = [&trainingFeatures, &trainingLabels](const std::vector<double>& features) {
            // Simplified ensemble prediction (average of 3 KNN models with different k)
            KNNModel knn1(3);
            KNNModel knn2(5);
            KNNModel knn3(7);
            
            knn1.train(trainingFeatures, trainingLabels);
            knn2.train(trainingFeatures, trainingLabels);
            knn3.train(trainingFeatures, trainingLabels);
            
            double pred1 = static_cast<double>(knn1.predict(features));
            double pred2 = static_cast<double>(knn2.predict(features));
            double pred3 = static_cast<double>(knn3.predict(features));
            
            return (pred1 + pred2 + pred3) / 3.0;
        };
    }
    
    // Generate grid of points to analyze
    std::vector<std::pair<Complex, double>> interestingPoints;
    double step = regionSize / resolution;
    
    for (int y = 0; y < resolution; y++) {
        double imag = -regionSize/2 + y * step;
        
        for (int x = 0; x < resolution; x++) {
            double real = -regionSize/2 + x * step;
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
            
            // Predict interestingness score
            double score = predictFn(features);
            
            // Store point and score if interesting enough
            if (score > 0.6) {  // Threshold for "interesting"
                interestingPoints.push_back({c, score});
            }
        }
    }
    
    // Sort points by interestingness score (descending)
    std::sort(interestingPoints.begin(), interestingPoints.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    return interestingPoints;
}