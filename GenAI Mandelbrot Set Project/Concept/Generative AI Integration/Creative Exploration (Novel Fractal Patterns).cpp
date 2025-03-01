#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <chrono>
#include <functional>
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"
#include "../../Code Structure (C++)/Inherited Classes/GenAI Class/Pattern Generation Logic.cpp"

/**
 * @brief Creative exploration of novel fractal patterns
 * 
 * This file demonstrates how to use generative techniques to discover
 * and explore novel fractal patterns beyond the standard Mandelbrot set.
 */

/**
 * @brief Simple implementation of a random seed generator
 * 
 * @return unsigned int Random seed
 */
unsigned int generateRandomSeed() {
    return static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );
}

/**
 * @brief Generate a PPM image file for a fractal
 * 
 * @param iterationData 2D array of iteration counts
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 * @param filename Output filename
 * @param colorMode Color mapping mode
 */
void generateFractalImage(
    const std::vector<std::vector<int>>& iterationData,
    int width,
    int height,
    int maxIterations,
    const std::string& filename,
    ColorMode colorMode = ColorMode::LOGARITHMIC
) {
    // Generate color palette
    std::vector<Color> palette = generateColorPalette(colorMode, maxIterations);
    
    // Open output file
    std::ofstream outFile(filename, std::ios::out);
    if (!outFile) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    // Write PPM header
    outFile << "P3\n"; // ASCII RGB format
    outFile << width << " " << height << "\n";
    outFile << "255\n"; // Max color value
    
    // Write pixel data
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iterations = iterationData[y][x];
            const Color& color = palette[iterations];
            
            outFile << static_cast<int>(color.r) << " "
                    << static_cast<int>(color.g) << " "
                    << static_cast<int>(color.b) << " ";
        }
        outFile << "\n";
    }
    
    outFile.close();
}

/**
 * @brief Generate a fractal pattern using the given iteration function
 * 
 * @param iterationFn Function that performs one iteration
 * @param xMin Minimum x value in complex plane
 * @param xMax Maximum x value in complex plane
 * @param yMin Minimum y value in complex plane
 * @param yMax Maximum y value in complex plane
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 * @param isJulia Whether this is a Julia set (fixed c, varying z0)
 * @param juliaC Julia set parameter (if isJulia is true)
 * @return std::vector<std::vector<int>> 2D array of iteration counts
 */
std::vector<std::vector<int>> generateFractalPattern(
    const std::function<Complex(Complex, Complex)>& iterationFn,
    double xMin,
    double xMax,
    double yMin,
    double yMax,
    int width,
    int height,
    int maxIterations,
    bool isJulia = false,
    Complex juliaC = Complex(0, 0)
) {
    // Initialize result array
    std::vector<std::vector<int>> iterationData(
        height, std::vector<int>(width, 0)
    );
    
    // Calculate steps
    double xStep = (xMax - xMin) / width;
    double yStep = (yMax - yMin) / height;
    
    // Threshold for escape
    double thresholdSquared = 4.0;  // |z| > 2.0
    
    // Generate fractal
    for (int y = 0; y < height; y++) {
        double imag = yMax - y * yStep;
        
        for (int x = 0; x < width; x++) {
            double real = xMin + x * xStep;
            
            Complex z, c;
            
            if (isJulia) {
                // For Julia sets, z varies and c is fixed
                z = Complex(real, imag);
                c = juliaC;
            } else {
                // For Mandelbrot-like sets, z starts at 0 and c varies
                z = Complex(0, 0);
                c = Complex(real, imag);
            }
            
            // Iterate until escape or max iterations
            int iterations = 0;
            
            while (iterations < maxIterations && z.magnitudeSquared() <= thresholdSquared) {
                z = iterationFn(z, c);
                iterations++;
            }
            
            iterationData[y][x] = iterations;
        }
    }
    
    return iterationData;
}

/**
 * @brief Explore a random variation of fractal patterns
 * 
 * @param count Number of patterns to generate
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 */
void exploreRandomPatterns(int count, int width, int height, int maxIterations) {
    std::cout << "=== Exploring Random Fractal Patterns ===" << std::endl;
    std::cout << "Generating " << count << " random patterns..." << std::endl;
    
    for (int i = 0; i < count; i++) {
        // Generate a random seed
        unsigned int seed = generateRandomSeed();
        
        // Choose a random pattern mode
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int> modeDist(0, static_cast<int>(PatternMode::RANDOMIZED));
        PatternMode mode = static_cast<PatternMode>(modeDist(gen));
        
        // Generate iteration function
        auto iterationFn = generateIterationFunction(mode, seed);
        
        // Determine if this is a Julia set
        bool isJulia = (mode == PatternMode::JULIA);
        Complex juliaC;
        
        if (isJulia) {
            juliaC = generateJuliaParameter(seed);
        }
        
        // Generate fractal pattern
        std::vector<std::vector<int>> iterationData = generateFractalPattern(
            iterationFn,
            -2.0, 2.0,   // X range
            -1.5, 1.5,   // Y range
            width, height,
            maxIterations,
            isJulia,
            juliaC
        );
        
        // Generate image filename
        std::string filename = "fractal_" + std::to_string(i) + "_" + 
                               std::to_string(seed) + ".ppm";
        
        // Choose a random color mode
        std::uniform_int_distribution<int> colorDist(0, 3);
        ColorMode colorMode = static_cast<ColorMode>(colorDist(gen));
        
        // Generate image
        generateFractalImage(iterationData, width, height, maxIterations, filename, colorMode);
        
        // Get pattern description
        std::string description = generatePatternDescription(mode, seed, juliaC);
        
        std::cout << "Generated pattern " << (i+1) << ": " << filename << std::endl;
        std::cout << "  " << description << std::endl;
    }
}

/**
 * @brief Create variations on a theme by modifying parameters
 * 
 * @param baseMode Base pattern mode
 * @param baseSeed Base random seed
 * @param variationCount Number of variations to generate
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 */
void createVariationsOnTheme(
    PatternMode baseMode,
    unsigned int baseSeed,
    int variationCount,
    int width,
    int height,
    int maxIterations
) {
    std::cout << std::endl;
    std::cout << "=== Creating Variations on a Theme ===" << std::endl;
    std::string baseDescription = generatePatternDescription(baseMode, baseSeed);
    std::cout << "Base pattern: " << baseDescription << std::endl;
    std::cout << "Generating " << variationCount << " variations..." << std::endl;
    
    // Generate the base pattern
    auto baseIterationFn = generateIterationFunction(baseMode, baseSeed);
    bool isJulia = (baseMode == PatternMode::JULIA);
    Complex baseJuliaC;
    
    if (isJulia) {
        baseJuliaC = generateJuliaParameter(baseSeed);
    }
    
    // Generate variations
    for (int i = 0; i < variationCount; i++) {
        // Create a new seed based on the base seed and variation index
        unsigned int varSeed = baseSeed + i + 1;
        
        // Choose a modified version of the base pattern
        PatternMode varMode = baseMode;
        
        // If this is not already a randomized pattern, occasionally switch to randomized
        if (baseMode != PatternMode::RANDOMIZED && i % 3 == 2) {
            varMode = PatternMode::RANDOMIZED;
        }
        
        // Generate iteration function for this variation
        auto iterationFn = generateIterationFunction(varMode, varSeed);
        
        // For Julia sets, vary the parameter slightly
        Complex juliaC = baseJuliaC;
        
        if (isJulia) {
            // Create a small variation on the Julia parameter
            std::mt19937 gen(varSeed);
            std::uniform_real_distribution<double> smallVar(-0.05, 0.05);
            juliaC = Complex(baseJuliaC.getReal() + smallVar(gen), 
                            baseJuliaC.getImag() + smallVar(gen));
        }
        
        // Choose a different region to zoom in on
        std::mt19937 gen(varSeed);
        std::uniform_real_distribution<double> zoomDist(0.5, 2.0);
        std::uniform_real_distribution<double> centerVarDist(-0.5, 0.5);
        
        double zoomFactor = zoomDist(gen);
        double centerXVar = centerVarDist(gen);
        double centerYVar = centerVarDist(gen);
        
        // Calculate region
        double xMin = -2.0 / zoomFactor + centerXVar;
        double xMax = 2.0 / zoomFactor + centerXVar;
        double yMin = -1.5 / zoomFactor + centerYVar;
        double yMax = 1.5 / zoomFactor + centerYVar;
        
        // Generate fractal pattern
        std::vector<std::vector<int>> iterationData = generateFractalPattern(
            iterationFn,
            xMin, xMax,
            yMin, yMax,
            width, height,
            maxIterations,
            isJulia,
            juliaC
        );
        
        // Generate image filename
        std::string filename = "variation_" + std::to_string(i) + "_" + 
                               std::to_string(varSeed) + ".ppm";
        
        // Choose a variation on the color scheme
        std::uniform_int_distribution<int> colorDist(0, 3);
        ColorMode colorMode = static_cast<ColorMode>(colorDist(gen));
        
        // Generate image
        generateFractalImage(iterationData, width, height, maxIterations, filename, colorMode);
        
        // Get pattern description
        std::string description = generatePatternDescription(varMode, varSeed, juliaC);
        
        std::cout << "Generated variation " << (i+1) << ": " << filename << std::endl;
        std::cout << "  " << description << std::endl;
        std::cout << "  Region: [" << xMin << "," << xMax << "] x [" << yMin << "," << yMax << "]" << std::endl;
    }
}

/**
 * @brief Explore combinations of different fractal formulas
 * 
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 */
void exploreCombinations(int width, int height, int maxIterations) {
    std::cout << std::endl;
    std::cout << "=== Exploring Combinations of Fractal Formulas ===" << std::endl;
    
    // Define several iteration functions
    auto mandelbrot = [](const Complex& z, const Complex& c) -> Complex {
        return z * z + c;
    };
    
    auto burningShip = [](const Complex& z, const Complex& c) -> Complex {
        double absReal = std::abs(z.getReal());
        double absImag = std::abs(z.getImag());
        Complex absZ(absReal, absImag);
        return absZ * absZ + c;
    };
    
    auto tricorn = [](const Complex& z, const Complex& c) -> Complex {
        Complex conjugate(z.getReal(), -z.getImag());
        return conjugate * conjugate + c;
    };
    
    auto sine = [](const Complex& z, const Complex& c) -> Complex {
        double sinReal = std::sin(z.getReal());
        double sinImag = std::sin(z.getImag());
        Complex sinZ(sinReal, sinImag);
        return sinZ + c;
    };
    
    // Create combinations of these functions
    std::vector<std::function<Complex(Complex, Complex)>> combinations;
    std::vector<std::string> descriptions;
    
    // Basic functions
    combinations.push_back(mandelbrot);
    descriptions.push_back("Standard Mandelbrot: z = z² + c");
    
    combinations.push_back(burningShip);
    descriptions.push_back("Burning Ship: z = (|Re(z)| + |Im(z)|i)² + c");
    
    combinations.push_back(tricorn);
    descriptions.push_back("Tricorn: z = (z*)² + c where z* is complex conjugate");
    
    combinations.push_back(sine);
    descriptions.push_back("Sine: z = sin(z) + c");
    
    // Combinations
    combinations.push_back([&](const Complex& z, const Complex& c) -> Complex {
        return mandelbrot(sine(z, c), c);
    });
    descriptions.push_back("Mandelbrot of Sine: z = (sin(z) + c)² + c");
    
    combinations.push_back([&](const Complex& z, const Complex& c) -> Complex {
        return burningShip(mandelbrot(z, c), c);
    });
    descriptions.push_back("Burning Ship of Mandelbrot: z = (|Re(z² + c)| + |Im(z² + c)|i)² + c");
    
    combinations.push_back([&](const Complex& z, const Complex& c) -> Complex {
        return mandelbrot(z, c) * 0.5 + burningShip(z, c) * 0.5;
    });
    descriptions.push_back("Blend 50/50: z = 0.5(z² + c) + 0.5((|Re(z)| + |Im(z)|i)² + c)");
    
    combinations.push_back([&](const Complex& z, const Complex& c) -> Complex {
        Complex temp = mandelbrot(z, c);
        double mag = temp.magnitude();
        if (mag > 0) {
            // Normalize and apply exponential factor
            temp = Complex(
                temp.getReal() / mag * std::log(mag + 1), 
                temp.getImag() / mag * std::log(mag + 1)
            );
        }
        return temp;
    });
    descriptions.push_back("Logarithmic normalization: z = log(|z² + c| + 1) * (z² + c)/|z² + c|");
    
    // Generate an image for each combination
    for (size_t i = 0; i < combinations.size(); i++) {
        std::cout << "Generating " << descriptions[i] << "..." << std::endl;
        
        // Generate fractal pattern
        std::vector<std::vector<int>> iterationData = generateFractalPattern(
            combinations[i],
            -2.0, 2.0,   // X range
            -1.5, 1.5,   // Y range
            width, height,
            maxIterations,
            false,       // Not a Julia set
            Complex(0, 0)
        );
        
        // Generate image filename
        std::string filename = "combo_" + std::to_string(i) + ".ppm";
        
        // Generate image with a consistent color scheme
        generateFractalImage(iterationData, width, height, maxIterations, filename, ColorMode::LOGARITHMIC);
        
        std::cout << "  Generated: " << filename << std::endl;
    }
}

/**
 * @brief Demonstrate the creation of animation keyframes
 * 
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 * @param frameCount Number of frames to generate
 */
void generateAnimationKeyframes(int width, int height, int maxIterations, int frameCount) {
    std::cout << std::endl;
    std::cout << "=== Generating Animation Keyframes ===" << std::endl;
    std::cout << "Creating a zoom sequence with " << frameCount << " frames..." << std::endl;
    
    // Select a pattern mode and seed
    PatternMode mode = PatternMode::MANDELBROT;
    unsigned int seed = 42;
    
    // Generate iteration function
    auto iterationFn = generateIterationFunction(mode, seed);
    
    // Define zoom target (a point of interest)
    double targetX = -0.743643887037158704752191506114774;
    double targetY = 0.131825904205311970493132056385139;
    
    // Create exponential zoom sequence
    double startZoom = 1.0;
    double endZoom = 1000.0;
    
    for (int frame = 0; frame < frameCount; frame++) {
        // Calculate current zoom factor (exponential interpolation)
        double t = static_cast<double>(frame) / (frameCount - 1);
        double zoomFactor = startZoom * std::pow(endZoom / startZoom, t);
        
        // Calculate region based on zoom
        double regionSize = 3.0 / zoomFactor;
        double xMin = targetX - regionSize/2;
        double xMax = targetX + regionSize/2;
        double yMin = targetY - regionSize * height / (2 * width);
        double yMax = targetY + regionSize * height / (2 * width);
        
        // Generate fractal pattern
        std::vector<std::vector<int>> iterationData = generateFractalPattern(
            iterationFn,
            xMin, xMax,
            yMin, yMax,
            width, height,
            maxIterations,
            false,
            Complex(0, 0)
        );
        
        // Generate image filename with padding for proper sorting
        std::ostringstream oss;
        oss << "zoom_" << std::setw(4) << std::setfill('0') << frame << ".ppm";
        std::string filename = oss.str();
        
        // Generate image
        generateFractalImage(
            iterationData, 
            width, height, 
            maxIterations, 
            filename, 
            ColorMode::SMOOTH
        );
        
        std::cout << "Generated frame " << frame + 1 << "/" << frameCount << ": "
                  << filename << " (zoom: " << zoomFactor << "x)" << std::endl;
    }
    
    std::cout << "Animation frames completed. You can convert these to a video using:" << std::endl;
    std::cout << "ffmpeg -framerate 30 -i zoom_%04d.ppm -c:v libx264 -pix_fmt yuv420p zoom_animation.mp4" << std::endl;
}

/**
 * @brief Main function to demonstrate creative exploration
 */
int main() {
    // Define parameters for the explorations
    int width = 800;
    int height = 600;
    int maxIterations = 500;
    
    // Explore random fractal patterns
    exploreRandomPatterns(3, width, height, maxIterations);
    
    // Create variations on a theme
    createVariationsOnTheme(
        PatternMode::MANDELBROT,
        42,
        3,
        width, height, maxIterations
    );
    
    // Explore combinations of different fractal formulas
    exploreCombinations(width, height, maxIterations);
    
    // Generate animation keyframes (zoom sequence)
    generateAnimationKeyframes(width, height, maxIterations, 5);
    
    return 0;
}