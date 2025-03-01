#include <iostream>
#include <string>
#include <vector>
#include "../Base Class/Properties/Complex Number Handling (z, c).h"
#include "../Base Class/Properties/Iteration Parameters.h"
#include "../Base Class/Methods/Iterate Equation.cpp"
#include "../Base Class/Methods/Check Divergence.cpp"
#include "../Visualization Class/Resolution Settings.cpp"
#include "../Visualization Class/Color Assignment Logic.cpp"
#include "../Visualization Class/Render Output.cpp"
#include "../Inherited Classes/GenAI Class/Integration with Fractal Computation.cpp"
#include "../Inherited Classes/GenAI Class/ML Prediction Module.cpp"
#include "../Inherited Classes/GenAI Class/Pattern Generation Logic.cpp"

/**
 * @brief Generate visual representations of Mandelbrot and related fractals
 * 
 * This module provides convenient functions to generate visual outputs
 * from the fractal computation.
 */

/**
 * @brief Generate a standard Mandelbrot set image
 * 
 * @param outputPath Path to save the output file
 * @param width Image width in pixels
 * @param height Image height in pixels
 * @param maxIterations Maximum number of iterations
 * @param region Optional region to render (default is the main cardioid view)
 * @param format Output image format
 * @return std::string Path to the generated image
 */
std::string generateMandelbrotImage(
    int width = 800,
    int height = 600,
    int maxIterations = 1000,
    ComplexRegion region = ComplexRegion(-2.5, 1.0, -1.5, 1.5),
    ImageFormat format = ImageFormat::PPM
) {
    // Create render settings
    Resolution resolution(width, height);
    IterationParameters iterParams(maxIterations, 2.0);
    RenderSettings settings(
        resolution,
        region,
        iterParams,
        ColorMode::LOGARITHMIC,
        format
    );
    
    // Render the Mandelbrot set
    return renderMandelbrotSet(settings);
}

/**
 * @brief Generate a Julia set image for a given parameter
 * 
 * @param c Complex parameter for the Julia set
 * @param width Image width in pixels
 * @param height Image height in pixels
 * @param maxIterations Maximum number of iterations
 * @param region Optional region to render (default is centered at origin)
 * @param format Output image format
 * @return std::string Path to the generated image
 */
std::string generateJuliaImage(
    const Complex& c,
    int width = 800,
    int height = 600,
    int maxIterations = 1000,
    ComplexRegion region = ComplexRegion(-2.0, 2.0, -2.0, 2.0),
    ImageFormat format = ImageFormat::PPM
) {
    // Create render settings
    Resolution resolution(width, height);
    IterationParameters iterParams(maxIterations, 2.0);
    RenderSettings settings(
        resolution,
        region,
        iterParams,
        ColorMode::SMOOTH,
        format
    );
    
    // Create a custom computation function for Julia sets
    auto computeJulia = [&c](const Complex& z, int maxIter, double threshold) -> int {
        Complex currentZ = z;
        int iterations = 0;
        double thresholdSquared = threshold * threshold;
        
        while (iterations < maxIter && currentZ.magnitudeSquared() <= thresholdSquared) {
            // Julia set formula: z = z² + c (with fixed c)
            currentZ = currentZ * currentZ + c;
            iterations++;
        }
        
        return iterations;
    };
    
    // Compute the Julia set
    const Resolution& res = settings.resolution;
    const ComplexRegion& reg = settings.region;
    const IterationParameters& params = settings.iterParams;
    
    // Initialize result array
    std::vector<std::vector<int>> iterationData(
        res.height, std::vector<int>(res.width, 0)
    );
    
    // Compute Julia set for each pixel
    for (int y = 0; y < res.height; y++) {
        for (int x = 0; x < res.width; x++) {
            // Map pixel to complex plane
            auto [real, imag] = pixelToComplex(res, reg, x, y);
            
            // Create complex number for this point
            Complex z(real, imag);
            
            // Compute iterations using the Julia function
            iterationData[y][x] = computeJulia(z, params.maxIterations, params.threshold);
        }
    }
    
    // Generate filename
    std::string filename = generateFilename(settings.format);
    
    // Augment filename to indicate Julia set
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        filename.insert(dotPos, "_julia");
    }
    
    // Write to file based on format
    switch (settings.format) {
        case ImageFormat::PPM:
            writePPMImage(filename, iterationData, settings);
            break;
        case ImageFormat::BMP:
            writeBMPImage(filename, iterationData, settings);
            break;
        case ImageFormat::PNG:
        case ImageFormat::JPG:
            throw std::runtime_error(
                "PNG and JPG formats require external libraries and are not implemented in this example"
            );
    }
    
    return filename;
}

/**
 * @brief Generate an AI-enhanced fractal image
 * 
 * @param patternMode Fractal pattern mode
 * @param width Image width in pixels
 * @param height Image height in pixels
 * @param seed Random seed for reproducibility
 * @param format Output image format
 * @return std::string Path to the generated image
 */
std::string generateAIEnhancedFractal(
    PatternMode patternMode = PatternMode::MANDELBROT,
    int width = 800,
    int height = 600,
    unsigned int seed = 42,
    ImageFormat format = ImageFormat::PPM
) {
    // Create render settings
    Resolution resolution(width, height);
    IterationParameters iterParams(1000, 2.0);
    RenderSettings renderSettings(
        resolution,
        ComplexRegion(-2.5, 1.0, -1.5, 1.5),  // Will be overridden by AI
        iterParams,
        ColorMode::LOGARITHMIC,  // Will be overridden by AI
        format
    );
    
    // Create GenAI settings
    GenAISettings genSettings(
        patternMode,
        MLModelType::KNN,
        seed,
        true,   // Use ML prediction
        false,  // Don't generate multiple sets
        1       // Number of sets
    );
    
    // Generate the AI-enhanced fractal
    return generateAIFractal(genSettings, renderSettings);
}

/**
 * @brief Generate a sequence of AI-enhanced fractals with different patterns
 * 
 * @param width Image width in pixels
 * @param height Image height in pixels
 * @param baseFileName Base filename for output
 * @return std::vector<std::string> Paths to the generated images
 */
std::vector<std::string> generateFractalSequence(
    int width = 800,
    int height = 600
) {
    std::vector<std::string> filenames;
    
    // Create render settings
    Resolution resolution(width, height);
    IterationParameters iterParams(1000, 2.0);
    RenderSettings renderSettings(
        resolution,
        ComplexRegion(-2.5, 1.0, -1.5, 1.5),
        iterParams,
        ColorMode::LOGARITHMIC,
        ImageFormat::PPM
    );
    
    // Generate different pattern types
    PatternMode patterns[] = {
        PatternMode::MANDELBROT,
        PatternMode::JULIA,
        PatternMode::BURNING_SHIP,
        PatternMode::MULTIBROT,
        PatternMode::HYBRID,
        PatternMode::RANDOMIZED
    };
    
    for (const auto& pattern : patterns) {
        // Create GenAI settings for this pattern
        GenAISettings genSettings(
            pattern,
            MLModelType::KNN,
            42 + static_cast<int>(pattern),  // Different seed for each pattern
            true,   // Use ML prediction
            false,  // Don't generate multiple sets
            1       // Number of sets
        );
        
        // Generate the fractal
        std::string filename = generateAIFractal(genSettings, renderSettings);
        filenames.push_back(filename);
    }
    
    return filenames;
}