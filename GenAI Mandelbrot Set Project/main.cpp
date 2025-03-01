#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <memory>
#include <algorithm>

// Include base components
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Base Class/Properties/Iteration Parameters.h"

// Include visualization components
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Visualization Class/Color Assignment Logic.cpp"
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Visualization Class/Resolution Settings.cpp"
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Visualization Class/Render Output.cpp"

// Include main program components
#include "GenAI Mandelbrot Set Project/Code Structure (C++)/Main Program/Initialize Fractal Parameters.cpp"

// Include output components
#include "GenAI Mandelbrot Set Project/Output/Image File/Format (e.g., PNG, BMP).cpp"

/**
 * @brief Main application entry point
 * 
 * This integrates various components of the AlmndBrd-Mandelbrot project
 * to provide a comprehensive fractal exploration tool.
 */
int main(int argc, char* argv[]) {
    std::cout << "AlmndBrd-Mandelbrot Explorer" << std::endl;
    std::cout << "Version 0.1.0" << std::endl;
    std::cout << "================================" << std::endl << std::endl;
    
    // Default parameters
    double xCenter = -0.75;
    double yCenter = 0.0;
    double zoom = 1.0;
    int width = 800;
    int height = 600;
    int maxIterations = 500;
    std::string outputFilename = "mandelbrot_output";
    
    // Parse command line arguments (simple version)
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --x=VALUE       Center x-coordinate (default: -0.75)" << std::endl;
            std::cout << "  --y=VALUE       Center y-coordinate (default: 0.0)" << std::endl;
            std::cout << "  --zoom=VALUE    Zoom level (default: 1.0)" << std::endl;
            std::cout << "  --width=VALUE   Image width (default: 800)" << std::endl;
            std::cout << "  --height=VALUE  Image height (default: 600)" << std::endl;
            std::cout << "  --iter=VALUE    Maximum iterations (default: 500)" << std::endl;
            std::cout << "  --out=FILENAME  Output filename (default: mandelbrot_output)" << std::endl;
            std::cout << "  --preset=NAME   Use a predefined preset (full, cardioid, seahorse, spiral)" << std::endl;
            return 0;
        }
        
        // Parse key=value arguments
        size_t equalPos = arg.find('=');
        if (equalPos != std::string::npos) {
            std::string key = arg.substr(0, equalPos);
            std::string value = arg.substr(equalPos + 1);
            
            if (key == "--x") xCenter = std::stod(value);
            else if (key == "--y") yCenter = std::stod(value);
            else if (key == "--zoom") zoom = std::stod(value);
            else if (key == "--width") width = std::stoi(value);
            else if (key == "--height") height = std::stoi(value);
            else if (key == "--iter") maxIterations = std::stoi(value);
            else if (key == "--out") outputFilename = value;
            else if (key == "--preset") {
                // Handle presets
                if (value == "full") {
                    xCenter = -0.5;
                    yCenter = 0.0;
                    zoom = 1.0;
                } else if (value == "cardioid") {
                    xCenter = -0.75;
                    yCenter = 0.0;
                    zoom = 4.0;
                } else if (value == "seahorse") {
                    xCenter = -0.75;
                    yCenter = -0.1;
                    zoom = 8.0;
                } else if (value == "spiral") {
                    xCenter = -0.745;
                    yCenter = 0.105;
                    zoom = 150.0;
                }
            }
        }
    }
    
    // Calculate region from center and zoom
    double regionWidth = 3.0 / zoom;
    double regionHeight = 2.0 / zoom;
    double xMin = xCenter - regionWidth / 2;
    double xMax = xCenter + regionWidth / 2;
    double yMin = yCenter - regionHeight / 2;
    double yMax = yCenter + regionHeight / 2;
    
    // Create region
    ComplexRegion region(xMin, xMax, yMin, yMax);
    
    // Create resolution
    Resolution resolution(width, height);
    
    // Create iteration parameters
    IterationParameters iterParams(maxIterations, 2.0);
    
    // Display settings
    std::cout << "Rendering Settings:" << std::endl;
    std::cout << "  Region: [" << xMin << ", " << xMax << "] x [" << yMin << ", " << yMax << "]" << std::endl;
    std::cout << "  Resolution: " << width << "x" << height << std::endl;
    std::cout << "  Maximum Iterations: " << maxIterations << std::endl;
    std::cout << "  Output Filename: " << outputFilename << std::endl;
    std::cout << std::endl;
    
    // Create render settings
    RenderSettings settings = createRenderSettings(
        region, resolution, iterParams, ColorMode::LOGARITHMIC, ImageFormatType::PPM
    );
    
    std::cout << "Calculating fractal..." << std::endl;
    
    // Measure execution time
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Initialize result array
    std::vector<std::vector<int>> iterationData(
        height, std::vector<int>(width, 0)
    );
    
    // Calculate step sizes
    double xStep = (region.xMax - region.xMin) / width;
    double yStep = (region.yMax - region.yMin) / height;
    
    // Compute Mandelbrot set for each pixel
    #pragma omp parallel for collapse(2) if(OpenMP_CXX_FOUND)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Map pixel to complex plane
            double real = region.xMin + x * xStep;
            double imag = region.yMax - y * yStep;
            
            // Create complex number for this point
            Complex c(real, imag);
            Complex z(0.0, 0.0);
            
            // Iterate until divergence or max iterations
            int iterations = 0;
            double thresholdSquared = 4.0;  // 2.0^2
            
            while (iterations < maxIterations && z.magnitudeSquared() <= thresholdSquared) {
                z = z * z + c;
                iterations++;
            }
            
            iterationData[y][x] = iterations;
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    std::cout << "Calculation completed in " << duration.count() << " ms" << std::endl;
    
    // Convert iteration data to colored pixel data
    std::cout << "Coloring fractal..." << std::endl;
    
    std::vector<std::vector<Color>> coloredData(height, std::vector<Color>(width));
    std::vector<Color> palette = generateColorPalette(ColorMode::LOGARITHMIC, maxIterations);
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iterations = iterationData[y][x];
            
            if (iterations == maxIterations) {
                // Point is in the set
                coloredData[y][x] = Color(0, 0, 0);
            } else {
                // Point is outside the set, color based on iterations
                coloredData[y][x] = palette[iterations];
            }
        }
    }
    
    // Save to file
    std::string filename = outputFilename + ".ppm";
    std::cout << "Saving to " << filename << "..." << std::endl;
    
    bool saved = saveColoredImagePPM(filename, coloredData);
    
    if (saved) {
        std::cout << "Image saved successfully." << std::endl;
    } else {
        std::cerr << "Error saving image." << std::endl;
        return 1;
    }
    
    std::cout << "Done." << std::endl;
    return 0;
}