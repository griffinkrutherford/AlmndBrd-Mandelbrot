#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <memory>
#include "../../Code Structure (C++)/Visualization Class/Color Assignment Logic.cpp"
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @file Blacked-Out Areas (Convergent).cpp
 * @brief Implementation of visualization techniques for convergent areas in fractal sets
 * 
 * This file provides methods for visualizing the interior of fractal sets
 * (areas where the iterative sequence does not escape to infinity).
 */

/**
 * @brief Method for visualizing the interior (convergent) areas of the set
 */
enum class InteriorVisualizationMethod {
    SOLID_BLACK,            // Traditional solid black
    LEVEL_SETS,             // Banding based on distance to boundary
    POTENTIAL,              // Potential function visualization
    DISTANCE_ESTIMATION,    // Distance estimation function
    PERIOD_DETECTION,       // Period detection for interior coloring
    BINARY_DECOMPOSITION,   // Binary decomposition coloring
    NORMALIZED_ITERATIONS,  // Normalized iteration count
    CURVATURE,              // Based on curvature of series
    INTERIOR_ORBIT_TRAPS,   // Orbit traps for interior points
    EXTERNAL_RAYS           // External rays and internal regions
};

/**
 * @brief Structure for interior visualization options
 */
struct InteriorVisualizationOptions {
    InteriorVisualizationMethod method;   // Visualization method
    int numLevels;                        // Number of level sets or bands
    double contrast;                      // Contrast adjustment
    double brightness;                    // Brightness adjustment
    bool invertColors;                    // Whether to invert colors
    int colorShift;                       // Color shift amount
    bool showBoundary;                    // Whether to highlight the boundary
    
    // Constructor with defaults
    InteriorVisualizationOptions(
        InteriorVisualizationMethod m = InteriorVisualizationMethod::SOLID_BLACK,
        int levels = 8,
        double c = 1.0,
        double b = 1.0,
        bool invert = false,
        int shift = 0,
        bool boundary = true
    ) : method(m),
        numLevels(levels),
        contrast(c),
        brightness(b),
        invertColors(invert),
        colorShift(shift),
        showBoundary(boundary) {
    }
};

/**
 * @brief Calculate distance estimation for a point in the Mandelbrot set
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum iterations
 * @param threshold Escape threshold
 * @return std::pair<double, bool> Distance estimate and whether point is in set
 */
std::pair<double, bool> calculateDistanceEstimation(
    const Complex& c,
    int maxIterations,
    double threshold
) {
    Complex z(0.0, 0.0);
    Complex dz(0.0, 0.0);  // Derivative
    
    for (int i = 0; i < maxIterations; i++) {
        // Calculate derivative: dz = 2 * z * dz + 1
        dz = dz * z;
        dz = dz * 2.0 + Complex(1.0, 0.0);
        
        // Update z: z = z² + c
        z = z * z + c;
        
        // Check for escape
        if (z.magnitudeSquared() > threshold * threshold) {
            // Point is outside the set
            double zMag = z.magnitude();
            double dzMag = dz.magnitude();
            double de = 2.0 * zMag * std::log(zMag) / dzMag;
            return {de, false};
        }
    }
    
    // Point is inside the set (or at least didn't escape)
    return {0.0, true};
}

/**
 * @brief Detect period of orbit for points inside the set
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum iterations
 * @param threshold Escape threshold
 * @return std::pair<int, bool> Period and whether detection was successful
 */
std::pair<int, bool> detectPeriod(
    const Complex& c,
    int maxIterations,
    double threshold
) {
    Complex z(0.0, 0.0);
    
    // Floyd's cycle-finding algorithm (tortoise and hare)
    Complex tortoise = z;
    Complex hare = z;
    
    // Initial step to avoid trivial cycles
    tortoise = tortoise * tortoise + c;
    hare = hare * hare + c;
    hare = hare * hare + c;
    
    int period = 1;
    
    while (period < maxIterations) {
        // Check for escape
        if (tortoise.magnitudeSquared() > threshold * threshold || 
            hare.magnitudeSquared() > threshold * threshold) {
            return {0, false};  // Point escapes, not in the set
        }
        
        // Check for cycle
        Complex diff = tortoise + hare * Complex(-1, 0);
        if (diff.magnitudeSquared() < 1e-10) {
            // Cycle found, find exact period
            z = tortoise;
            int cycleLength = 1;
            Complex next = z * z + c;
            
            while ((next + z * Complex(-1, 0)).magnitudeSquared() > 1e-10) {
                z = next;
                next = z * z + c;
                cycleLength++;
                
                if (cycleLength > maxIterations) {
                    return {0, false};  // Couldn't determine period
                }
            }
            
            return {cycleLength, true};
        }
        
        // Move tortoise one step, hare two steps
        tortoise = tortoise * tortoise + c;
        hare = hare * hare + c;
        hare = hare * hare + c;
        
        period++;
    }
    
    return {0, false};  // Couldn't determine period
}

/**
 * @brief Calculate potential function for interior points
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum iterations
 * @param threshold Escape threshold
 * @return double Potential value (0.0-1.0)
 */
double calculatePotential(
    const Complex& c,
    int maxIterations,
    double threshold
) {
    Complex z(0.0, 0.0);
    
    for (int i = 0; i < maxIterations; i++) {
        z = z * z + c;
        
        if (z.magnitudeSquared() > threshold * threshold) {
            // Point escaped, calculate potential
            double zMag = z.magnitude();
            return i - std::log2(std::log(zMag));
        }
    }
    
    // Point is inside the set, use special value
    return 0.0;
}

/**
 * @brief Calculate binary decomposition value
 * 
 * @param c Complex parameter
 * @param maxIterations Maximum iterations
 * @param threshold Escape threshold
 * @return int Binary decomposition value
 */
int calculateBinaryDecomposition(
    const Complex& c,
    int maxIterations,
    double threshold
) {
    Complex z(0.0, 0.0);
    int decomposition = 0;
    int bit = 1;
    
    for (int i = 0; i < std::min(maxIterations, 20); i++) {  // Limit to 20 bits
        z = z * z + c;
        
        if (z.magnitudeSquared() > threshold * threshold) {
            return 0;  // Point escapes, not in the set
        }
        
        // Record whether real part is positive
        if (z.getReal() > 0) {
            decomposition |= bit;
        }
        
        bit <<= 1;
    }
    
    return decomposition;
}

/**
 * @brief Generate a color palette for interior visualization
 * 
 * @param numColors Number of colors
 * @param options Visualization options
 * @return std::vector<Color> Color palette
 */
std::vector<Color> generateInteriorPalette(
    int numColors,
    const InteriorVisualizationOptions& options
) {
    std::vector<Color> palette(numColors);
    
    switch (options.method) {
        case InteriorVisualizationMethod::SOLID_BLACK: {
            // Just black
            for (int i = 0; i < numColors; i++) {
                palette[i] = Color(0, 0, 0);
            }
            break;
        }
        
        case InteriorVisualizationMethod::LEVEL_SETS: {
            // Gradients between black and dark colors
            for (int i = 0; i < numColors; i++) {
                double position = static_cast<double>(i) / (numColors - 1);
                
                // Apply contrast and brightness
                position = 0.5 + (position - 0.5) * options.contrast;
                position = std::min(1.0, std::max(0.0, position));
                position = position * options.brightness;
                
                // Shift hue based on colorShift
                double hue = (position * 360.0 + options.colorShift) / options.numLevels;
                hue = std::fmod(hue, 360.0);
                
                // Get dark color based on hue
                double saturation = 0.7;
                double value = 0.3 * position;  // Keep it dark
                
                if (options.invertColors) {
                    value = 0.3 * (1.0 - position);
                }
                
                // Convert HSV to RGB
                double c = value * saturation;
                double x = c * (1 - std::abs(std::fmod(hue / 60.0, 2) - 1));
                double m = value - c;
                
                double r, g, b;
                if (hue < 60) {
                    r = c; g = x; b = 0;
                } else if (hue < 120) {
                    r = x; g = c; b = 0;
                } else if (hue < 180) {
                    r = 0; g = c; b = x;
                } else if (hue < 240) {
                    r = 0; g = x; b = c;
                } else if (hue < 300) {
                    r = x; g = 0; b = c;
                } else {
                    r = c; g = 0; b = x;
                }
                
                palette[i] = Color(
                    static_cast<unsigned char>((r + m) * 255),
                    static_cast<unsigned char>((g + m) * 255),
                    static_cast<unsigned char>((b + m) * 255)
                );
            }
            break;
        }
        
        case InteriorVisualizationMethod::PERIOD_DETECTION: {
            // Distinct colors for different periods
            for (int i = 0; i < numColors; i++) {
                // Each period gets a distinct color
                double hue = (i * 137.5 + options.colorShift) / options.numLevels;
                hue = std::fmod(hue, 360.0);
                
                double saturation = 0.8;
                double value = 0.5;
                
                if (options.invertColors) {
                    hue = 360.0 - hue;
                }
                
                // Convert HSV to RGB
                double c = value * saturation;
                double x = c * (1 - std::abs(std::fmod(hue / 60.0, 2) - 1));
                double m = value - c;
                
                double r, g, b;
                if (hue < 60) {
                    r = c; g = x; b = 0;
                } else if (hue < 120) {
                    r = x; g = c; b = 0;
                } else if (hue < 180) {
                    r = 0; g = c; b = x;
                } else if (hue < 240) {
                    r = 0; g = x; b = c;
                } else if (hue < 300) {
                    r = x; g = 0; b = c;
                } else {
                    r = c; g = 0; b = x;
                }
                
                palette[i] = Color(
                    static_cast<unsigned char>((r + m) * 255),
                    static_cast<unsigned char>((g + m) * 255),
                    static_cast<unsigned char>((b + m) * 255)
                );
            }
            break;
        }
        
        case InteriorVisualizationMethod::BINARY_DECOMPOSITION: {
            // Colors based on binary patterns
            for (int i = 0; i < numColors; i++) {
                // Use bit patterns to create colors
                int r = 0, g = 0, b = 0;
                
                // Extract 3 bytes for RGB from decomposition value
                for (int bit = 0; bit < 8; bit++) {
                    if (i & (1 << bit)) {
                        if (bit < 3) r |= (1 << bit);
                        else if (bit < 6) g |= (1 << (bit - 3));
                        else b |= (1 << (bit - 6));
                    }
                }
                
                // Scale to reasonable range (not too bright, not too dark)
                r = 40 + r * 80 / 7;
                g = 40 + g * 80 / 7;
                b = 40 + b * 80 / 7;
                
                if (options.invertColors) {
                    r = 255 - r;
                    g = 255 - g;
                    b = 255 - b;
                }
                
                palette[i] = Color(
                    static_cast<unsigned char>(r),
                    static_cast<unsigned char>(g),
                    static_cast<unsigned char>(b)
                );
            }
            break;
        }
        
        case InteriorVisualizationMethod::POTENTIAL:
        case InteriorVisualizationMethod::DISTANCE_ESTIMATION:
        case InteriorVisualizationMethod::NORMALIZED_ITERATIONS:
        case InteriorVisualizationMethod::CURVATURE:
        case InteriorVisualizationMethod::INTERIOR_ORBIT_TRAPS:
        case InteriorVisualizationMethod::EXTERNAL_RAYS:
        default: {
            // Default to dark grayscale for unimplemented methods
            for (int i = 0; i < numColors; i++) {
                double position = static_cast<double>(i) / (numColors - 1);
                position = position * 0.3;  // Keep it dark
                
                if (options.invertColors) {
                    position = 0.3 - position;
                }
                
                unsigned char value = static_cast<unsigned char>(position * 255);
                palette[i] = Color(value, value, value);
            }
            break;
        }
    }
    
    return palette;
}

/**
 * @brief Visualize the interior (convergent) areas of a fractal set
 * 
 * @param width Image width
 * @param height Image height
 * @param xMin Minimum x-coordinate in complex plane
 * @param xMax Maximum x-coordinate in complex plane
 * @param yMin Minimum y-coordinate in complex plane
 * @param yMax Maximum y-coordinate in complex plane
 * @param maxIterations Maximum iterations
 * @param options Visualization options
 * @return std::vector<std::vector<Color>> Colored pixel data
 */
std::vector<std::vector<Color>> visualizeInterior(
    int width,
    int height,
    double xMin,
    double xMax,
    double yMin,
    double yMax,
    int maxIterations,
    const InteriorVisualizationOptions& options
) {
    std::vector<std::vector<Color>> result(height, std::vector<Color>(width));
    
    // Generate interior palette
    int numColors = std::max(16, options.numLevels);
    
    if (options.method == InteriorVisualizationMethod::BINARY_DECOMPOSITION) {
        numColors = 256;  // For binary decomposition, use more colors
    }
    
    std::vector<Color> interiorPalette = generateInteriorPalette(numColors, options);
    
    // Calculate step sizes
    double xStep = (xMax - xMin) / width;
    double yStep = (yMax - yMin) / height;
    
    // Process each pixel
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Map pixel to complex plane
            double real = xMin + x * xStep;
            double imag = yMax - y * yStep;
            Complex c(real, imag);
            
            // Default to black (for interior points)
            result[y][x] = Color(0, 0, 0);
            
            // Choose visualization method
            switch (options.method) {
                case InteriorVisualizationMethod::SOLID_BLACK: {
                    // Simple escape-time algorithm
                    Complex z(0.0, 0.0);
                    int i = 0;
                    
                    for (; i < maxIterations; i++) {
                        z = z * z + c;
                        
                        if (z.magnitudeSquared() > 4.0) {
                            break;
                        }
                    }
                    
                    if (i == maxIterations) {
                        // Interior point - black
                        result[y][x] = Color(0, 0, 0);
                    } else {
                        // Exterior point - white
                        result[y][x] = Color(255, 255, 255);
                    }
                    break;
                }
                
                case InteriorVisualizationMethod::DISTANCE_ESTIMATION: {
                    // Distance estimation
                    auto [distance, isInSet] = calculateDistanceEstimation(c, maxIterations, 2.0);
                    
                    if (isInSet) {
                        // Interior point - use level set coloring
                        double boundaryDistance = std::min(1.0, distance * 10.0);  // Scale for visibility
                        int level = static_cast<int>(boundaryDistance * options.numLevels);
                        level = std::min(options.numLevels - 1, std::max(0, level));
                        
                        result[y][x] = interiorPalette[level];
                    } else {
                        // Exterior point
                        double intensity = 1.0 - std::min(1.0, distance * 0.1);  // Scale for visibility
                        
                        if (options.showBoundary && intensity > 0.9) {
                            // Highlight boundary
                            result[y][x] = Color(255, 0, 0);  // Red for boundary
                        } else {
                            // Regular exterior
                            unsigned char value = static_cast<unsigned char>(intensity * 255);
                            result[y][x] = Color(value, value, value);
                        }
                    }
                    break;
                }
                
                case InteriorVisualizationMethod::LEVEL_SETS: {
                    // Interior level sets
                    Complex z(0.0, 0.0);
                    int i = 0;
                    
                    for (; i < maxIterations; i++) {
                        z = z * z + c;
                        
                        if (z.magnitudeSquared() > 4.0) {
                            break;
                        }
                    }
                    
                    if (i == maxIterations) {
                        // Interior point - use level sets
                        double zMag = z.magnitude();
                        int level = static_cast<int>(zMag * options.numLevels) % options.numLevels;
                        
                        result[y][x] = interiorPalette[level];
                    } else {
                        // Exterior point - white with boundary highlight
                        double smoothed = i + 1 - std::log2(std::log(z.magnitude()));
                        double normalized = smoothed / maxIterations;
                        
                        if (options.showBoundary && normalized < 0.05) {
                            // Highlight boundary
                            result[y][x] = Color(255, 0, 0);  // Red for boundary
                        } else {
                            // Regular exterior
                            result[y][x] = Color(255, 255, 255);
                        }
                    }
                    break;
                }
                
                case InteriorVisualizationMethod::PERIOD_DETECTION: {
                    // Period detection coloring
                    auto [period, detected] = detectPeriod(c, maxIterations, 2.0);
                    
                    if (detected && period > 0) {
                        // Interior point with detected period
                        int colorIndex = (period - 1) % numColors;
                        result[y][x] = interiorPalette[colorIndex];
                    } else {
                        // Either exterior or no period detected
                        Complex z(0.0, 0.0);
                        int i = 0;
                        
                        for (; i < maxIterations; i++) {
                            z = z * z + c;
                            
                            if (z.magnitudeSquared() > 4.0) {
                                break;
                            }
                        }
                        
                        if (i == maxIterations) {
                            // Interior point but no period detected - dark gray
                            result[y][x] = Color(30, 30, 30);
                        } else {
                            // Exterior point
                            result[y][x] = Color(255, 255, 255);
                        }
                    }
                    break;
                }
                
                case InteriorVisualizationMethod::BINARY_DECOMPOSITION: {
                    // Binary decomposition coloring
                    int decomposition = calculateBinaryDecomposition(c, maxIterations, 2.0);
                    
                    if (decomposition > 0) {
                        // Interior point with decomposition value
                        result[y][x] = interiorPalette[decomposition % numColors];
                    } else {
                        // Exterior point
                        Complex z(0.0, 0.0);
                        int i = 0;
                        
                        for (; i < maxIterations; i++) {
                            z = z * z + c;
                            
                            if (z.magnitudeSquared() > 4.0) {
                                break;
                            }
                        }
                        
                        if (i == maxIterations) {
                            // Interior point but no decomposition - dark gray
                            result[y][x] = Color(30, 30, 30);
                        } else {
                            // Exterior point
                            result[y][x] = Color(255, 255, 255);
                        }
                    }
                    break;
                }
                
                case InteriorVisualizationMethod::POTENTIAL:
                case InteriorVisualizationMethod::NORMALIZED_ITERATIONS:
                case InteriorVisualizationMethod::CURVATURE:
                case InteriorVisualizationMethod::INTERIOR_ORBIT_TRAPS:
                case InteriorVisualizationMethod::EXTERNAL_RAYS:
                default: {
                    // Unimplemented methods fall back to basic interior/exterior
                    Complex z(0.0, 0.0);
                    int i = 0;
                    
                    for (; i < maxIterations; i++) {
                        z = z * z + c;
                        
                        if (z.magnitudeSquared() > 4.0) {
                            break;
                        }
                    }
                    
                    if (i == maxIterations) {
                        // Interior point - black
                        result[y][x] = Color(0, 0, 0);
                    } else {
                        // Exterior point - white
                        result[y][x] = Color(255, 255, 255);
                    }
                    break;
                }
            }
        }
    }
    
    return result;
}

/**
 * @brief Save visualization to PPM file
 * 
 * @param filename Output filename
 * @param coloredData Colored pixel data
 * @return bool Success indicator
 */
bool saveVisualizationPPM(
    const std::string& filename,
    const std::vector<std::vector<Color>>& coloredData
) {
    try {
        int height = coloredData.size();
        int width = (height > 0) ? coloredData[0].size() : 0;
        
        if (height == 0 || width == 0) {
            return false;
        }
        
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            return false;
        }
        
        // Write PPM header
        file << "P6\n";
        file << width << " " << height << "\n";
        file << "255\n";
        
        // Write pixel data
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                const Color& color = coloredData[y][x];
                file.write(reinterpret_cast<const char*>(&color.r), 1);
                file.write(reinterpret_cast<const char*>(&color.g), 1);
                file.write(reinterpret_cast<const char*>(&color.b), 1);
            }
        }
        
        return true;
    } catch (...) {
        return false;
    }
}

/**
 * @brief Demonstrate interior visualization methods
 */
void demonstrateInteriorVisualizations() {
    // Settings for visualization
    const int WIDTH = 800;
    const int HEIGHT = 600;
    const int MAX_ITERATIONS = 1000;
    
    // Mandelbrot set view
    const double X_MIN = -2.0;
    const double X_MAX = 1.0;
    const double Y_MIN = -1.2;
    const double Y_MAX = 1.2;
    
    // Zoom into an interesting region with interior structure
    const double ZOOMED_X_MIN = -0.75 - 0.15;
    const double ZOOMED_X_MAX = -0.75 + 0.15;
    const double ZOOMED_Y_MIN = 0.0 - 0.15;
    const double ZOOMED_Y_MAX = 0.0 + 0.15;
    
    // Define visualization methods to demonstrate
    struct DemoMethod {
        InteriorVisualizationMethod method;
        std::string name;
        bool useZoomedView;
    };
    
    std::vector<DemoMethod> methods = {
        {InteriorVisualizationMethod::SOLID_BLACK, "Solid_Black", false},
        {InteriorVisualizationMethod::LEVEL_SETS, "Level_Sets", false},
        {InteriorVisualizationMethod::DISTANCE_ESTIMATION, "Distance_Estimation", false},
        {InteriorVisualizationMethod::PERIOD_DETECTION, "Period_Detection", true},
        {InteriorVisualizationMethod::BINARY_DECOMPOSITION, "Binary_Decomposition", true}
    };
    
    // Generate and save visualizations
    for (const auto& demo : methods) {
        std::cout << "Generating " << demo.name << " visualization..." << std::endl;
        
        // Set up options
        InteriorVisualizationOptions options(demo.method);
        
        // Choose view
        double xMin, xMax, yMin, yMax;
        if (demo.useZoomedView) {
            xMin = ZOOMED_X_MIN;
            xMax = ZOOMED_X_MAX;
            yMin = ZOOMED_Y_MIN;
            yMax = ZOOMED_Y_MAX;
        } else {
            xMin = X_MIN;
            xMax = X_MAX;
            yMin = Y_MIN;
            yMax = Y_MAX;
        }
        
        // Generate visualization
        auto visualization = visualizeInterior(
            WIDTH, HEIGHT, xMin, xMax, yMin, yMax, MAX_ITERATIONS, options
        );
        
        // Save the image
        std::string filename = "interior_" + demo.name + ".ppm";
        saveVisualizationPPM(filename, visualization);
        
        std::cout << "  Saved " << filename << std::endl;
    }
}

/**
 * @brief Demonstrate parameter variations for interior visualization
 */
void demonstrateParameterVariations() {
    // Settings for visualization
    const int WIDTH = 800;
    const int HEIGHT = 600;
    const int MAX_ITERATIONS = 1000;
    
    // Zoomed view for parameter testing
    const double X_MIN = -1.0 - 0.2;
    const double X_MAX = -1.0 + 0.2;
    const double Y_MIN = 0.0 - 0.2;
    const double Y_MAX = 0.0 + 0.2;
    
    // Define parameter variations
    struct ParamVariation {
        std::string name;
        int numLevels;
        double contrast;
        double brightness;
        bool invertColors;
        int colorShift;
    };
    
    std::vector<ParamVariation> variations = {
        {"Default", 8, 1.0, 1.0, false, 0},
        {"MoreLevels", 16, 1.0, 1.0, false, 0},
        {"HighContrast", 8, 2.0, 1.0, false, 0},
        {"Bright", 8, 1.0, 1.5, false, 0},
        {"Inverted", 8, 1.0, 1.0, true, 0},
        {"ColorShift", 8, 1.0, 1.0, false, 120}
    };
    
    // Use level sets for parameter testing
    InteriorVisualizationMethod method = InteriorVisualizationMethod::LEVEL_SETS;
    
    // Generate and save variations
    for (const auto& var : variations) {
        std::cout << "Generating parameter variation: " << var.name << "..." << std::endl;
        
        // Set up options
        InteriorVisualizationOptions options(
            method,
            var.numLevels,
            var.contrast,
            var.brightness,
            var.invertColors,
            var.colorShift
        );
        
        // Generate visualization
        auto visualization = visualizeInterior(
            WIDTH, HEIGHT, X_MIN, X_MAX, Y_MIN, Y_MAX, MAX_ITERATIONS, options
        );
        
        // Save the image
        std::string filename = "param_" + var.name + ".ppm";
        saveVisualizationPPM(filename, visualization);
        
        std::cout << "  Saved " << filename << std::endl;
    }
}

/**
 * @brief Main function to demonstrate interior visualization
 */
int main() {
    std::cout << "Demonstrating interior visualization methods..." << std::endl;
    demonstrateInteriorVisualizations();
    
    std::cout << std::endl << "Demonstrating parameter variations..." << std::endl;
    demonstrateParameterVariations();
    
    std::cout << std::endl << "All demonstrations complete." << std::endl;
    return 0;
}