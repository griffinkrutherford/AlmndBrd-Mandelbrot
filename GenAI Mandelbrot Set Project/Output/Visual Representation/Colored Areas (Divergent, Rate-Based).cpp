#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include <random>
#include "../../Code Structure (C++)/Visualization Class/Color Assignment Logic.cpp"
#include "../../Code Structure (C++)/Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @file Colored Areas (Divergent, Rate-Based).cpp
 * @brief Implementation of coloring techniques for divergent areas in fractal visualizations
 * 
 * This file provides methods for coloring the escape-time visualizations 
 * of fractal sets, focusing on the areas outside the set (divergent regions)
 * with color schemes based on the rate of divergence.
 */

/**
 * @brief Color mapping method for divergent areas
 */
enum class DivergentColorMethod {
    LINEAR,           // Linear mapping from iteration count to color
    LOGARITHMIC,      // Logarithmic mapping (reduces banding)
    SMOOTH,           // Smooth coloring based on continuous escape time
    HISTOGRAM,        // Histogram equalization (better distribution)
    ANGLE_BASED,      // Based on angle of last z value
    TRIANGLE_INEQUALITY,  // Using triangle inequality average
    ORBIT_TRAP,       // Based on distances to geometric shapes
    CURVATURE,        // Based on trajectory curvature
    BIOMORPHIC        // Biomorphic-inspired coloring
};

/**
 * @brief Color palette style for divergent coloring
 */
enum class PaletteStyle {
    RAINBOW,          // Full rainbow spectrum
    GRADIENT,         // Two-color gradient
    MULTI_GRADIENT,   // Multi-color gradient
    BANDED,           // Distinct color bands
    FIRE,             // Fire-like colors (black to yellow through red)
    ELECTRIC,         // Electric-like colors (purples, blues, cyans)
    TERRAIN,          // Terrain-inspired colors (greens, browns)
    OCEAN,            // Ocean-inspired colors (blues, teals)
    GRAYSCALE,        // Black to white
    CUSTOM            // User-defined palette
};

/**
 * @brief Structure for customizing color mapping
 */
struct ColorMappingOptions {
    DivergentColorMethod method;     // Coloring method
    PaletteStyle style;              // Color palette style
    double colorCycles;              // Number of color cycles (for cyclic palettes)
    bool invertPalette;              // Whether to invert the color palette
    double brightness;               // Overall brightness adjustment (0.0-2.0)
    double contrast;                 // Overall contrast adjustment (0.0-2.0)
    double saturation;               // Overall saturation adjustment (0.0-2.0)
    int paletteShift;                // Shift the palette start point
    
    // Constructor with defaults
    ColorMappingOptions(
        DivergentColorMethod m = DivergentColorMethod::SMOOTH,
        PaletteStyle s = PaletteStyle::RAINBOW,
        double c = 1.0,
        bool i = false,
        double b = 1.0,
        double con = 1.0,
        double sat = 1.0,
        int shift = 0
    ) : method(m),
        style(s),
        colorCycles(c),
        invertPalette(i),
        brightness(b),
        contrast(con),
        saturation(sat),
        paletteShift(shift) {
    }
};

/**
 * @brief Generate a custom color palette
 * 
 * @param paletteStyle Style of palette to generate
 * @param colorCount Number of colors in the palette
 * @param options Color mapping options
 * @return std::vector<Color> Generated color palette
 */
std::vector<Color> generateCustomPalette(
    PaletteStyle paletteStyle,
    int colorCount,
    const ColorMappingOptions& options
) {
    std::vector<Color> palette(colorCount);
    
    switch (paletteStyle) {
        case PaletteStyle::RAINBOW: {
            // Rainbow palette (HSV color wheel)
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / colorCount;
                
                // Apply color cycles
                double hue = position * 360.0 * options.colorCycles;
                hue = std::fmod(hue + options.paletteShift, 360.0);
                
                // Invert if requested
                if (options.invertPalette) {
                    hue = 360.0 - hue;
                }
                
                // Convert HSV to RGB
                double saturation = 1.0 * options.saturation;
                double value = 1.0 * options.brightness;
                
                // Clamp values
                saturation = std::min(1.0, std::max(0.0, saturation));
                value = std::min(1.0, std::max(0.0, value));
                
                // HSV to RGB conversion
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
                
                // Apply contrast around 0.5
                auto applyContrast = [&options](double value) {
                    double normalized = value - 0.5;
                    normalized *= options.contrast;
                    return std::min(1.0, std::max(0.0, normalized + 0.5));
                };
                
                r = applyContrast(r + m);
                g = applyContrast(g + m);
                b = applyContrast(b + m);
                
                palette[i] = Color(
                    static_cast<unsigned char>(r * 255),
                    static_cast<unsigned char>(g * 255),
                    static_cast<unsigned char>(b * 255)
                );
            }
            break;
        }
            
        case PaletteStyle::GRADIENT: {
            // Two-color gradient
            Color startColor(0, 0, 0);      // Black
            Color endColor(255, 255, 255);  // White
            
            // Can be customized with specific colors
            switch (options.paletteShift % 6) {
                case 0: // Black to white
                    startColor = Color(0, 0, 0);
                    endColor = Color(255, 255, 255);
                    break;
                case 1: // Blue to yellow
                    startColor = Color(0, 0, 255);
                    endColor = Color(255, 255, 0);
                    break;
                case 2: // Purple to green
                    startColor = Color(128, 0, 128);
                    endColor = Color(0, 255, 0);
                    break;
                case 3: // Red to cyan
                    startColor = Color(255, 0, 0);
                    endColor = Color(0, 255, 255);
                    break;
                case 4: // Dark blue to light blue
                    startColor = Color(0, 0, 64);
                    endColor = Color(0, 128, 255);
                    break;
                case 5: // Dark red to yellow
                    startColor = Color(64, 0, 0);
                    endColor = Color(255, 255, 0);
                    break;
            }
            
            if (options.invertPalette) {
                std::swap(startColor, endColor);
            }
            
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / (colorCount - 1);
                
                // Apply contrast
                position = 0.5 + (position - 0.5) * options.contrast;
                position = std::min(1.0, std::max(0.0, position));
                
                // Linear interpolation
                unsigned char r = static_cast<unsigned char>(
                    startColor.r + position * (endColor.r - startColor.r)
                );
                unsigned char g = static_cast<unsigned char>(
                    startColor.g + position * (endColor.g - startColor.g)
                );
                unsigned char b = static_cast<unsigned char>(
                    startColor.b + position * (endColor.b - startColor.b)
                );
                
                // Apply brightness
                auto adjustBrightness = [&options](unsigned char value) {
                    double normalized = value / 255.0 * options.brightness;
                    return static_cast<unsigned char>(std::min(255.0, std::max(0.0, normalized * 255)));
                };
                
                r = adjustBrightness(r);
                g = adjustBrightness(g);
                b = adjustBrightness(b);
                
                palette[i] = Color(r, g, b);
            }
            break;
        }
            
        case PaletteStyle::FIRE: {
            // Fire-like colors (black -> red -> orange -> yellow)
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / (colorCount - 1);
                
                // Apply color cycles and shifts
                position = std::fmod(position * options.colorCycles + options.paletteShift / 360.0, 1.0);
                
                // Invert if requested
                if (options.invertPalette) {
                    position = 1.0 - position;
                }
                
                // Apply contrast
                position = 0.5 + (position - 0.5) * options.contrast;
                position = std::min(1.0, std::max(0.0, position));
                
                // Fire palette formula
                double r = position < 0.5 ? position * 2 : 1.0;
                double g = position < 0.25 ? 0 : position < 0.75 ? (position - 0.25) * 2 : 1.0;
                double b = position < 0.5 ? 0 : (position - 0.5) * 2;
                
                // Apply brightness and saturation
                auto adjustColor = [&options](double value) {
                    // Brightness
                    value *= options.brightness;
                    
                    // Saturation (pull toward middle gray)
                    if (options.saturation < 1.0) {
                        value = 0.5 + (value - 0.5) * options.saturation;
                    }
                    
                    return static_cast<unsigned char>(std::min(255.0, std::max(0.0, value * 255)));
                };
                
                palette[i] = Color(
                    adjustColor(r),
                    adjustColor(g),
                    adjustColor(b)
                );
            }
            break;
        }
            
        case PaletteStyle::OCEAN: {
            // Ocean colors (dark blue -> teal -> light blue)
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / (colorCount - 1);
                
                // Apply color cycles and shifts
                position = std::fmod(position * options.colorCycles + options.paletteShift / 360.0, 1.0);
                
                // Invert if requested
                if (options.invertPalette) {
                    position = 1.0 - position;
                }
                
                // Apply contrast
                position = 0.5 + (position - 0.5) * options.contrast;
                position = std::min(1.0, std::max(0.0, position));
                
                // Ocean palette formula
                double r = position < 0.25 ? 0 : position < 0.5 ? (position - 0.25) * 4 : position < 0.75 ? 1.0 : 1.0;
                double g = position < 0.25 ? 0 : position < 0.5 ? (position - 0.25) * 4 : position < 0.75 ? 1.0 : (1.0 - position) * 4;
                double b = position < 0.25 ? position * 4 : 1.0;
                
                // Apply brightness and saturation
                auto adjustColor = [&options](double value) {
                    // Brightness
                    value *= options.brightness;
                    
                    // Saturation (pull toward middle gray)
                    if (options.saturation < 1.0) {
                        value = 0.5 + (value - 0.5) * options.saturation;
                    }
                    
                    return static_cast<unsigned char>(std::min(255.0, std::max(0.0, value * 255)));
                };
                
                palette[i] = Color(
                    adjustColor(r),
                    adjustColor(g),
                    adjustColor(b)
                );
            }
            break;
        }
            
        case PaletteStyle::GRAYSCALE: {
            // Grayscale (black to white)
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / (colorCount - 1);
                
                // Apply color cycles
                position = std::fmod(position * options.colorCycles + options.paletteShift / 360.0, 1.0);
                
                // Invert if requested
                if (options.invertPalette) {
                    position = 1.0 - position;
                }
                
                // Apply contrast
                position = 0.5 + (position - 0.5) * options.contrast;
                position = std::min(1.0, std::max(0.0, position));
                
                // Apply brightness
                position *= options.brightness;
                position = std::min(1.0, std::max(0.0, position));
                
                unsigned char value = static_cast<unsigned char>(position * 255);
                palette[i] = Color(value, value, value);
            }
            break;
        }
            
        case PaletteStyle::MULTI_GRADIENT:
        case PaletteStyle::BANDED:
        case PaletteStyle::ELECTRIC:
        case PaletteStyle::TERRAIN:
        case PaletteStyle::CUSTOM:
        default: {
            // Default to rainbow if not implemented
            // (This would be replaced with actual implementations)
            for (int i = 0; i < colorCount; i++) {
                double position = static_cast<double>(i) / colorCount;
                double hue = position * 360.0;
                
                // Simple HSV to RGB conversion
                double saturation = 1.0;
                double value = 1.0;
                
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
    }
    
    return palette;
}

/**
 * @brief Calculate continuous escape time for smooth coloring
 * 
 * @param iterationCount Discrete iteration count
 * @param maxIterations Maximum iterations
 * @param lastZ Final z value
 * @return double Continuous escape time
 */
double calculateContinuousEscapeTime(
    int iterationCount,
    int maxIterations,
    const Complex& lastZ
) {
    if (iterationCount == maxIterations) {
        return static_cast<double>(maxIterations);
    }
    
    // Calculate normalized iteration count
    double zn = lastZ.magnitudeSquared();
    double nu = std::log(std::log(zn) / std::log(2.0)) / std::log(2.0);
    return static_cast<double>(iterationCount) + 1.0 - nu;
}

/**
 * @brief Calculate angle-based coloring
 * 
 * @param lastZ Final z value
 * @return double Angle-based color parameter (0.0-1.0)
 */
double calculateAngleColoring(const Complex& lastZ) {
    // Calculate the angle of the last z value
    double angle = std::atan2(lastZ.getImag(), lastZ.getReal());
    
    // Normalize to 0.0-1.0 range
    return (angle + M_PI) / (2.0 * M_PI);
}

/**
 * @brief Calculate triangle inequality average coloring
 * 
 * @param orbitHistory History of z values during iteration
 * @return double Triangle inequality parameter (0.0-1.0)
 */
double calculateTriangleInequalityAverage(const std::vector<Complex>& orbitHistory) {
    if (orbitHistory.size() < 3) {
        return 0.0;
    }
    
    double sum = 0.0;
    int count = 0;
    
    for (size_t i = 2; i < orbitHistory.size(); i++) {
        // Calculate sides of triangle formed by origin, current, and previous point
        double a = orbitHistory[i].magnitude();
        double b = orbitHistory[i-1].magnitude();
        double c = (orbitHistory[i] + orbitHistory[i-1] * Complex(-1, 0)).magnitude();
        
        // Skip if any side is too small
        if (a < 1e-10 || b < 1e-10 || c < 1e-10) {
            continue;
        }
        
        // Calculate triangle inequality average
        double value = (a + b + c) / (3.0 * std::max({a, b, c}));
        sum += value;
        count++;
    }
    
    return (count > 0) ? (sum / count) : 0.0;
}

/**
 * @brief Calculate orbit trap coloring
 * 
 * @param orbitHistory History of z values during iteration
 * @param trapType Type of trap (0=point, 1=circle, 2=line, 3=cross)
 * @return double Orbit trap parameter (0.0-1.0)
 */
double calculateOrbitTrap(
    const std::vector<Complex>& orbitHistory,
    int trapType = 0
) {
    if (orbitHistory.empty()) {
        return 0.0;
    }
    
    double minDistance = std::numeric_limits<double>::max();
    
    for (const auto& z : orbitHistory) {
        double distance;
        
        switch (trapType) {
            case 0: {
                // Point trap at origin
                distance = z.magnitude();
                break;
            }
            case 1: {
                // Circle trap with radius 1
                double mag = z.magnitude();
                distance = std::abs(mag - 1.0);
                break;
            }
            case 2: {
                // Line trap along real axis
                distance = std::abs(z.getImag());
                break;
            }
            case 3: {
                // Cross trap (real and imaginary axes)
                distance = std::min(std::abs(z.getReal()), std::abs(z.getImag()));
                break;
            }
            default:
                distance = z.magnitude();
        }
        
        minDistance = std::min(minDistance, distance);
    }
    
    // Normalize to 0.0-1.0 (using a reasonable maximum distance)
    return std::min(1.0, minDistance / 2.0);
}

/**
 * @brief Map iteration data to colors using specified coloring method
 * 
 * @param iterationData 2D array of iteration counts
 * @param maxIterations Maximum iterations
 * @param zValues 2D array of final z values (optional, used for some methods)
 * @param orbitHistories Vector of orbit histories (optional, used for some methods)
 * @param options Color mapping options
 * @return std::vector<std::vector<Color>> Colored pixel data
 */
std::vector<std::vector<Color>> colorDivergentAreas(
    const std::vector<std::vector<int>>& iterationData,
    int maxIterations,
    const std::vector<std::vector<Complex>>& zValues = {},
    const std::vector<std::vector<std::vector<Complex>>>& orbitHistories = {},
    const ColorMappingOptions& options = ColorMappingOptions()
) {
    // Get dimensions
    int height = iterationData.size();
    int width = (height > 0) ? iterationData[0].size() : 0;
    
    // Create output array
    std::vector<std::vector<Color>> coloredData(height, std::vector<Color>(width));
    
    // Create palette (with extra colors for interpolation)
    int paletteSize = maxIterations + 1;
    std::vector<Color> palette = generateCustomPalette(options.style, paletteSize, options);
    
    // Helper function for histogram equalization
    auto calculateHistogram = [&]() -> std::vector<double> {
        // Count iterations
        std::vector<int> histogram(maxIterations + 1, 0);
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int iterations = iterationData[y][x];
                if (iterations < maxIterations) {  // Only count escaping points
                    histogram[iterations]++;
                }
            }
        }
        
        // Calculate cumulative distribution
        std::vector<double> cdf(maxIterations + 1, 0.0);
        int total = 0;
        
        for (int i = 0; i <= maxIterations; i++) {
            total += histogram[i];
        }
        
        if (total > 0) {
            int sum = 0;
            for (int i = 0; i <= maxIterations; i++) {
                sum += histogram[i];
                cdf[i] = static_cast<double>(sum) / total;
            }
        }
        
        return cdf;
    };
    
    // Prepare data for advanced coloring methods
    bool hasZValues = !zValues.empty() && zValues.size() == height && 
                      !zValues.empty() && zValues[0].size() == width;
                      
    bool hasOrbitHistories = !orbitHistories.empty() && orbitHistories.size() == height && 
                             !orbitHistories.empty() && orbitHistories[0].size() == width;
    
    // Calculate histogram if needed
    std::vector<double> cdf;
    if (options.method == DivergentColorMethod::HISTOGRAM) {
        cdf = calculateHistogram();
    }
    
    // Map iteration counts to colors
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iterations = iterationData[y][x];
            
            // Points in the set (didn't escape) are colored black
            if (iterations == maxIterations) {
                coloredData[y][x] = Color(0, 0, 0);
                continue;
            }
            
            // Calculate color index based on coloring method
            double colorIndex = 0.0;
            
            switch (options.method) {
                case DivergentColorMethod::LINEAR: {
                    // Simple linear mapping
                    colorIndex = static_cast<double>(iterations) / maxIterations;
                    break;
                }
                
                case DivergentColorMethod::LOGARITHMIC: {
                    // Logarithmic mapping to reduce banding
                    colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    break;
                }
                
                case DivergentColorMethod::SMOOTH: {
                    // Smooth coloring requires z values
                    if (hasZValues) {
                        double smooth = calculateContinuousEscapeTime(
                            iterations, maxIterations, zValues[y][x]
                        );
                        colorIndex = smooth / maxIterations;
                    } else {
                        // Fall back to logarithmic if z values not available
                        colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    }
                    break;
                }
                
                case DivergentColorMethod::HISTOGRAM: {
                    // Use cumulative distribution
                    colorIndex = cdf[iterations];
                    break;
                }
                
                case DivergentColorMethod::ANGLE_BASED: {
                    // Angle-based coloring requires z values
                    if (hasZValues) {
                        double angle = calculateAngleColoring(zValues[y][x]);
                        // Mix iteration count with angle
                        double iterPart = static_cast<double>(iterations) / maxIterations;
                        colorIndex = 0.7 * iterPart + 0.3 * angle;
                    } else {
                        // Fall back to logarithmic if z values not available
                        colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    }
                    break;
                }
                
                case DivergentColorMethod::TRIANGLE_INEQUALITY: {
                    // Triangle inequality average requires orbit histories
                    if (hasOrbitHistories) {
                        double tia = calculateTriangleInequalityAverage(orbitHistories[y][x]);
                        // Mix iteration count with triangle inequality
                        double iterPart = static_cast<double>(iterations) / maxIterations;
                        colorIndex = 0.6 * iterPart + 0.4 * tia;
                    } else {
                        // Fall back to logarithmic if orbit histories not available
                        colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    }
                    break;
                }
                
                case DivergentColorMethod::ORBIT_TRAP: {
                    // Orbit trap coloring requires orbit histories
                    if (hasOrbitHistories) {
                        double trap = calculateOrbitTrap(orbitHistories[y][x], options.paletteShift % 4);
                        // Mix iteration count with orbit trap
                        double iterPart = static_cast<double>(iterations) / maxIterations;
                        colorIndex = 0.5 * iterPart + 0.5 * trap;
                    } else {
                        // Fall back to logarithmic if orbit histories not available
                        colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    }
                    break;
                }
                
                case DivergentColorMethod::CURVATURE: 
                case DivergentColorMethod::BIOMORPHIC:
                default: {
                    // Fall back to logarithmic for unimplemented methods
                    colorIndex = std::log(1.0 + iterations) / std::log(1.0 + maxIterations);
                    break;
                }
            }
            
            // Apply color cycles
            colorIndex = std::fmod(colorIndex * options.colorCycles, 1.0);
            
            // Map to palette
            int index = static_cast<int>(colorIndex * (paletteSize - 1));
            index = std::min(paletteSize - 1, std::max(0, index));
            
            coloredData[y][x] = palette[index];
        }
    }
    
    return coloredData;
}

/**
 * @brief Create test pattern for demonstrating coloring methods
 * 
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iterations
 * @return std::pair<std::vector<std::vector<int>>, std::vector<std::vector<Complex>>> 
 *         Pair of iteration counts and final z values
 */
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<Complex>>>
createTestPattern(int width, int height, int maxIterations) {
    std::vector<std::vector<int>> iterations(height, std::vector<int>(width));
    std::vector<std::vector<Complex>> zValues(height, std::vector<Complex>(width));
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Create a simple pattern
            double centerX = width / 2.0;
            double centerY = height / 2.0;
            double dx = x - centerX;
            double dy = y - centerY;
            double distance = std::sqrt(dx*dx + dy*dy);
            double angle = std::atan2(dy, dx);
            
            // Make iterations vary with distance from center
            int iteration = static_cast<int>(distance / (std::max(width, height) / 2.0) * maxIterations);
            iteration = std::min(maxIterations - 1, iteration);  // Keep all points "escaping"
            
            // Make z values vary with angle
            Complex z(std::cos(angle) * 2.0, std::sin(angle) * 2.0);
            
            iterations[y][x] = iteration;
            zValues[y][x] = z;
        }
    }
    
    return {iterations, zValues};
}

/**
 * @brief Save colored image to PPM file
 * 
 * @param filename Output filename
 * @param coloredData Colored pixel data
 * @return bool Success indicator
 */
bool saveColoredImagePPM(
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
 * @brief Demonstrate different coloring methods
 * 
 * @param width Image width
 * @param height Image height
 */
void demonstrateColoringMethods(int width, int height) {
    const int MAX_ITERATIONS = 100;
    
    // Create test pattern
    auto [iterations, zValues] = createTestPattern(width, height, MAX_ITERATIONS);
    
    // Create artificial orbit histories (simplified for demonstration)
    std::vector<std::vector<std::vector<Complex>>> orbitHistories(
        height, std::vector<std::vector<Complex>>(width)
    );
    
    // Fill with simple orbits
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iter = iterations[y][x];
            
            // Create a simple orbit history
            std::vector<Complex> orbit;
            orbit.push_back(Complex(0, 0));  // Start at origin
            
            // Create artificial orbit points
            double angle = std::atan2(y - height/2.0, x - width/2.0);
            for (int i = 1; i <= iter; i++) {
                double scale = static_cast<double>(i) / MAX_ITERATIONS * 2.0;
                orbit.push_back(Complex(
                    std::cos(angle) * scale,
                    std::sin(angle) * scale
                ));
            }
            
            orbitHistories[y][x] = orbit;
        }
    }
    
    // Define coloring methods to demonstrate
    struct DemoMethod {
        DivergentColorMethod method;
        PaletteStyle palette;
        std::string name;
    };
    
    std::vector<DemoMethod> methods = {
        {DivergentColorMethod::LINEAR, PaletteStyle::RAINBOW, "Linear_Rainbow"},
        {DivergentColorMethod::LOGARITHMIC, PaletteStyle::RAINBOW, "Log_Rainbow"},
        {DivergentColorMethod::SMOOTH, PaletteStyle::RAINBOW, "Smooth_Rainbow"},
        {DivergentColorMethod::LINEAR, PaletteStyle::GRADIENT, "Linear_Gradient"},
        {DivergentColorMethod::LINEAR, PaletteStyle::FIRE, "Linear_Fire"},
        {DivergentColorMethod::LINEAR, PaletteStyle::OCEAN, "Linear_Ocean"},
        {DivergentColorMethod::LOGARITHMIC, PaletteStyle::GRADIENT, "Log_Gradient"},
        {DivergentColorMethod::ANGLE_BASED, PaletteStyle::RAINBOW, "Angle_Rainbow"},
        {DivergentColorMethod::TRIANGLE_INEQUALITY, PaletteStyle::RAINBOW, "Triangle_Rainbow"},
        {DivergentColorMethod::ORBIT_TRAP, PaletteStyle::RAINBOW, "OrbitTrap_Rainbow"}
    };
    
    // Generate and save images for each method
    for (const auto& demo : methods) {
        // Create options
        ColorMappingOptions options(demo.method, demo.palette);
        
        // Color the pattern
        auto coloredData = colorDivergentAreas(
            iterations, MAX_ITERATIONS, zValues, orbitHistories, options
        );
        
        // Save the image
        std::string filename = "colored_" + demo.name + ".ppm";
        saveColoredImagePPM(filename, coloredData);
        
        std::cout << "Generated " << filename << std::endl;
    }
}

/**
 * @brief Demonstrate color parameter effects
 * 
 * @param width Image width
 * @param height Image height
 */
void demonstrateColorParameters(int width, int height) {
    const int MAX_ITERATIONS = 100;
    
    // Create test pattern
    auto [iterations, zValues] = createTestPattern(width, height, MAX_ITERATIONS);
    
    // Define parameter variations
    struct ParamDemo {
        std::string name;
        double colorCycles;
        bool invertPalette;
        double brightness;
        double contrast;
        double saturation;
        int paletteShift;
    };
    
    std::vector<ParamDemo> demos = {
        {"Default", 1.0, false, 1.0, 1.0, 1.0, 0},
        {"DoubleCycle", 2.0, false, 1.0, 1.0, 1.0, 0},
        {"Inverted", 1.0, true, 1.0, 1.0, 1.0, 0},
        {"HighContrast", 1.0, false, 1.0, 2.0, 1.0, 0},
        {"LowContrast", 1.0, false, 1.0, 0.5, 1.0, 0},
        {"Bright", 1.0, false, 1.5, 1.0, 1.0, 0},
        {"Dark", 1.0, false, 0.5, 1.0, 1.0, 0},
        {"Saturated", 1.0, false, 1.0, 1.0, 1.5, 0},
        {"Desaturated", 1.0, false, 1.0, 1.0, 0.5, 0},
        {"Shifted", 1.0, false, 1.0, 1.0, 1.0, 120}
    };
    
    // Generate and save images for each parameter set
    for (const auto& demo : demos) {
        // Create options
        ColorMappingOptions options(
            DivergentColorMethod::SMOOTH,
            PaletteStyle::RAINBOW,
            demo.colorCycles,
            demo.invertPalette,
            demo.brightness,
            demo.contrast,
            demo.saturation,
            demo.paletteShift
        );
        
        // Color the pattern
        auto coloredData = colorDivergentAreas(
            iterations, MAX_ITERATIONS, zValues, {}, options
        );
        
        // Save the image
        std::string filename = "param_" + demo.name + ".ppm";
        saveColoredImagePPM(filename, coloredData);
        
        std::cout << "Generated " << filename << std::endl;
    }
}

/**
 * @brief Main function to demonstrate divergent area coloring
 */
int main() {
    // Demonstrate different coloring methods
    std::cout << "Demonstrating coloring methods..." << std::endl;
    demonstrateColoringMethods(600, 400);
    
    // Demonstrate color parameter effects
    std::cout << std::endl << "Demonstrating color parameters..." << std::endl;
    demonstrateColorParameters(600, 400);
    
    std::cout << std::endl << "All demonstrations complete." << std::endl;
    return 0;
}