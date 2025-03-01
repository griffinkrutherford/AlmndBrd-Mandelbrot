#include <vector>
#include <algorithm>
#include <cmath>
#include "../../Base Class/Properties/Iteration Parameters.h"

/**
 * @brief RGB Color structure
 */
struct Color {
    unsigned char r, g, b;
    
    Color(unsigned char red = 0, unsigned char green = 0, unsigned char blue = 0)
        : r(red), g(green), b(blue) {}
};

/**
 * @brief Color mapping modes for Mandelbrot visualization
 */
enum class ColorMode {
    LINEAR,         // Linear mapping from iterations to color
    LOGARITHMIC,    // Logarithmic mapping for smoother transitions
    SMOOTH,         // Smooth coloring with fractional iterations
    HISTOGRAM,      // Histogram equalization for better distribution
    CYCLE           // Cyclic color palette
};

/**
 * @brief Maps iteration count to RGB color using linear gradient
 * 
 * @param iteration Iteration count
 * @param maxIterations Maximum iterations
 * @return Color RGB color
 */
Color linearColorMap(int iteration, int maxIterations) {
    if (iteration == maxIterations) {
        // Point is in the set (didn't escape)
        return Color(0, 0, 0);
    }
    
    // Simple linear mapping from iteration to color
    double ratio = static_cast<double>(iteration) / maxIterations;
    
    // Create a blue-to-white gradient
    unsigned char value = static_cast<unsigned char>(255 * ratio);
    return Color(value, value, 255);
}

/**
 * @brief Maps iteration count to RGB color using logarithmic scale
 * 
 * @param iteration Iteration count
 * @param maxIterations Maximum iterations
 * @return Color RGB color
 */
Color logarithmicColorMap(int iteration, int maxIterations) {
    if (iteration == maxIterations) {
        // Point is in the set (didn't escape)
        return Color(0, 0, 0);
    }
    
    // Logarithmic mapping for smoother color transitions
    double ratio = log(static_cast<double>(iteration)) / log(static_cast<double>(maxIterations));
    ratio = std::min(std::max(ratio, 0.0), 1.0); // Clamp to [0, 1]
    
    // Create a color based on position in spectrum
    double hue = 360.0 * ratio;
    double saturation = 1.0;
    double value = 1.0;
    
    // Convert HSV to RGB
    double c = value * saturation;
    double x = c * (1 - std::abs(fmod(hue / 60.0, 2) - 1));
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
    
    return Color(
        static_cast<unsigned char>((r + m) * 255),
        static_cast<unsigned char>((g + m) * 255),
        static_cast<unsigned char>((b + m) * 255)
    );
}

/**
 * @brief Maps iteration count to RGB color using smooth coloring technique
 * 
 * @param iteration Iteration count
 * @param lastZ Final z value when iteration ended
 * @param maxIterations Maximum iterations
 * @return Color RGB color
 */
Color smoothColorMap(int iteration, const Complex& lastZ, int maxIterations) {
    if (iteration == maxIterations) {
        // Point is in the set (didn't escape)
        return Color(0, 0, 0);
    }
    
    // Smooth coloring formula based on the continuous iteration count
    double zn = lastZ.magnitudeSquared();
    double nu = log(log(zn) / 2.0) / log(2.0);
    double smoothed = iteration + 1 - nu;
    
    // Normalize
    double ratio = smoothed / maxIterations;
    ratio = std::min(std::max(ratio, 0.0), 1.0); // Clamp to [0, 1]
    
    // Create colors based on smooth iteration value
    // This creates a rainbow-like pattern
    double r = 9 * (1 - ratio) * ratio * ratio * ratio * 255;
    double g = 15 * (1 - ratio) * (1 - ratio) * ratio * ratio * 255;
    double b = 8.5 * (1 - ratio) * (1 - ratio) * (1 - ratio) * ratio * 255;
    
    return Color(
        static_cast<unsigned char>(r),
        static_cast<unsigned char>(g),
        static_cast<unsigned char>(b)
    );
}

/**
 * @brief Create an array of colors for use in visualization
 * 
 * @param colorMode Color mapping mode
 * @param maxIterations Maximum iterations
 * @return std::vector<Color> Array of colors
 */
std::vector<Color> generateColorPalette(ColorMode colorMode, int maxIterations) {
    std::vector<Color> palette;
    palette.reserve(maxIterations + 1);
    
    Complex dummy(0, 0); // Dummy complex for smooth coloring
    
    for (int i = 0; i <= maxIterations; i++) {
        switch (colorMode) {
            case ColorMode::LINEAR:
                palette.push_back(linearColorMap(i, maxIterations));
                break;
            case ColorMode::LOGARITHMIC:
                palette.push_back(logarithmicColorMap(i, maxIterations));
                break;
            case ColorMode::SMOOTH:
                palette.push_back(smoothColorMap(i, dummy, maxIterations));
                break;
            case ColorMode::CYCLE: {
                // Cycle through colors
                double ratio = static_cast<double>(i % 16) / 16.0;
                double r = sin(ratio * 2 * M_PI) * 127 + 128;
                double g = sin((ratio * 2 * M_PI) + 2 * M_PI / 3) * 127 + 128;
                double b = sin((ratio * 2 * M_PI) + 4 * M_PI / 3) * 127 + 128;
                palette.push_back(Color(r, g, b));
                break;
            }
            default:
                palette.push_back(linearColorMap(i, maxIterations));
                break;
        }
    }
    
    return palette;
}