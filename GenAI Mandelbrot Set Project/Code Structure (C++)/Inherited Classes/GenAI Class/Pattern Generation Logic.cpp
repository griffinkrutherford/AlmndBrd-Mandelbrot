#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include <string>
#include <functional>
#include "../../Base Class/Properties/Complex Number Handling (z, c).h"

/**
 * @brief Pattern generation modes for fractal exploration
 */
enum class PatternMode {
    MANDELBROT,          // Standard Mandelbrot set
    JULIA,               // Julia set
    BURNING_SHIP,        // Burning Ship fractal
    MULTIBROT,           // Generalized Mandelbrot with higher powers
    HYBRID,              // Hybrid of multiple fractals
    RANDOMIZED           // Random modifications to fractal equations
};

/**
 * @brief Generate a modified iteration function for fractal exploration
 * 
 * @param mode Pattern generation mode
 * @param seed Random seed for reproducibility
 * @return std::function<Complex(Complex, Complex)> Modified iteration function
 */
std::function<Complex(Complex, Complex)> 
generateIterationFunction(PatternMode mode, unsigned int seed = 42) {
    // Default Mandelbrot iteration function: z = z² + c
    auto mandelbrotFn = [](const Complex& z, const Complex& c) -> Complex {
        return z * z + c;
    };
    
    // Initialize random number generator
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    
    switch (mode) {
        case PatternMode::MANDELBROT:
            return mandelbrotFn;
            
        case PatternMode::JULIA:
            return mandelbrotFn; // Same formula, different usage
            
        case PatternMode::BURNING_SHIP: {
            // Burning Ship uses absolute values of real and imaginary parts
            return [](const Complex& z, const Complex& c) -> Complex {
                double absReal = std::abs(z.getReal());
                double absImag = std::abs(z.getImag());
                Complex absZ(absReal, absImag);
                return absZ * absZ + c;
            };
        }
            
        case PatternMode::MULTIBROT: {
            // Generate a random power between 3 and 6
            std::uniform_int_distribution<int> powerDist(3, 6);
            int power = powerDist(gen);
            
            return [power](const Complex& z, const Complex& c) -> Complex {
                Complex result = z;
                for (int i = 1; i < power; i++) {
                    result = result * z;
                }
                return result + c;
            };
        }
            
        case PatternMode::HYBRID: {
            // Randomly combine different fractal formulas
            std::uniform_int_distribution<int> typeDist(0, 2);
            int type = typeDist(gen);
            
            if (type == 0) {
                // Combination of Mandelbrot and Burning Ship
                return [](const Complex& z, const Complex& c) -> Complex {
                    double absReal = std::abs(z.getReal());
                    double imag = z.getImag();
                    Complex modZ(absReal, imag);
                    return modZ * modZ + c;
                };
            } 
            else if (type == 1) {
                // Sine transformation
                return [](const Complex& z, const Complex& c) -> Complex {
                    double sinReal = std::sin(z.getReal());
                    double sinImag = std::sin(z.getImag());
                    Complex sinZ(sinReal, sinImag);
                    return sinZ * z + c;
                };
            }
            else {
                // Cosine transformation
                return [](const Complex& z, const Complex& c) -> Complex {
                    double cosReal = std::cos(z.getReal());
                    double cosImag = std::cos(z.getImag());
                    Complex cosZ(cosReal, cosImag);
                    return z * z + cosZ * c;
                };
            }
        }
            
        case PatternMode::RANDOMIZED: {
            // Generate random coefficients
            double a = dist(gen);
            double b = dist(gen);
            double c = dist(gen);
            double d = dist(gen);
            
            std::uniform_int_distribution<int> opDist(0, 3);
            int operation = opDist(gen);
            
            return [a, b, c, d, operation](const Complex& z, const Complex& c) -> Complex {
                Complex z2 = z * z;
                Complex scaled = Complex(z.getReal() * a, z.getImag() * b);
                Complex mod = Complex(c.getReal() * c, c.getImag() * d);
                
                switch (operation) {
                    case 0:
                        return z2 + scaled + mod;
                    case 1:
                        return z2 * scaled + mod;
                    case 2: {
                        double mag = z.magnitude();
                        double angle = std::atan2(z.getImag(), z.getReal());
                        double newReal = mag * std::cos(angle * a) + c.getReal();
                        double newImag = mag * std::sin(angle * b) + c.getImag();
                        return Complex(newReal, newImag);
                    }
                    case 3:
                    default: {
                        // Inverse squared plus modification
                        double denom = z.magnitudeSquared();
                        if (denom < 1e-10) denom = 1e-10;
                        Complex inv(z.getReal() / denom, -z.getImag() / denom);
                        return inv * c + z2;
                    }
                }
            };
        }
            
        default:
            return mandelbrotFn;
    }
}

/**
 * @brief Compute fractal using the given iteration function
 * 
 * @param iterationFn Iteration function for the fractal
 * @param c Complex constant (for Mandelbrot and similar fractals)
 * @param z0 Initial value of z (for Julia sets)
 * @param maxIterations Maximum number of iterations
 * @param threshold Escape threshold
 * @return std::pair<int, std::vector<Complex>> 
 *         Pair of iteration count and history
 */
std::pair<int, std::vector<Complex>> 
computeFractal(
    const std::function<Complex(Complex, Complex)>& iterationFn,
    const Complex& c,
    const Complex& z0,
    int maxIterations,
    double threshold
) {
    std::vector<Complex> history;
    Complex z = z0;
    history.push_back(z);
    
    int iterations = 0;
    double thresholdSquared = threshold * threshold;
    
    while (iterations < maxIterations && z.magnitudeSquared() <= thresholdSquared) {
        z = iterationFn(z, c);
        history.push_back(z);
        iterations++;
    }
    
    return {iterations, history};
}

/**
 * @brief Generate parameters for Julia set exploration
 * 
 * @param seed Random seed for reproducibility
 * @return Complex Julia set parameter c
 */
Complex generateJuliaParameter(unsigned int seed = 42) {
    std::mt19937 gen(seed);
    
    // Generate parameters near known interesting Julia sets
    std::uniform_int_distribution<int> typeDist(0, 5);
    int type = typeDist(gen);
    
    std::uniform_real_distribution<double> smallVariation(-0.05, 0.05);
    double variation1 = smallVariation(gen);
    double variation2 = smallVariation(gen);
    
    switch (type) {
        case 0:  // Near dendrite Julia set
            return Complex(-0.75 + variation1, 0.1 + variation2);
            
        case 1:  // Near Douady's rabbit
            return Complex(-0.123 + variation1, 0.745 + variation2);
            
        case 2:  // Near San Marco fractal
            return Complex(-0.75 + variation1, 0.0 + variation2);
            
        case 3:  // Near Siegel disk
            return Complex(-0.391 + variation1, -0.587 + variation2);
            
        case 4:  // Near dendrite
            return Complex(0.285 + variation1, 0.0 + variation2);
            
        case 5:  // Random interesting point
        default: {
            std::uniform_real_distribution<double> realDist(-1.0, 0.5);
            std::uniform_real_distribution<double> imagDist(-1.0, 1.0);
            return Complex(realDist(gen), imagDist(gen));
        }
    }
}

/**
 * @brief Generate a description of the pattern generation parameters
 * 
 * @param mode Pattern generation mode
 * @param seed Random seed used
 * @param juliaParam Julia set parameter (if applicable)
 * @return std::string Description of the pattern parameters
 */
std::string generatePatternDescription(
    PatternMode mode,
    unsigned int seed,
    const Complex& juliaParam = Complex(0, 0)
) {
    std::string description;
    
    switch (mode) {
        case PatternMode::MANDELBROT:
            description = "Standard Mandelbrot Set: z = z² + c, starting with z = 0";
            break;
            
        case PatternMode::JULIA:
            description = "Julia Set with parameter c = " + 
                          std::to_string(juliaParam.getReal()) + " + " + 
                          std::to_string(juliaParam.getImag()) + "i";
            break;
            
        case PatternMode::BURNING_SHIP:
            description = "Burning Ship Fractal: z = (|Re(z)| + |Im(z)|i)² + c";
            break;
            
        case PatternMode::MULTIBROT: {
            // Determine power from seed
            std::mt19937 tempGen(seed);
            std::uniform_int_distribution<int> powerDist(3, 6);
            int power = powerDist(tempGen);
            
            description = "Multibrot Set of degree " + std::to_string(power) + 
                          ": z = z^" + std::to_string(power) + " + c";
            break;
        }
            
        case PatternMode::HYBRID:
            description = "Hybrid Fractal (Seed: " + std::to_string(seed) + 
                          "): Combination of multiple fractal formulas";
            break;
            
        case PatternMode::RANDOMIZED:
            description = "Randomized Fractal Pattern (Seed: " + std::to_string(seed) + 
                          "): Custom equation with random coefficients";
            break;
            
        default:
            description = "Unknown Pattern Type";
    }
    
    return description;
}