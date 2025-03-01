#include <iostream>
#include <vector>
#include <cmath>

/**
 * Simple Complex number class for the Mandelbrot set
 */
class Complex {
private:
    double real;
    double imag;

public:
    Complex(double r = 0.0, double i = 0.0) : real(r), imag(i) {}
    
    double getReal() const { return real; }
    double getImag() const { return imag; }
    
    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }
    
    Complex operator*(const Complex& other) const {
        return Complex(
            real * other.real - imag * other.imag,
            real * other.imag + imag * other.real
        );
    }
    
    double magnitudeSquared() const {
        return real * real + imag * imag;
    }
};

/**
 * Simple Mandelbrot set renderer
 */
int main() {
    std::cout << "AlmndBrd-Mandelbrot Explorer (Simplified)" << std::endl;
    
    // Parameters
    const int WIDTH = 80;
    const int HEIGHT = 40;
    const int MAX_ITERATIONS = 100;
    const double X_MIN = -2.0;
    const double X_MAX = 1.0;
    const double Y_MIN = -1.5;
    const double Y_MAX = 1.5;
    
    // Calculate Mandelbrot set
    std::vector<std::vector<int>> iterations(HEIGHT, std::vector<int>(WIDTH));
    
    for (int y = 0; y < HEIGHT; y++) {
        double imag = Y_MAX - y * (Y_MAX - Y_MIN) / HEIGHT;
        
        for (int x = 0; x < WIDTH; x++) {
            double real = X_MIN + x * (X_MAX - X_MIN) / WIDTH;
            
            Complex c(real, imag);
            Complex z(0, 0);
            int i = 0;
            
            while (i < MAX_ITERATIONS && z.magnitudeSquared() < 4.0) {
                z = z * z + c;
                i++;
            }
            
            iterations[y][x] = i;
        }
    }
    
    // Display as ASCII art
    const char* CHARSET = " .:-=+*#%@";
    const int CHARSET_LENGTH = 10;
    
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int i = iterations[y][x];
            
            if (i == MAX_ITERATIONS) {
                std::cout << CHARSET[CHARSET_LENGTH - 1];
            } else {
                int index = i * (CHARSET_LENGTH - 1) / MAX_ITERATIONS;
                std::cout << CHARSET[index];
            }
        }
        std::cout << std::endl;
    }
    
    return 0;
}
