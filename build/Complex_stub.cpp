#include "Complex Number Handling (z, c).h"
#include <cmath>

Complex::Complex(double r, double i) : real(r), imag(i) {}

double Complex::getReal() const { 
    return real; 
}

double Complex::getImag() const { 
    return imag; 
}

void Complex::setReal(double r) {
    real = r;
}

void Complex::setImag(double i) {
    imag = i;
}

double Complex::magnitudeSquared() const { 
    return real*real + imag*imag; 
}

double Complex::magnitude() const {
    return std::sqrt(magnitudeSquared());
}

Complex Complex::operator+(const Complex& other) const { 
    return Complex(real + other.real, imag + other.imag); 
}

Complex Complex::operator*(const Complex& other) const { 
    return Complex(real * other.real - imag * other.imag, 
                   real * other.imag + imag * other.real); 
}