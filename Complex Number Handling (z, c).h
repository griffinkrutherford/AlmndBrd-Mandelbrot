#ifndef COMPLEX_NUMBER_HANDLING_H
#define COMPLEX_NUMBER_HANDLING_H

class Complex {
private:
    double real;
    double imag;

public:
    Complex(double r = 0.0, double i = 0.0);
    
    // Accessors
    double getReal() const;
    double getImag() const;
    
    // Mutators
    void setReal(double r);
    void setImag(double i);
    
    // Operations
    double magnitudeSquared() const;
    double magnitude() const;
    Complex operator+(const Complex& other) const;
    Complex operator*(const Complex& other) const;
};

#endif // COMPLEX_NUMBER_HANDLING_H
