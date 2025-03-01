#ifndef COMPLEX_H
#define COMPLEX_H

/**
 * @brief Complex number class for Mandelbrot set calculations
 * 
 * Implements complex number operations required for iterating the Mandelbrot equation z = z² + c
 */
class Complex {
private:
    double real;  // Real component
    double imag;  // Imaginary component

public:
    /**
     * @brief Construct a new Complex object
     * 
     * @param r Real component
     * @param i Imaginary component
     */
    Complex(double r = 0.0, double i = 0.0);

    /**
     * @brief Get the real component
     * 
     * @return double Real component value
     */
    double getReal() const;

    /**
     * @brief Get the imaginary component
     * 
     * @return double Imaginary component value
     */
    double getImag() const;

    /**
     * @brief Set the real component
     * 
     * @param r New real component value
     */
    void setReal(double r);

    /**
     * @brief Set the imaginary component
     * 
     * @param i New imaginary component value
     */
    void setImag(double i);

    /**
     * @brief Add two complex numbers
     * 
     * @param other Complex number to add
     * @return Complex Result of addition
     */
    Complex operator+(const Complex& other) const;

    /**
     * @brief Multiply two complex numbers
     * 
     * @param other Complex number to multiply with
     * @return Complex Result of multiplication
     */
    Complex operator*(const Complex& other) const;

    /**
     * @brief Calculate magnitude (absolute value) of complex number
     * 
     * @return double Magnitude of complex number
     */
    double magnitude() const;

    /**
     * @brief Calculate magnitude squared (avoids square root calculation)
     * 
     * @return double Square of magnitude
     */
    double magnitudeSquared() const;
};

#endif // COMPLEX_H