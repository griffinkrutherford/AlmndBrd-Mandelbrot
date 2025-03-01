cmake_minimum_required(VERSION 3.14)
project(AlmndBrd-Mandelbrot VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable more warnings to catch potential issues
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# macOS specific - Find Homebrew packages
if(APPLE)
    # Try to find Homebrew prefix
    execute_process(
        COMMAND brew --prefix
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    # Use explicit paths for Homebrew packages
    set(SDL2_PATH "${HOMEBREW_PREFIX}/opt/sdl2")
    set(SDL2_TTF_PATH "${HOMEBREW_PREFIX}/opt/sdl2_ttf")
    
    message(STATUS "Homebrew prefix: ${HOMEBREW_PREFIX}")
    message(STATUS "SDL2 path: ${SDL2_PATH}")
    message(STATUS "SDL2_TTF path: ${SDL2_TTF_PATH}")
    
    # Add Homebrew include directories directly
    include_directories(
        "${SDL2_PATH}/include"
        "${SDL2_PATH}/include/SDL2"
        "${SDL2_TTF_PATH}/include"
        "${SDL2_TTF_PATH}/include/SDL2"
    )
    
    # Add Homebrew lib directories
    link_directories(
        "${SDL2_PATH}/lib"
        "${SDL2_TTF_PATH}/lib"
    )
endif()

# Find SDL2 (using find_package as fallback)
find_package(SDL2 QUIET)
if(NOT SDL2_FOUND)
    # If find_package fails, use direct link approach
    set(SDL2_LIBRARIES SDL2)
endif()

# Find SDL2_ttf (using find_package as fallback)
find_package(SDL2_ttf QUIET)
if(NOT SDL2_TTF_FOUND)
    # If find_package fails, use direct link approach
    set(SDL2_TTF_LIBRARIES SDL2_ttf)
endif()

# Include directories
include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/GenAI Mandelbrot Set Project/Code Structure (C++)/Base Class/Properties"
)

# Create a variable for your source files
set(SOURCE_FILES
    main.cpp
)

# Create the Complex Number header if it doesn't exist or replace with updated one
file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/Complex Number Handling (z, c).h"
"#ifndef COMPLEX_NUMBER_HANDLING_H
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
    
    // Mutators (setters)
    void setReal(double r);
    void setImag(double i);
    
    // Operations
    double magnitudeSquared() const;
    double magnitude() const;
    Complex operator+(const Complex& other) const;
    Complex operator*(const Complex& other) const;
};

#endif // COMPLEX_NUMBER_HANDLING_H
")

# Always use stub implementation to avoid mismatches
message(STATUS "Creating stub implementation for Complex class")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/Complex_stub.cpp" 
"#include \"Complex Number Handling (z, c).h\"
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
}")
list(APPEND SOURCE_FILES "${CMAKE_CURRENT_BINARY_DIR}/Complex_stub.cpp")

# Create ML prediction stub implementation
message(STATUS "Creating stub implementation for ML prediction")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ML_stub.cpp" 
"#include \"Complex Number Handling (z, c).h\"
#include <vector>
#include <utility>

enum class MLModelType { KNN, RANDOM_FOREST, NEURAL_NET, HYBRID, EVOLUTIONARY };

std::vector<std::pair<Complex, double>> predictInterestingRegions(
    MLModelType modelType,
    int trainingSize,
    double regionSize,
    int resolution
) { 
    // Return empty vector as stub implementation
    return {}; 
}")
list(APPEND SOURCE_FILES "${CMAKE_CURRENT_BINARY_DIR}/ML_stub.cpp")

add_executable(mandelbrot_app ${SOURCE_FILES})

# For macOS, directly link the libraries and frameworks
if(APPLE)
    target_link_libraries(mandelbrot_app PRIVATE
        SDL2
        SDL2_ttf
        "-framework Cocoa"
        "-framework Metal"
        "-framework CoreVideo"
        "-framework IOKit"
    )
else()
    target_link_libraries(mandelbrot_app PRIVATE
        ${SDL2_LIBRARIES}
        ${SDL2_TTF_LIBRARIES}
    )
endif()

# Optional: Display found libraries for debugging
message(STATUS "SDL2 libraries: ${SDL2_LIBRARIES}")
message(STATUS "SDL2_TTF libraries: ${SDL2_TTF_LIBRARIES}")