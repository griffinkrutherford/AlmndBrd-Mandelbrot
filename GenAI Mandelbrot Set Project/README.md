# AlmndBrd-Mandelbrot

A Mandelbrot set explorer with Generative AI integration.

## Project Overview

AlmndBrd-Mandelbrot is a C++ application that combines traditional Mandelbrot set exploration with innovative Generative AI techniques. It provides tools for visualizing, analyzing, and generating novel fractal patterns.

## Features

- High-performance Mandelbrot set calculation
- Multiple visualization methods for both interior and exterior regions
- Resolution optimization based on system capabilities
- Machine learning integration for interesting pattern detection
- Support for various image output formats
- Customizable color mapping with numerous color schemes

## Building the Project

### Prerequisites

- CMake 3.14 or higher
- C++17 compatible compiler (GCC 9+, Clang 10+, or MSVC 2019+)
- Optional: OpenMP for parallelization
- Optional: libpng, libjpeg for PNG/JPEG output

### Build Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/AlmndBrd-Mandelbrot.git
   cd AlmndBrd-Mandelbrot
   ```

2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

3. Configure with CMake:
   ```bash
   cmake ..
   ```

4. Build the project:
   ```bash
   cmake --build .
   ```

### Optional CMake Flags

- `-DUSE_OPENMP=OFF` to disable OpenMP parallelization
- `-DBUILD_TESTS=ON` to build the test suite

## Usage

Run the main application:

```bash
./almnd_brd_app --x=-0.75 --y=0.0 --zoom=4.0 --width=1920 --height=1080 --iter=1000 --out=my_fractal
```

Or use a preset:

```bash
./almnd_brd_app --preset=spiral --width=1920 --height=1080 --iter=2000
```

## Available Presets

- `full`: Full Mandelbrot set view
- `cardioid`: Close-up of the main cardioid
- `seahorse`: Seahorse valley
- `spiral`: Deep zoom into a spiral pattern

## Project Structure

- **Base Class**: Core fractal computation
- **Visualization Class**: Rendering and color mapping
- **GenAI Class**: AI integration for fractal exploration
- **Main Program**: Application entry points
- **Output**: Image output and customization
- **Concept**: Mathematical basis and AI integration concepts

## Running Tests

If you built with `-DBUILD_TESTS=ON`, you can run the tests:

```bash
ctest
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The Mandelbrot set was discovered by Benoît Mandelbrot in 1980
- Implementation inspired by numerous fractal visualization techniques and tools