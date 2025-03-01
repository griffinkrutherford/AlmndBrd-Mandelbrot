#ifndef ITERATION_PARAMETERS_H
#define ITERATION_PARAMETERS_H

/**
 * @brief Configuration parameters for Mandelbrot set iteration
 * 
 * Contains threshold and iteration limit parameters for the fractal computation
 */
struct IterationParameters {
    int maxIterations;    // Maximum number of iterations to perform
    double threshold;     // Escape threshold for divergence detection
    
    /**
     * @brief Construct iteration parameters with default values
     * 
     * @param maxIter Maximum number of iterations (default: 100)
     * @param thresh Escape threshold (default: 2.0)
     */
    IterationParameters(int maxIter = 100, double thresh = 2.0) 
        : maxIterations(maxIter), threshold(thresh) {}
};

#endif // ITERATION_PARAMETERS_H