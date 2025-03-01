#include <utility>
#include <stdexcept>

/**
 * @brief Resolution settings for Mandelbrot set visualization
 */
struct Resolution {
    int width;      // Width in pixels
    int height;     // Height in pixels
    
    Resolution(int w = 800, int h = 600) : width(w), height(h) {
        validateResolution();
    }
    
    /**
     * Validate that resolution is within acceptable bounds
     */
    void validateResolution() const {
        if (width <= 0 || height <= 0) {
            throw std::invalid_argument("Resolution dimensions must be positive");
        }
        
        if (width > 10000 || height > 10000) {
            throw std::invalid_argument("Resolution is too large (max 10000x10000)");
        }
    }
};

/**
 * @brief Region of complex plane to visualize
 */
struct ComplexRegion {
    double xMin;    // Minimum real value
    double xMax;    // Maximum real value
    double yMin;    // Minimum imaginary value
    double yMax;    // Maximum imaginary value
    
    ComplexRegion(double x1 = -2.0, double x2 = 1.0, double y1 = -1.5, double y2 = 1.5)
        : xMin(x1), xMax(x2), yMin(y1), yMax(y2) {
        validateRegion();
    }
    
    /**
     * Validate that region is properly defined
     */
    void validateRegion() const {
        if (xMin >= xMax || yMin >= yMax) {
            throw std::invalid_argument("Invalid region bounds (min must be less than max)");
        }
    }
    
    /**
     * Calculate width of region
     */
    double width() const {
        return xMax - xMin;
    }
    
    /**
     * Calculate height of region
     */
    double height() const {
        return yMax - yMin;
    }
    
    /**
     * Calculate aspect ratio of region
     */
    double aspectRatio() const {
        return width() / height();
    }
};

/**
 * @brief Calculate pixel to complex plane mapping
 * 
 * @param resolution Screen resolution
 * @param region Complex plane region
 * @param pixelX X coordinate in pixels
 * @param pixelY Y coordinate in pixels
 * @return std::pair<double, double> Corresponding point (real, imag) in complex plane
 */
std::pair<double, double> pixelToComplex(
    const Resolution& resolution,
    const ComplexRegion& region,
    int pixelX,
    int pixelY
) {
    // Map pixel coordinates to complex plane coordinates
    double real = region.xMin + (pixelX * region.width() / resolution.width);
    double imag = region.yMax - (pixelY * region.height() / resolution.height);
    
    return std::make_pair(real, imag);
}

/**
 * @brief Adjust resolution to match the aspect ratio of the complex region
 * 
 * @param resolution Resolution to adjust (will be modified)
 * @param region Complex region
 * @param constrainWidth If true, keep width and adjust height, otherwise keep height and adjust width
 */
void adjustResolutionToRegion(
    Resolution& resolution,
    const ComplexRegion& region,
    bool constrainWidth = true
) {
    double regionAspect = region.aspectRatio();
    double resolutionAspect = static_cast<double>(resolution.width) / resolution.height;
    
    if (constrainWidth) {
        // Keep width, adjust height
        resolution.height = static_cast<int>(resolution.width / regionAspect);
    } else {
        // Keep height, adjust width
        resolution.width = static_cast<int>(resolution.height * regionAspect);
    }
    
    // Ensure we still have valid resolution
    resolution.validateResolution();
}

/**
 * @brief Zoom into a specific region of the complex plane
 * 
 * @param region Current region (will be modified)
 * @param centerX Center point X in current region coordinates
 * @param centerY Center point Y in current region coordinates
 * @param zoomFactor Factor to zoom by (> 1 for zoom in, < 1 for zoom out)
 */
void zoomRegion(
    ComplexRegion& region,
    double centerX,
    double centerY,
    double zoomFactor
) {
    if (zoomFactor <= 0) {
        throw std::invalid_argument("Zoom factor must be positive");
    }
    
    // Calculate new dimensions
    double newWidth = region.width() / zoomFactor;
    double newHeight = region.height() / zoomFactor;
    
    // Calculate new bounds centered on the given point
    region.xMin = centerX - newWidth / 2;
    region.xMax = centerX + newWidth / 2;
    region.yMin = centerY - newHeight / 2;
    region.yMax = centerY + newHeight / 2;
    
    // Validate the new region
    region.validateRegion();
}