#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "../../Code Structure (C++)/Visualization Class/Resolution Settings.cpp"
#include "../../Code Structure (C++)/Main Program/Initialize Fractal Parameters.cpp"

/**
 * @file Resolution (Customizable).cpp
 * @brief Implementation of customizable resolution handling for fractal visualization
 * 
 * This file provides functionality for managing, optimizing, and selecting
 * appropriate resolutions for fractal rendering based on different requirements.
 */

/**
 * @brief Preset resolution profiles for different use cases
 */
enum class ResolutionProfile {
    PREVIEW,        // Quick previews during exploration
    STANDARD,       // Standard quality for normal viewing
    HIGH_QUALITY,   // High quality for detailed examination
    PRINT,          // Print-ready resolutions
    WALLPAPER,      // Common wallpaper resolutions
    ANIMATION,      // Resolutions suitable for animations
    CUSTOM          // User-defined custom resolution
};

/**
 * @brief Resolution constraints for optimization
 */
struct ResolutionConstraints {
    int minWidth;           // Minimum allowed width
    int minHeight;          // Minimum allowed height
    int maxWidth;           // Maximum allowed width
    int maxHeight;          // Maximum allowed height
    double minAspectRatio;  // Minimum aspect ratio (width/height)
    double maxAspectRatio;  // Maximum aspect ratio (width/height)
    int64_t maxPixels;      // Maximum total pixels (width * height)
    
    // Constructor with sensible defaults
    ResolutionConstraints(
        int minW = 320,
        int minH = 240,
        int maxW = 7680,    // 8K UHD
        int maxH = 4320,    // 8K UHD
        double minAR = 0.5, // Portrait orientation limit
        double maxAR = 3.0, // Wide panorama limit
        int64_t maxP = 33177600  // 8K UHD (7680x4320)
    ) : minWidth(minW),
        minHeight(minH),
        maxWidth(maxW),
        maxHeight(maxH),
        minAspectRatio(minAR),
        maxAspectRatio(maxAR),
        maxPixels(maxP) {
    }
};

/**
 * @brief Get a set of common resolutions for a profile
 * 
 * @param profile Resolution profile
 * @return std::vector<Resolution> List of common resolutions for the profile
 */
std::vector<Resolution> getProfileResolutions(ResolutionProfile profile) {
    std::vector<Resolution> resolutions;
    
    switch (profile) {
        case ResolutionProfile::PREVIEW:
            resolutions.push_back(Resolution(320, 240));   // QVGA
            resolutions.push_back(Resolution(640, 480));   // VGA
            resolutions.push_back(Resolution(800, 600));   // SVGA
            break;
            
        case ResolutionProfile::STANDARD:
            resolutions.push_back(Resolution(1024, 768));  // XGA
            resolutions.push_back(Resolution(1280, 720));  // HD 720p
            resolutions.push_back(Resolution(1366, 768));  // Common laptop
            resolutions.push_back(Resolution(1280, 1024)); // SXGA
            resolutions.push_back(Resolution(1600, 900));  // HD+
            resolutions.push_back(Resolution(1920, 1080)); // Full HD
            break;
            
        case ResolutionProfile::HIGH_QUALITY:
            resolutions.push_back(Resolution(1920, 1080)); // Full HD
            resolutions.push_back(Resolution(2560, 1440)); // QHD
            resolutions.push_back(Resolution(3440, 1440)); // Ultrawide QHD
            resolutions.push_back(Resolution(3840, 2160)); // 4K UHD
            break;
            
        case ResolutionProfile::PRINT:
            // Print at 300 DPI
            resolutions.push_back(Resolution(2480, 3508)); // A4 portrait
            resolutions.push_back(Resolution(3508, 2480)); // A4 landscape
            resolutions.push_back(Resolution(4960, 7016)); // A3 portrait (2x A4)
            resolutions.push_back(Resolution(7016, 4960)); // A3 landscape (2x A4)
            break;
            
        case ResolutionProfile::WALLPAPER:
            resolutions.push_back(Resolution(1920, 1080)); // Full HD
            resolutions.push_back(Resolution(2560, 1440)); // QHD
            resolutions.push_back(Resolution(3440, 1440)); // Ultrawide
            resolutions.push_back(Resolution(3840, 2160)); // 4K UHD
            resolutions.push_back(Resolution(5120, 1440)); // Super Ultrawide
            resolutions.push_back(Resolution(7680, 4320)); // 8K UHD
            break;
            
        case ResolutionProfile::ANIMATION:
            resolutions.push_back(Resolution(640, 360));   // nHD
            resolutions.push_back(Resolution(854, 480));   // FWVGA
            resolutions.push_back(Resolution(1280, 720));  // HD
            resolutions.push_back(Resolution(1920, 1080)); // Full HD
            break;
            
        case ResolutionProfile::CUSTOM:
        default:
            // Return empty for custom profile
            break;
    }
    
    return resolutions;
}

/**
 * @brief Optimize resolution based on constraints
 * 
 * @param desiredWidth Desired width
 * @param desiredHeight Desired height
 * @param constraints Resolution constraints
 * @return Resolution Optimized resolution
 */
Resolution optimizeResolution(
    int desiredWidth,
    int desiredHeight,
    const ResolutionConstraints& constraints
) {
    // Calculate aspect ratio
    double aspectRatio = static_cast<double>(desiredWidth) / desiredHeight;
    
    // Clamp aspect ratio to constraints
    aspectRatio = std::max(constraints.minAspectRatio, std::min(aspectRatio, constraints.maxAspectRatio));
    
    // Clamp width and height to min/max constraints
    int width = std::max(constraints.minWidth, std::min(desiredWidth, constraints.maxWidth));
    int height = std::max(constraints.minHeight, std::min(desiredHeight, constraints.maxHeight));
    
    // Check if total pixels exceed maximum
    int64_t totalPixels = static_cast<int64_t>(width) * height;
    
    if (totalPixels > constraints.maxPixels) {
        // Scale down while preserving aspect ratio
        double scale = std::sqrt(static_cast<double>(constraints.maxPixels) / totalPixels);
        width = static_cast<int>(width * scale);
        height = static_cast<int>(height * scale);
        
        // Ensure minimum dimensions are still met
        width = std::max(constraints.minWidth, width);
        height = std::max(constraints.minHeight, height);
    }
    
    // Ensure the width and height are even numbers (some algorithms prefer this)
    width = (width / 2) * 2;
    height = (height / 2) * 2;
    
    return Resolution(width, height);
}

/**
 * @brief Find the best matching resolution from a profile
 * 
 * @param desiredWidth Desired width
 * @param desiredHeight Desired height
 * @param profile Resolution profile
 * @return Resolution Best matching resolution
 */
Resolution findBestMatchingResolution(
    int desiredWidth,
    int desiredHeight,
    ResolutionProfile profile
) {
    std::vector<Resolution> profileResolutions = getProfileResolutions(profile);
    
    if (profileResolutions.empty()) {
        // Return custom resolution if profile has no presets
        return Resolution(desiredWidth, desiredHeight);
    }
    
    // Find the closest resolution by minimizing the difference in pixel count
    int64_t desiredPixels = static_cast<int64_t>(desiredWidth) * desiredHeight;
    double desiredAspectRatio = static_cast<double>(desiredWidth) / desiredHeight;
    
    // Best match so far
    Resolution bestMatch = profileResolutions[0];
    double bestMatchScore = std::numeric_limits<double>::max();
    
    for (const auto& res : profileResolutions) {
        int64_t pixelDiff = std::abs(static_cast<int64_t>(res.width) * res.height - desiredPixels);
        double aspectRatio = static_cast<double>(res.width) / res.height;
        double aspectRatioDiff = std::abs(aspectRatio - desiredAspectRatio);
        
        // Create a combined score giving more weight to aspect ratio than total pixels
        double score = aspectRatioDiff * 10.0 + pixelDiff / 1000000.0;
        
        if (score < bestMatchScore) {
            bestMatchScore = score;
            bestMatch = res;
        }
    }
    
    return bestMatch;
}

/**
 * @brief Adjust resolution to match aspect ratio of the complex region
 * 
 * @param resolution Resolution to adjust
 * @param region Complex region
 * @param constrainWidth If true, keep width and adjust height, otherwise keep height and adjust width
 * @return Resolution Adjusted resolution
 */
Resolution adjustResolutionToRegion(
    const Resolution& resolution,
    const ComplexRegion& region,
    bool constrainWidth = true
) {
    double regionAspect = region.aspectRatio();
    double resolutionAspect = static_cast<double>(resolution.width) / resolution.height;
    
    Resolution adjusted = resolution;
    
    if (constrainWidth) {
        // Keep width, adjust height
        adjusted.height = static_cast<int>(adjusted.width / regionAspect);
    } else {
        // Keep height, adjust width
        adjusted.width = static_cast<int>(adjusted.height * regionAspect);
    }
    
    // Ensure values are even
    adjusted.width = (adjusted.width / 2) * 2;
    adjusted.height = (adjusted.height / 2) * 2;
    
    return adjusted;
}

/**
 * @brief Calculate optimal resolution based on available memory
 * 
 * @param availableMemoryMB Available memory in MB
 * @param bytesPerPixel Memory usage per pixel (typically 4 for RGBA)
 * @param memoryUsageRatio Ratio of memory to allocate (0.0-1.0)
 * @param preferredAspectRatio Preferred aspect ratio
 * @return Resolution Optimal resolution based on memory constraints
 */
Resolution calculateMemoryConstrainedResolution(
    int availableMemoryMB,
    int bytesPerPixel = 4,
    double memoryUsageRatio = 0.5,
    double preferredAspectRatio = 16.0/9.0
) {
    // Calculate total available pixels
    int64_t availableBytes = static_cast<int64_t>(availableMemoryMB) * 1024 * 1024 * memoryUsageRatio;
    int64_t maxPixels = availableBytes / bytesPerPixel;
    
    // Calculate optimal dimensions while preserving aspect ratio
    int width = static_cast<int>(std::sqrt(maxPixels * preferredAspectRatio));
    int height = static_cast<int>(width / preferredAspectRatio);
    
    // Ensure the width and height are even numbers
    width = (width / 2) * 2;
    height = (height / 2) * 2;
    
    return Resolution(width, height);
}

/**
 * @brief Create a custom resolution based on specific output requirements
 * 
 * @param outputType Type of output (e.g., "screen", "print", "poster")
 * @param width Width in output units (pixels, inches, etc.)
 * @param height Height in output units
 * @param dpi DPI for print media (ignored for screen)
 * @return Resolution Custom resolution
 */
Resolution createCustomResolution(
    const std::string& outputType,
    double width,
    double height,
    int dpi = 300
) {
    if (outputType == "screen") {
        // Screen resolution is already in pixels
        return Resolution(
            static_cast<int>(width),
            static_cast<int>(height)
        );
    } else if (outputType == "print" || outputType == "poster") {
        // Convert inches to pixels using DPI
        return Resolution(
            static_cast<int>(std::ceil(width * dpi)),
            static_cast<int>(std::ceil(height * dpi))
        );
    } else {
        throw std::invalid_argument("Unknown output type: " + outputType);
    }
}

/**
 * @brief Suggest a resolution for rendering a specific region
 * 
 * @param region Complex region to render
 * @param profile Resolution profile
 * @param detailLevel Detail level settings
 * @return Resolution Suggested resolution
 */
Resolution suggestResolutionForRegion(
    const ComplexRegion& region,
    ResolutionProfile profile,
    DetailLevel detailLevel
) {
    // Base resolution depending on detail level
    Resolution baseResolution;
    
    switch (detailLevel) {
        case DetailLevel::DRAFT:
            baseResolution = getResolutionPreset(ResolutionPreset::TINY);
            break;
        case DetailLevel::LOW:
            baseResolution = getResolutionPreset(ResolutionPreset::SMALL);
            break;
        case DetailLevel::MEDIUM:
            baseResolution = getResolutionPreset(ResolutionPreset::MEDIUM);
            break;
        case DetailLevel::HIGH:
            baseResolution = getResolutionPreset(ResolutionPreset::HD);
            break;
        case DetailLevel::ULTRA:
            baseResolution = getResolutionPreset(ResolutionPreset::FULL_HD);
            break;
        default:
            baseResolution = getResolutionPreset(ResolutionPreset::MEDIUM);
    }
    
    // Adjust for complexity of the region (smaller regions need higher resolution)
    double regionSize = region.width() * region.height();
    double complexityFactor = 1.0;
    
    // Standard size for comparison (full Mandelbrot view)
    const double STANDARD_SIZE = 7.5; // (-2.5 to 1.0) * (-1.5 to 1.5)
    
    if (regionSize < STANDARD_SIZE) {
        // Increase resolution for smaller regions (zoomed in)
        complexityFactor = std::min(3.0, STANDARD_SIZE / regionSize);
    }
    
    // Scale the base resolution
    int scaledWidth = static_cast<int>(baseResolution.width * std::sqrt(complexityFactor));
    int scaledHeight = static_cast<int>(baseResolution.height * std::sqrt(complexityFactor));
    
    // Find the best matching resolution from the profile
    Resolution matchedResolution = findBestMatchingResolution(
        scaledWidth, scaledHeight, profile
    );
    
    // Adjust to match the aspect ratio of the region
    Resolution finalResolution = adjustResolutionToRegion(
        matchedResolution, region
    );
    
    return finalResolution;
}

/**
 * @brief Estimate time to render based on resolution and detail level
 * 
 * @param resolution Image resolution
 * @param detailLevel Detail level settings
 * @return std::pair<double, std::string> Estimated time in seconds and humanized string
 */
std::pair<double, std::string> estimateRenderTime(
    const Resolution& resolution,
    DetailLevel detailLevel
) {
    // Constants for estimation (these would be calibrated based on system)
    double pixelsPerSecondBase = 1000000.0; // 1 million pixels per second baseline
    
    // Detail level factor
    double detailFactor;
    switch (detailLevel) {
        case DetailLevel::DRAFT:   detailFactor = 0.2; break;
        case DetailLevel::LOW:     detailFactor = 0.5; break;
        case DetailLevel::MEDIUM:  detailFactor = 1.0; break;
        case DetailLevel::HIGH:    detailFactor = 2.0; break;
        case DetailLevel::ULTRA:   detailFactor = 5.0; break;
        default:                   detailFactor = 1.0;
    }
    
    // Calculate pixels and time
    int64_t totalPixels = static_cast<int64_t>(resolution.width) * resolution.height;
    double timeSeconds = totalPixels / (pixelsPerSecondBase / detailFactor);
    
    // Generate human-readable time string
    std::string timeString;
    if (timeSeconds < 1.0) {
        // Milliseconds
        int ms = static_cast<int>(timeSeconds * 1000);
        timeString = std::to_string(ms) + " milliseconds";
    } else if (timeSeconds < 60.0) {
        // Seconds
        int sec = static_cast<int>(timeSeconds);
        timeString = std::to_string(sec) + " seconds";
    } else if (timeSeconds < 3600.0) {
        // Minutes and seconds
        int min = static_cast<int>(timeSeconds / 60);
        int sec = static_cast<int>(timeSeconds) % 60;
        timeString = std::to_string(min) + " minutes " + std::to_string(sec) + " seconds";
    } else {
        // Hours and minutes
        int hour = static_cast<int>(timeSeconds / 3600);
        int min = (static_cast<int>(timeSeconds) % 3600) / 60;
        timeString = std::to_string(hour) + " hours " + std::to_string(min) + " minutes";
    }
    
    return {timeSeconds, timeString};
}

/**
 * @brief Print resolution information table
 * 
 * @param resolutions Vector of resolutions
 * @param showDetails Whether to show detailed information for each resolution
 */
void printResolutionTable(
    const std::vector<Resolution>& resolutions,
    bool showDetails = true
) {
    std::cout << "=== Available Resolutions ===" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::left << std::setw(5) << "ID" 
              << std::setw(15) << "Dimensions" 
              << std::setw(15) << "Aspect Ratio" 
              << std::setw(15) << "Total Pixels" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (size_t i = 0; i < resolutions.size(); i++) {
        const auto& res = resolutions[i];
        double aspectRatio = static_cast<double>(res.width) / res.height;
        int64_t totalPixels = static_cast<int64_t>(res.width) * res.height;
        
        std::cout << std::left << std::setw(5) << i
                  << std::setw(15) << (std::to_string(res.width) + "x" + std::to_string(res.height))
                  << std::setw(15) << std::fixed << std::setprecision(2) << aspectRatio
                  << std::setw(15) << totalPixels << std::endl;
        
        if (showDetails) {
            // Show estimated render times for different detail levels
            auto draftTime = estimateRenderTime(res, DetailLevel::DRAFT);
            auto mediumTime = estimateRenderTime(res, DetailLevel::MEDIUM);
            auto ultraTime = estimateRenderTime(res, DetailLevel::ULTRA);
            
            std::cout << "    Render time estimates:" << std::endl;
            std::cout << "    - Draft: " << draftTime.second << std::endl;
            std::cout << "    - Medium: " << mediumTime.second << std::endl;
            std::cout << "    - Ultra: " << ultraTime.second << std::endl;
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Main function to demonstrate resolution customization
 */
int main() {
    // Show available resolution profiles
    std::cout << "=== Resolution Profiles ===" << std::endl;
    std::cout << "1. Preview" << std::endl;
    std::cout << "2. Standard" << std::endl;
    std::cout << "3. High Quality" << std::endl;
    std::cout << "4. Print" << std::endl;
    std::cout << "5. Wallpaper" << std::endl;
    std::cout << "6. Animation" << std::endl;
    std::cout << std::endl;
    
    // Get profile resolutions
    auto previewResolutions = getProfileResolutions(ResolutionProfile::PREVIEW);
    auto standardResolutions = getProfileResolutions(ResolutionProfile::STANDARD);
    auto highQualityResolutions = getProfileResolutions(ResolutionProfile::HIGH_QUALITY);
    
    // Print resolution tables
    std::cout << "Preview Resolutions:" << std::endl;
    printResolutionTable(previewResolutions, false);
    std::cout << std::endl;
    
    std::cout << "Standard Resolutions:" << std::endl;
    printResolutionTable(standardResolutions);
    std::cout << std::endl;
    
    // Demonstrate resolution optimization
    std::cout << "=== Resolution Optimization Examples ===" << std::endl;
    
    // Define constraints
    ResolutionConstraints constraints;
    constraints.maxPixels = 4000000;  // Limit to 4MP
    
    // Try to optimize a large resolution
    Resolution desired(3840, 2160);  // 4K UHD
    Resolution optimized = optimizeResolution(desired.width, desired.height, constraints);
    
    std::cout << "Desired: " << desired.width << "x" << desired.height 
              << " (" << (static_cast<int64_t>(desired.width) * desired.height) << " pixels)" << std::endl;
    std::cout << "Optimized: " << optimized.width << "x" << optimized.height 
              << " (" << (static_cast<int64_t>(optimized.width) * optimized.height) << " pixels)" << std::endl;
    std::cout << std::endl;
    
    // Demonstrate region-specific resolution suggestion
    std::cout << "=== Region-Specific Resolution Suggestions ===" << std::endl;
    
    // Define test regions
    ComplexRegion fullView(-2.5, 1.0, -1.5, 1.5);
    ComplexRegion cardioid(-1.2, -0.4, -0.4, 0.4);
    ComplexRegion deepZoom(-0.743643887037158, -0.743643887037157, 0.131825904205311, 0.131825904205312);
    
    std::cout << "Full View Region:" << std::endl;
    std::cout << "  Medium Detail: " 
              << suggestResolutionForRegion(fullView, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).width << "x"
              << suggestResolutionForRegion(fullView, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).height << std::endl;
    std::cout << "  Ultra Detail: " 
              << suggestResolutionForRegion(fullView, ResolutionProfile::STANDARD, DetailLevel::ULTRA).width << "x"
              << suggestResolutionForRegion(fullView, ResolutionProfile::STANDARD, DetailLevel::ULTRA).height << std::endl;
              
    std::cout << "Main Cardioid Region:" << std::endl;
    std::cout << "  Medium Detail: " 
              << suggestResolutionForRegion(cardioid, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).width << "x"
              << suggestResolutionForRegion(cardioid, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).height << std::endl;
              
    std::cout << "Deep Zoom Region:" << std::endl;
    std::cout << "  Medium Detail: " 
              << suggestResolutionForRegion(deepZoom, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).width << "x"
              << suggestResolutionForRegion(deepZoom, ResolutionProfile::STANDARD, DetailLevel::MEDIUM).height << std::endl;
    std::cout << "  High Detail: " 
              << suggestResolutionForRegion(deepZoom, ResolutionProfile::HIGH_QUALITY, DetailLevel::HIGH).width << "x"
              << suggestResolutionForRegion(deepZoom, ResolutionProfile::HIGH_QUALITY, DetailLevel::HIGH).height << std::endl;
    std::cout << std::endl;
    
    // Demonstrate memory-constrained resolution calculation
    std::cout << "=== Memory-Constrained Resolutions ===" << std::endl;
    
    std::vector<int> memorySizes = {1024, 2048, 4096, 8192};  // MB
    
    for (int memory : memorySizes) {
        Resolution memoryConstrained = calculateMemoryConstrainedResolution(memory);
        std::cout << "Available Memory: " << memory << " MB" << std::endl;
        std::cout << "  Optimal Resolution: " << memoryConstrained.width << "x" << memoryConstrained.height << std::endl;
        std::cout << "  Total Pixels: " << (static_cast<int64_t>(memoryConstrained.width) * memoryConstrained.height) << std::endl;
        
        // Estimate render time
        auto renderTime = estimateRenderTime(memoryConstrained, DetailLevel::MEDIUM);
        std::cout << "  Estimated Render Time (Medium): " << renderTime.second << std::endl;
        std::cout << std::endl;
    }
    
    // Demonstrate print resolution calculation
    std::cout << "=== Print Resolution Examples ===" << std::endl;
    
    Resolution a4Portrait = createCustomResolution("print", 8.27, 11.69);  // A4 paper in inches
    Resolution posterLarge = createCustomResolution("print", 24, 36, 150);  // Poster with lower DPI
    
    std::cout << "A4 Print (300 DPI): " << a4Portrait.width << "x" << a4Portrait.height << std::endl;
    std::cout << "Large Poster (150 DPI): " << posterLarge.width << "x" << posterLarge.height << std::endl;
    
    return 0;
}