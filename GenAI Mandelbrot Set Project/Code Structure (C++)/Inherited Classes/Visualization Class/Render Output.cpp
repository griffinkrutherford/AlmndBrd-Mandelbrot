#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <cstdlib>
#include "Resolution Settings.cpp"
#include "Color Assignment Logic.cpp"
#include "../../Base Class/Methods/Iterate Equation.cpp"
#include "../../Base Class/Methods/Check Divergence.cpp"

/**
 * @brief Supported image output formats
 */
enum class ImageFormat {
    PPM,    // Portable Pixmap Format (simple text-based format)
    BMP,    // Bitmap format
    PNG,    // PNG format (requires external library in full implementation)
    JPG     // JPEG format (requires external library in full implementation)
};

/**
 * @brief Settings for rendering a Mandelbrot set
 */
struct RenderSettings {
    Resolution resolution;             // Image resolution
    ComplexRegion region;              // Region of complex plane
    IterationParameters iterParams;    // Iteration parameters
    ColorMode colorMode;               // Color mapping mode
    ImageFormat format;                // Output format
    
    RenderSettings(
        Resolution res = Resolution(800, 600),
        ComplexRegion reg = ComplexRegion(-2.0, 1.0, -1.5, 1.5),
        IterationParameters iter = IterationParameters(1000, 2.0),
        ColorMode mode = ColorMode::LOGARITHMIC,
        ImageFormat fmt = ImageFormat::PPM
    ) : resolution(res), region(reg), iterParams(iter), colorMode(mode), format(fmt) {}
};

/**
 * @brief Generate filename for the output image
 * 
 * @param format Image format
 * @return std::string Filename with timestamp
 */
std::string generateFilename(ImageFormat format) {
    // Get current time for unique filename
    std::time_t now = std::time(nullptr);
    char timestamp[16];
    std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", std::localtime(&now));
    
    std::string extension;
    switch (format) {
        case ImageFormat::PPM: extension = ".ppm"; break;
        case ImageFormat::BMP: extension = ".bmp"; break;
        case ImageFormat::PNG: extension = ".png"; break;
        case ImageFormat::JPG: extension = ".jpg"; break;
    }
    
    return "mandelbrot_" + std::string(timestamp) + extension;
}

/**
 * @brief Compute the Mandelbrot set and return iteration counts for each pixel
 * 
 * @param settings Render settings
 * @return std::vector<std::vector<int>> 2D array of iteration counts
 */
std::vector<std::vector<int>> computeMandelbrotData(const RenderSettings& settings) {
    const Resolution& res = settings.resolution;
    const ComplexRegion& region = settings.region;
    const IterationParameters& params = settings.iterParams;
    
    // Initialize result array
    std::vector<std::vector<int>> iterationData(
        res.height, std::vector<int>(res.width, 0)
    );
    
    // Compute Mandelbrot set for each pixel
    for (int y = 0; y < res.height; y++) {
        for (int x = 0; x < res.width; x++) {
            // Map pixel to complex plane
            auto [real, imag] = pixelToComplex(res, region, x, y);
            
            // Create complex number for this point
            Complex c(real, imag);
            
            // Compute iterations until divergence or max iterations
            iterationData[y][x] = computeFull(c, params.maxIterations, params.threshold);
        }
    }
    
    return iterationData;
}

/**
 * @brief Write iteration data to a PPM image file
 * 
 * @param filename Output filename
 * @param iterationData 2D array of iteration counts
 * @param settings Render settings
 */
void writePPMImage(
    const std::string& filename,
    const std::vector<std::vector<int>>& iterationData,
    const RenderSettings& settings
) {
    const Resolution& res = settings.resolution;
    
    // Generate color palette based on color mode
    std::vector<Color> palette = generateColorPalette(
        settings.colorMode, settings.iterParams.maxIterations
    );
    
    // Open output file
    std::ofstream outFile(filename, std::ios::out);
    if (!outFile) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    // Write PPM header
    outFile << "P3\n"; // ASCII RGB format
    outFile << res.width << " " << res.height << "\n";
    outFile << "255\n"; // Max color value
    
    // Write pixel data
    for (int y = 0; y < res.height; y++) {
        for (int x = 0; x < res.width; x++) {
            int iterations = iterationData[y][x];
            const Color& color = palette[iterations];
            
            outFile << static_cast<int>(color.r) << " "
                    << static_cast<int>(color.g) << " "
                    << static_cast<int>(color.b) << " ";
        }
        outFile << "\n";
    }
    
    outFile.close();
}

/**
 * @brief Write iteration data to a BMP image file
 * 
 * @param filename Output filename
 * @param iterationData 2D array of iteration counts
 * @param settings Render settings
 */
void writeBMPImage(
    const std::string& filename,
    const std::vector<std::vector<int>>& iterationData,
    const RenderSettings& settings
) {
    const Resolution& res = settings.resolution;
    
    // Generate color palette based on color mode
    std::vector<Color> palette = generateColorPalette(
        settings.colorMode, settings.iterParams.maxIterations
    );
    
    // Open output file
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    // Calculate row size and padding
    int rowSize = res.width * 3; // 3 bytes per pixel (RGB)
    int padding = (4 - (rowSize % 4)) % 4; // Row size must be multiple of 4
    int paddedRowSize = rowSize + padding;
    
    // Calculate file size
    int fileSize = 54 + (paddedRowSize * res.height); // 54 bytes for header
    
    // BMP Header (14 bytes)
    unsigned char bmpHeader[14] = {
        'B', 'M',               // Signature
        0, 0, 0, 0,             // File size (bytes 2-5)
        0, 0,                   // Reserved
        0, 0,                   // Reserved
        54, 0, 0, 0             // Pixel data offset
    };
    
    // Update file size in header
    bmpHeader[2] = fileSize & 0xFF;
    bmpHeader[3] = (fileSize >> 8) & 0xFF;
    bmpHeader[4] = (fileSize >> 16) & 0xFF;
    bmpHeader[5] = (fileSize >> 24) & 0xFF;
    
    // DIB Header (40 bytes)
    unsigned char dibHeader[40] = {
        40, 0, 0, 0,            // DIB header size
        0, 0, 0, 0,             // Width (bytes 4-7)
        0, 0, 0, 0,             // Height (bytes 8-11)
        1, 0,                   // Color planes
        24, 0,                  // Bits per pixel (24 for RGB)
        0, 0, 0, 0,             // No compression
        0, 0, 0, 0,             // Image size (can be 0 for RGB)
        0, 0, 0, 0,             // X pixels per meter
        0, 0, 0, 0,             // Y pixels per meter
        0, 0, 0, 0,             // Colors in color table
        0, 0, 0, 0              // Important color count
    };
    
    // Update width and height in DIB header
    dibHeader[4] = res.width & 0xFF;
    dibHeader[5] = (res.width >> 8) & 0xFF;
    dibHeader[6] = (res.width >> 16) & 0xFF;
    dibHeader[7] = (res.width >> 24) & 0xFF;
    
    dibHeader[8] = res.height & 0xFF;
    dibHeader[9] = (res.height >> 8) & 0xFF;
    dibHeader[10] = (res.height >> 16) & 0xFF;
    dibHeader[11] = (res.height >> 24) & 0xFF;
    
    // Write headers
    outFile.write(reinterpret_cast<char*>(bmpHeader), sizeof(bmpHeader));
    outFile.write(reinterpret_cast<char*>(dibHeader), sizeof(dibHeader));
    
    // Write pixel data (BMP stores rows from bottom to top)
    unsigned char paddingBytes[3] = {0, 0, 0};
    
    for (int y = res.height - 1; y >= 0; y--) {
        for (int x = 0; x < res.width; x++) {
            int iterations = iterationData[y][x];
            const Color& color = palette[iterations];
            
            // BMP stores colors as BGR
            unsigned char pixel[3] = {
                color.b,
                color.g,
                color.r
            };
            
            outFile.write(reinterpret_cast<char*>(pixel), 3);
        }
        
        // Write padding bytes
        if (padding > 0) {
            outFile.write(reinterpret_cast<char*>(paddingBytes), padding);
        }
    }
    
    outFile.close();
}

/**
 * @brief Render Mandelbrot set and save to image file
 * 
 * @param settings Render settings
 * @return std::string Path to generated image file
 */
std::string renderMandelbrotSet(const RenderSettings& settings) {
    // Generate filename
    std::string filename = generateFilename(settings.format);
    
    // Compute Mandelbrot set
    std::vector<std::vector<int>> iterationData = computeMandelbrotData(settings);
    
    // Write to file based on format
    switch (settings.format) {
        case ImageFormat::PPM:
            writePPMImage(filename, iterationData, settings);
            break;
        case ImageFormat::BMP:
            writeBMPImage(filename, iterationData, settings);
            break;
        case ImageFormat::PNG:
        case ImageFormat::JPG:
            throw std::runtime_error(
                "PNG and JPG formats require external libraries and are not implemented in this example"
            );
    }
    
    return filename;
}

/**
 * @brief Display text-based representation of Mandelbrot set to console
 * 
 * @param settings Render settings (resolution should be small for console output)
 */
void displayTextMandelbrot(const RenderSettings& settings) {
    // Character set from dense to sparse
    const std::string CHARSET = " .:-=+*#%@";
    
    // Compute Mandelbrot set
    std::vector<std::vector<int>> iterationData = computeMandelbrotData(settings);
    
    // Display as text
    for (int y = 0; y < settings.resolution.height; y++) {
        for (int x = 0; x < settings.resolution.width; x++) {
            int iterations = iterationData[y][x];
            
            char displayChar;
            if (iterations == settings.iterParams.maxIterations) {
                // Point is likely in the Mandelbrot set
                displayChar = CHARSET.back();
            } else {
                // Map iterations to character
                int index = static_cast<int>(
                    static_cast<double>(iterations) / settings.iterParams.maxIterations * (CHARSET.length() - 1)
                );
                displayChar = CHARSET[index];
            }
            
            std::cout << displayChar;
        }
        std::cout << std::endl;
    }
}