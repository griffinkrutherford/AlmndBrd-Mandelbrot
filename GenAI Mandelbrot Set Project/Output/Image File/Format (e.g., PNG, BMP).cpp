#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <stdexcept>
#include "../../Code Structure (C++)/Visualization Class/Resolution Settings.cpp"
#include "../../Code Structure (C++)/Visualization Class/Color Assignment Logic.cpp"

/**
 * @file Format (e.g., PNG, BMP).cpp
 * @brief Implementation of image format conversion and output for fractal visualization
 * 
 * This file provides functionality for saving fractal visualization data to
 * different image formats, focusing on BMP and PPM with notes about PNG integration.
 */

/**
 * @brief Available image format types
 */
enum class ImageFormatType {
    BMP,        // Windows Bitmap
    PNG,        // Portable Network Graphics
    JPG,        // JPEG
    PPM,        // Portable Pixmap (text-based)
    TIFF,       // Tagged Image File Format
    EXR         // OpenEXR (HDR format)
};

/**
 * @brief Format-specific options for image output
 */
struct FormatOptions {
    // General options
    bool includeMetadata;       // Whether to include fractal parameters as metadata
    bool progressiveLoading;    // For formats that support progressive loading
    
    // Format-specific options
    union {
        struct {
            int compressionLevel;   // 0-9, where 9 is maximum compression
            bool interlaced;        // Whether to use interlaced encoding
        } png;
        
        struct {
            int quality;            // 0-100, where 100 is highest quality
            bool progressive;       // Whether to use progressive encoding
        } jpg;
        
        struct {
            bool rle;               // Whether to use RLE compression
        } bmp;
        
        struct {
            bool binary;            // Whether to use binary (P6) or ASCII (P3) encoding
        } ppm;
        
        struct {
            int compressionType;    // 0=none, 1=RLE, 2=LZW, etc.
        } tiff;
        
        struct {
            int compressionType;    // 0=none, 1=RLE, 2=ZIP, 3=PIZ, etc.
            bool halfFloat;         // Whether to use 16-bit half float
        } exr;
    } specific;
    
    // Constructor with defaults
    FormatOptions() : includeMetadata(true), progressiveLoading(false) {
        memset(&specific, 0, sizeof(specific));
        specific.png.compressionLevel = 6;
        specific.png.interlaced = false;
        specific.jpg.quality = 85;
        specific.jpg.progressive = false;
        specific.bmp.rle = false;
        specific.ppm.binary = true;
        specific.tiff.compressionType = 1;  // RLE
        specific.exr.compressionType = 3;   // PIZ
        specific.exr.halfFloat = true;
    }
};

/**
 * @brief Write BMP file header
 * 
 * @param file Output stream
 * @param width Image width
 * @param height Image height
 * @param paddedRowSize Row size including padding
 */
void writeBMPHeader(std::ofstream& file, int width, int height, int paddedRowSize) {
    int fileSize = 54 + (paddedRowSize * height); // 54 bytes for header
    
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
    dibHeader[4] = width & 0xFF;
    dibHeader[5] = (width >> 8) & 0xFF;
    dibHeader[6] = (width >> 16) & 0xFF;
    dibHeader[7] = (width >> 24) & 0xFF;
    
    dibHeader[8] = height & 0xFF;
    dibHeader[9] = (height >> 8) & 0xFF;
    dibHeader[10] = (height >> 16) & 0xFF;
    dibHeader[11] = (height >> 24) & 0xFF;
    
    // Write headers
    file.write(reinterpret_cast<char*>(bmpHeader), sizeof(bmpHeader));
    file.write(reinterpret_cast<char*>(dibHeader), sizeof(dibHeader));
}

/**
 * @brief Save fractal data as BMP file
 * 
 * @param filename Output filename
 * @param pixelData RGB pixel data (row-major order)
 * @param width Image width
 * @param height Image height
 * @param options Format options
 * @return bool Success indicator
 */
bool saveBMP(
    const std::string& filename,
    const std::vector<std::vector<Color>>& pixelData,
    int width,
    int height,
    const FormatOptions& options
) {
    try {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        
        // Calculate row size and padding
        int rowSize = width * 3; // 3 bytes per pixel (RGB)
        int padding = (4 - (rowSize % 4)) % 4; // Row size must be multiple of 4
        int paddedRowSize = rowSize + padding;
        
        // Write BMP header
        writeBMPHeader(file, width, height, paddedRowSize);
        
        // Write pixel data (BMP stores rows from bottom to top)
        std::vector<unsigned char> paddingBytes(3, 0);
        
        for (int y = height - 1; y >= 0; y--) {
            for (int x = 0; x < width; x++) {
                const Color& color = pixelData[y][x];
                
                // BMP stores colors as BGR
                unsigned char pixel[3] = {
                    color.b,
                    color.g,
                    color.r
                };
                
                file.write(reinterpret_cast<char*>(pixel), 3);
            }
            
            // Write padding bytes
            if (padding > 0) {
                file.write(reinterpret_cast<char*>(paddingBytes.data()), padding);
            }
        }
        
        std::cout << "BMP image saved successfully to " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error saving BMP file: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief Save fractal data as PPM file
 * 
 * @param filename Output filename
 * @param pixelData RGB pixel data (row-major order)
 * @param width Image width
 * @param height Image height
 * @param options Format options
 * @return bool Success indicator
 */
bool savePPM(
    const std::string& filename,
    const std::vector<std::vector<Color>>& pixelData,
    int width,
    int height,
    const FormatOptions& options
) {
    try {
        std::ofstream file(filename, options.specific.ppm.binary ? 
                          std::ios::binary : std::ios::out);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        
        // Write PPM header
        if (options.specific.ppm.binary) {
            file << "P6\n"; // Binary RGB format
        } else {
            file << "P3\n"; // ASCII RGB format
        }
        
        // Add metadata as comments if requested
        if (options.includeMetadata) {
            file << "# Generated by AlmndBrd-Mandelbrot\n";
            file << "# Width: " << width << ", Height: " << height << "\n";
        }
        
        file << width << " " << height << "\n";
        file << "255\n"; // Max color value
        
        // Write pixel data
        if (options.specific.ppm.binary) {
            // Binary format (more compact)
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    const Color& color = pixelData[y][x];
                    unsigned char pixel[3] = {color.r, color.g, color.b};
                    file.write(reinterpret_cast<char*>(pixel), 3);
                }
            }
        } else {
            // ASCII format (human-readable)
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    const Color& color = pixelData[y][x];
                    file << static_cast<int>(color.r) << " "
                         << static_cast<int>(color.g) << " "
                         << static_cast<int>(color.b) << " ";
                }
                file << "\n";
            }
        }
        
        std::cout << "PPM image saved successfully to " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error saving PPM file: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief Get file extension for the specified format
 * 
 * @param format Image format type
 * @return std::string File extension including dot
 */
std::string getFileExtension(ImageFormatType format) {
    switch (format) {
        case ImageFormatType::BMP: return ".bmp";
        case ImageFormatType::PNG: return ".png";
        case ImageFormatType::JPG: return ".jpg";
        case ImageFormatType::PPM: return ".ppm";
        case ImageFormatType::TIFF: return ".tiff";
        case ImageFormatType::EXR: return ".exr";
        default: return ".bmp";
    }
}

/**
 * @brief Save fractal visualization to the specified image format
 * 
 * @param filename Output filename (extension will be added if not present)
 * @param iterationData 2D array of iteration counts
 * @param width Image width
 * @param height Image height
 * @param maxIterations Maximum iteration value
 * @param colorMode Color mapping mode
 * @param format Target image format
 * @param options Format-specific options
 * @return bool Success indicator
 */
bool saveFractalImage(
    const std::string& filename,
    const std::vector<std::vector<int>>& iterationData,
    int width,
    int height,
    int maxIterations,
    ColorMode colorMode = ColorMode::LOGARITHMIC,
    ImageFormatType format = ImageFormatType::BMP,
    const FormatOptions& options = FormatOptions()
) {
    // Check dimensions
    if (iterationData.size() != height || 
        (height > 0 && iterationData[0].size() != width)) {
        std::cerr << "Error: Iteration data dimensions do not match specified width and height" << std::endl;
        return false;
    }
    
    // Generate color palette
    std::vector<Color> palette = generateColorPalette(colorMode, maxIterations);
    
    // Create pixel data array
    std::vector<std::vector<Color>> pixelData(height, std::vector<Color>(width));
    
    // Map iteration values to colors
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iterations = iterationData[y][x];
            pixelData[y][x] = palette[iterations];
        }
    }
    
    // Ensure filename has correct extension
    std::string outputFilename = filename;
    std::string extension = getFileExtension(format);
    
    if (outputFilename.size() < extension.size() || 
        outputFilename.substr(outputFilename.size() - extension.size()) != extension) {
        outputFilename += extension;
    }
    
    // Save in the requested format
    switch (format) {
        case ImageFormatType::BMP:
            return saveBMP(outputFilename, pixelData, width, height, options);
            
        case ImageFormatType::PPM:
            return savePPM(outputFilename, pixelData, width, height, options);
            
        case ImageFormatType::PNG:
            std::cout << "PNG format requires external libraries (e.g., libpng)." << std::endl;
            std::cout << "Implementation notes are provided below." << std::endl;
            return false;
            
        case ImageFormatType::JPG:
            std::cout << "JPEG format requires external libraries (e.g., libjpeg)." << std::endl;
            std::cout << "Implementation notes are provided below." << std::endl;
            return false;
            
        case ImageFormatType::TIFF:
            std::cout << "TIFF format requires external libraries (e.g., libtiff)." << std::endl;
            std::cout << "Implementation notes are provided below." << std::endl;
            return false;
            
        case ImageFormatType::EXR:
            std::cout << "EXR format requires external libraries (e.g., OpenEXR)." << std::endl;
            std::cout << "Implementation notes are provided below." << std::endl;
            return false;
            
        default:
            std::cerr << "Unsupported format requested." << std::endl;
            return false;
    }
}

/**
 * @brief Information about PNG implementation (not actual implementation)
 * 
 * This function provides notes about how PNG support would be implemented
 * using the libpng library. It's included for educational purposes.
 */
void pngImplementationNotes() {
    std::cout << "=== PNG Implementation Notes ===" << std::endl;
    std::cout << "PNG support can be added using libpng:" << std::endl;
    std::cout << std::endl;
    
    std::cout << "1. Dependencies:" << std::endl;
    std::cout << "   - libpng (http://www.libpng.org/)" << std::endl;
    std::cout << "   - zlib (for compression)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "2. Implementation steps:" << std::endl;
    std::cout << "   a. Initialize PNG structures (png_create_write_struct, png_create_info_struct)" << std::endl;
    std::cout << "   b. Set up error handling with setjmp" << std::endl;
    std::cout << "   c. Initialize I/O (png_init_io)" << std::endl;
    std::cout << "   d. Write header (png_set_IHDR)" << std::endl;
    std::cout << "   e. Set compression options if needed" << std::endl;
    std::cout << "   f. Write image data row by row" << std::endl;
    std::cout << "   g. Finalize (png_write_end)" << std::endl;
    std::cout << "   h. Clean up resources" << std::endl;
    std::cout << std::endl;
    
    std::cout << "3. Advantages of PNG for fractals:" << std::endl;
    std::cout << "   - Lossless compression (important for sharp fractal boundaries)" << std::endl;
    std::cout << "   - Alpha channel support (for transparent overlays)" << std::endl;
    std::cout << "   - Metadata support (for storing fractal parameters)" << std::endl;
    std::cout << "   - Color depth options (8-bit through 16-bit per channel)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "For full implementation, refer to libpng documentation and examples." << std::endl;
}

/**
 * @brief Convert iteration count array to colored pixel data
 * 
 * @param iterationData 2D array of iteration counts
 * @param maxIterations Maximum iterations
 * @param colorMode Color mapping mode
 * @return std::vector<std::vector<Color>> Colored pixel data
 */
std::vector<std::vector<Color>> convertIterationDataToPixels(
    const std::vector<std::vector<int>>& iterationData,
    int maxIterations,
    ColorMode colorMode
) {
    int height = iterationData.size();
    int width = (height > 0) ? iterationData[0].size() : 0;
    
    // Generate color palette
    std::vector<Color> palette = generateColorPalette(colorMode, maxIterations);
    
    // Create pixel data array
    std::vector<std::vector<Color>> pixelData(height, std::vector<Color>(width));
    
    // Map iteration values to colors
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int iterations = iterationData[y][x];
            if (iterations >= 0 && iterations <= maxIterations) {
                pixelData[y][x] = palette[iterations];
            } else {
                // Default color for invalid iteration values
                pixelData[y][x] = Color(0, 0, 0);
            }
        }
    }
    
    return pixelData;
}

/**
 * @brief Main function to test image format functionality
 */
int main() {
    // Parameters for test image
    const int WIDTH = 800;
    const int HEIGHT = 600;
    const int MAX_ITERATIONS = 100;
    
    // Create a simple test pattern (a gradient based on position)
    std::vector<std::vector<int>> testPattern(HEIGHT, std::vector<int>(WIDTH));
    
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // Create a radial pattern
            double dx = x - WIDTH/2;
            double dy = y - HEIGHT/2;
            double distance = std::sqrt(dx*dx + dy*dy);
            
            // Map distance to iteration count (for testing)
            int iterations = static_cast<int>(MAX_ITERATIONS * distance / (WIDTH/2));
            iterations = std::min(MAX_ITERATIONS, std::max(0, iterations));
            
            testPattern[y][x] = iterations;
        }
    }
    
    // Save in different formats for testing
    FormatOptions options;
    
    // BMP format
    saveFractalImage(
        "test_pattern_linear",
        testPattern,
        WIDTH, HEIGHT,
        MAX_ITERATIONS,
        ColorMode::LINEAR,
        ImageFormatType::BMP,
        options
    );
    
    // PPM format
    options.specific.ppm.binary = false; // ASCII format for inspection
    saveFractalImage(
        "test_pattern_logarithmic",
        testPattern,
        WIDTH, HEIGHT,
        MAX_ITERATIONS,
        ColorMode::LOGARITHMIC,
        ImageFormatType::PPM,
        options
    );
    
    // Notes about PNG implementation
    pngImplementationNotes();
    
    return 0;
}