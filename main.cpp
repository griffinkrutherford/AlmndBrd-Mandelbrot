#include <SDL.h>
#include <SDL_ttf.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <random>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>
#include <chrono>

#include "Complex Number Handling (z, c).h"

enum class MLModelType {
    KNN, RANDOM_FOREST, NEURAL_NET, HYBRID, EVOLUTIONARY
};

extern std::vector<std::pair<Complex, double>> predictInterestingRegions(
    MLModelType modelType,
    int trainingSize,
    double regionSize,
    int resolution
);

// Global variables to store screen dimensions, will be set during initialization
int SCREEN_WIDTH = 1920;  // Default values, will be updated based on display
int SCREEN_HEIGHT = 1080; // Default values, will be updated based on display
const Uint32 POPUP_DURATION = 2000;

// Mutex for renderer access
std::mutex rendererMutex;

// Number of threads to use for parallel computation
const int NUM_THREADS = std::thread::hardware_concurrency() > 0 
                      ? std::thread::hardware_concurrency() 
                      : 4;  // Default to 4 if we can't detect

// Dynamic iteration scaling
int getMaxIterationsForZoom(double zoomLevel) {
    // Base iterations
    const int BASE_ITERATIONS = 500;
    // Maximum iterations cap (to prevent excessive computation)
    const int MAX_ITERATIONS_CAP = 10000;
    
    // Scale iterations logarithmically with zoom level
    int iterations = static_cast<int>(BASE_ITERATIONS * std::log10(1.0 + zoomLevel));
    
    // Apply a minimum and maximum cap
    return std::max(BASE_ITERATIONS, std::min(iterations, MAX_ITERATIONS_CAP));
}

struct Viewport {
    double x_min = -2.0, x_max = 1.0;
    double y_min = -1.5, y_max = 1.5;
    
    void zoom(double cx, double cy, double factor) {
        double width = x_max - x_min;
        double height = y_max - y_min;
        double new_width = width * factor;
        double new_height = height * factor;
        double dx = (cx - x_min) / width;
        double dy = (cy - y_min) / height;
        x_min = cx - new_width * dx;
        x_max = cx + new_width * (1 - dx);
        y_min = cy - new_height * dy;
        y_max = cy + new_height * (1 - dy);
    }
    
    // Zoom to a specific rectangle with aspect ratio preservation
    void zoomToRect(int x1, int y1, int x2, int y2) {
        // Ensure we have a valid rectangle
        if (x1 > x2) std::swap(x1, x2);
        if (y1 > y2) std::swap(y1, y2);
        
        // Calculate dimensions
        int width = x2 - x1;
        int height = y2 - y1;
        
        // Calculate center of the selection in screen coordinates
        int centerX = x1 + width / 2;
        int centerY = y1 + height / 2;
        
        // Convert center to complex plane coordinates
        double centerReal = x_min + (centerX / (double)SCREEN_WIDTH) * (x_max - x_min);
        double centerImag = y_min + (centerY / (double)SCREEN_HEIGHT) * (y_max - y_min);
        
        // Calculate new dimensions based on the larger of width/height ratio
        double currentAspectRatio = (x_max - x_min) / (y_max - y_min);
        double selectionRatio = (double)width / height;
        
        double newWidth, newHeight;
        if (selectionRatio > currentAspectRatio) {
            // Selection is wider than current view
            newWidth = (x_max - x_min) * ((double)width / SCREEN_WIDTH);
            newHeight = newWidth / currentAspectRatio;
        } else {
            // Selection is taller than current view
            newHeight = (y_max - y_min) * ((double)height / SCREEN_HEIGHT);
            newWidth = newHeight * currentAspectRatio;
        }
        
        // Set new viewport bounds centered on the selection
        x_min = centerReal - newWidth / 2.0;
        x_max = centerReal + newWidth / 2.0;
        y_min = centerImag - newHeight / 2.0;
        y_max = centerImag + newHeight / 2.0;
    }
    
    void pan(double dx, double dy) {
        x_min += dx; x_max += dx;
        y_min += dy; y_max += dy;
    }
    
    double getZoomLevel() const {
        return 3.0 / (x_max - x_min);
    }
};

struct Color {
    Uint8 r, g, b;
};

Color hsv2rgb(double h, double s, double v) {
    double c = v * s;
    double h_prime = fmod(h / 60.0, 6);
    double x = c * (1 - fabs(fmod(h_prime, 2) - 1));
    double m = v - c;
    double r_, g_, b_;
    if (0 <= h_prime && h_prime < 1) { r_ = c; g_ = x; b_ = 0; }
    else if (1 <= h_prime && h_prime < 2) { r_ = x; g_ = c; b_ = 0; }
    else if (2 <= h_prime && h_prime < 3) { r_ = 0; g_ = c; b_ = x; }
    else if (3 <= h_prime && h_prime < 4) { r_ = 0; g_ = x; b_ = c; }
    else if (4 <= h_prime && h_prime < 5) { r_ = x; g_ = 0; b_ = c; }
    else { r_ = c; g_ = 0; b_ = x; }
    return { static_cast<Uint8>((r_ + m) * 255),
             static_cast<Uint8>((g_ + m) * 255),
             static_cast<Uint8>((b_ + m) * 255) };
}

double mandelbrotIterations(Complex c, int max_iter) {
    Complex z(0, 0);
    int iter = 0;
    while (z.magnitudeSquared() < 4.0 && iter < max_iter) {
        z = z * z + c;
        iter++;
    }
    if (iter == max_iter)
        return max_iter;
    double log_zn = std::log(z.magnitudeSquared()) / 2.0;
    double nu = std::log(log_zn / std::log(2)) / std::log(2);
    return iter + 1 - nu;
}

// Function to get dynamic coloring based on zoom level
Color getColorForIteration(double iter, int max_iter, double zoomLevel) {
    if (iter >= max_iter)
        return {0, 0, 0};  // Black for points in the set
    
    // Create different color schemes based on zoom level
    // This cycles through several color schemes as you zoom deeper
    int colorScheme = static_cast<int>(std::log10(zoomLevel)) % 5;
    
    double normalized = iter / max_iter;
    
    switch (colorScheme) {
        case 0: {
            // Rainbow coloring (default)
            double hue = 360.0 * normalized;
            return hsv2rgb(hue, 1.0, 1.0);
        }
        case 1: {
            // Blue-to-yellow gradient
            double hue = 240.0 + normalized * 60.0;
            return hsv2rgb(hue, 0.8, 1.0);
        }
        case 2: {
            // Purple-to-green gradient
            double hue = 280.0 + normalized * 100.0;
            return hsv2rgb(hue, 0.7, 0.9);
        }
        case 3: {
            // Red-to-blue-to-green cycling
            double hue = fmod(720.0 * normalized, 360.0);
            return hsv2rgb(hue, 0.9, 0.8);
        }
        case 4: {
            // Fire colors (red/yellow/orange)
            double hue = 0.0 + normalized * 60.0;
            return hsv2rgb(hue, 1.0, 1.0);
        }
        default:
            // Fallback to grayscale
            int gray = static_cast<int>(255.0 * normalized);
            return {static_cast<Uint8>(gray), static_cast<Uint8>(gray), static_cast<Uint8>(gray)};
    }
}

// Struct to store parameters for worker threads
struct ThreadParams {
    int yStart, yEnd;  // Range of rows to process
    int width, height;  // Image dimensions
    double x_min, x_max, y_min, y_max;  // Viewport
    int max_iter;  // Maximum iterations
    double zoomLevel;  // Current zoom level
    std::vector<std::vector<Color>>* resultBuffer;  // Output buffer
    std::atomic<int>* completedLines;  // Progress counter
};

// Worker thread function to compute a portion of the Mandelbrot set
void renderMandelbrotThread(ThreadParams params) {
    // Process assigned rows
    for (int y = params.yStart; y < params.yEnd; y++) {
        for (int x = 0; x < params.width; x++) {
            // Map pixel to complex plane
            double real = params.x_min + (x / (double)params.width) * (params.x_max - params.x_min);
            double imag = params.y_min + (y / (double)params.height) * (params.y_max - params.y_min);
            Complex c(real, imag);
            
            // Calculate iterations
            double iter = mandelbrotIterations(c, params.max_iter);
            
            // Get color based on iteration count and zoom level
            Color color = getColorForIteration(iter, params.max_iter, params.zoomLevel);
            
            // Store the result
            (*params.resultBuffer)[y][x] = color;
        }
        
        // Update progress
        if (params.completedLines) {
            (*params.completedLines)++;
        }
    }
}

// Multithreaded rendering of the Mandelbrot set to a texture
SDL_Texture* renderMandelbrotToTexture(SDL_Renderer* renderer, Viewport& viewport, int max_iter, std::atomic<int>* progressCounter = nullptr) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Create a texture to render to - using mutex for renderer access
    SDL_Texture* texture = nullptr;
    {
        std::lock_guard<std::mutex> lock(rendererMutex);
        texture = SDL_CreateTexture(
            renderer, 
            SDL_PIXELFORMAT_RGBA8888, 
            SDL_TEXTUREACCESS_STREAMING, 
            SCREEN_WIDTH, SCREEN_HEIGHT
        );
    }
    
    if (!texture) {
        std::cerr << "Failed to create texture: " << SDL_GetError() << "\n";
        return nullptr;
    }
    
    // Create color buffer
    std::vector<std::vector<Color>> colorBuffer(SCREEN_HEIGHT, std::vector<Color>(SCREEN_WIDTH));
    
    // Reset progress counter
    if (progressCounter) {
        *progressCounter = 0;
    }
    
    // Get current zoom level for dynamic coloring
    double zoomLevel = viewport.getZoomLevel();
    
    // Calculate rows per thread
    int rowsPerThread = SCREEN_HEIGHT / NUM_THREADS;
    std::vector<std::thread> threads;
    
    // Create and start worker threads
    for (int i = 0; i < NUM_THREADS; i++) {
        int startY = i * rowsPerThread;
        int endY = (i == NUM_THREADS - 1) ? SCREEN_HEIGHT : (i + 1) * rowsPerThread;
        
        ThreadParams params = {
            startY, endY,
            SCREEN_WIDTH, SCREEN_HEIGHT,
            viewport.x_min, viewport.x_max, viewport.y_min, viewport.y_max,
            max_iter,
            zoomLevel,
            &colorBuffer,
            progressCounter
        };
        
        threads.push_back(std::thread(renderMandelbrotThread, params));
    }
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    
    // Update the texture with the color buffer - using mutex for renderer access
    void* pixels;
    int pitch;
    
    {
        std::lock_guard<std::mutex> lock(rendererMutex);
        if (SDL_LockTexture(texture, nullptr, &pixels, &pitch) < 0) {
            std::cerr << "Failed to lock texture: " << SDL_GetError() << "\n";
            SDL_DestroyTexture(texture);
            return nullptr;
        }
        
        // Copy color data to texture
        for (int y = 0; y < SCREEN_HEIGHT; y++) {
            for (int x = 0; x < SCREEN_WIDTH; x++) {
                Uint32* pixel = (Uint32*)((Uint8*)pixels + y * pitch + x * 4);
                Color color = colorBuffer[y][x];
                *pixel = SDL_MapRGBA(SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888), color.r, color.g, color.b, 255);
            }
        }
        
        SDL_UnlockTexture(texture);
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Render time: " << duration.count() << " ms using " << NUM_THREADS << " threads" << std::endl;
    
    return texture;
}

void renderPopup(SDL_Renderer* renderer, TTF_Font* font, const std::string& message) {
    std::lock_guard<std::mutex> lock(rendererMutex);
    SDL_Color textColor = {255, 255, 255, 255};
    SDL_Surface* textSurface = TTF_RenderText_Blended(font, message.c_str(), textColor);
    if (!textSurface) {
        std::cerr << "Text Render Error: " << TTF_GetError() << "\n";
        return;
    }
    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
    int textW = 0, textH = 0;
    SDL_QueryTexture(textTexture, NULL, NULL, &textW, &textH);
    
    SDL_Rect bgRect = {10, 10, textW + 20, textH + 20};
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 150);
    SDL_RenderFillRect(renderer, &bgRect);
    
    SDL_Rect textRect = {20, 20, textW, textH};
    SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
    
    SDL_DestroyTexture(textTexture);
    SDL_FreeSurface(textSurface);
}

// Render progress bar to show rendering status
void renderProgressBar(SDL_Renderer* renderer, int progress, int total, int width, int height, int x, int y) {
    std::lock_guard<std::mutex> lock(rendererMutex);
    // Background
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 200);
    SDL_Rect bgRect = {x, y, width, height};
    SDL_RenderFillRect(renderer, &bgRect);
    
    // Border
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderDrawRect(renderer, &bgRect);
    
    // Progress
    if (total > 0 && progress > 0) {
        int fillWidth = (progress * (width - 4)) / total;
        SDL_Rect fillRect = {x + 2, y + 2, fillWidth, height - 4};
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
        SDL_RenderFillRect(renderer, &fillRect);
    }
}

// Draw selection rectangle with locked aspect ratio
void renderSelectionRect(SDL_Renderer* renderer, int startX, int startY, int endX, int endY) {
    std::lock_guard<std::mutex> lock(rendererMutex);
    // Calculate differences
    int width = abs(endX - startX);
    int height = abs(endY - startY);
    
    // Lock aspect ratio to match the window's aspect ratio
    double aspectRatio = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
    
    // Determine which dimension to adjust based on drag direction
    if (width > height * aspectRatio) {
        // Width is proportionally larger, adjust height
        height = width / aspectRatio;
    } else {
        // Height is proportionally larger, adjust width
        width = height * aspectRatio;
    }
    
    // Calculate final coordinates based on drag direction
    int x1, y1, x2, y2;
    if (endX >= startX) {
        x1 = startX;
        x2 = startX + width;
    } else {
        x1 = startX - width;
        x2 = startX;
    }
    
    if (endY >= startY) {
        y1 = startY;
        y2 = startY + height;
    } else {
        y1 = startY - height;
        y2 = startY;
    }
    
    // Clamp to screen boundaries
    x1 = std::max(0, std::min(SCREEN_WIDTH-1, x1));
    x2 = std::max(0, std::min(SCREEN_WIDTH-1, x2));
    y1 = std::max(0, std::min(SCREEN_HEIGHT-1, y1));
    y2 = std::max(0, std::min(SCREEN_HEIGHT-1, y2));
    
    // Draw semi-transparent fill
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(renderer, 200, 200, 255, 50);
    SDL_Rect fillRect = {std::min(x1, x2), std::min(y1, y2), 
                         abs(x2 - x1), abs(y2 - y1)};
    SDL_RenderFillRect(renderer, &fillRect);
    
    // Draw outline
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 200);
    SDL_Rect outlineRect = {std::min(x1, x2), std::min(y1, y2), 
                           abs(x2 - x1), abs(y2 - y1)};
    SDL_RenderDrawRect(renderer, &outlineRect);
}

// Function to get the highest available resolution
bool getHighestResolution(int& outWidth, int& outHeight) {
    // Initialize return values
    outWidth = 1920;  // Default fallback values
    outHeight = 1080;
    
    // Get the number of display modes
    int displayIndex = 0; // Primary display
    int numDisplayModes = SDL_GetNumDisplayModes(displayIndex);
    
    if (numDisplayModes < 1) {
        std::cerr << "No display modes available: " << SDL_GetError() << "\n";
        return false;
    }
    
    // Start with the first mode
    SDL_DisplayMode highestMode;
    if (SDL_GetDisplayMode(displayIndex, 0, &highestMode) != 0) {
        std::cerr << "Failed to get display mode: " << SDL_GetError() << "\n";
        return false;
    }
    
    // Find the mode with the highest resolution
    for (int i = 0; i < numDisplayModes; i++) {
        SDL_DisplayMode mode;
        if (SDL_GetDisplayMode(displayIndex, i, &mode) != 0) {
            continue;
        }
        
        // Check if this mode has higher resolution
        if (mode.w * mode.h > highestMode.w * highestMode.h) {
            highestMode = mode;
        }
    }
    
    // Set output values
    outWidth = highestMode.w;
    outHeight = highestMode.h;
    
    std::cout << "Highest resolution found: " << outWidth << "x" << outHeight << "\n";
    return true;
}

// GenAI integration: Navigate to an interesting region using ML prediction
void navigateToInterestingRegion(Viewport& viewport, MLModelType modelType = MLModelType::KNN) {
    try {
        // Get interesting regions from ML module
        auto suggestions = predictInterestingRegions(
            modelType,  // Model type
            1000,       // Training size
            4.0,        // Region size
            100         // Resolution
        );
        
        if (!suggestions.empty()) {
            // Sort by confidence (higher confidence first)
            std::sort(suggestions.begin(), suggestions.end(), 
                     [](const auto& a, const auto& b) {
                         return a.second > b.second;
                     });
            
            // Get the most interesting point
            Complex point = suggestions[0].first;
            double confidence = suggestions[0].second;
            
            // Navigate to the interesting region
            // Calculate region size based on confidence (higher confidence = closer zoom)
            double regionSize = 3.0 * (1.0 - confidence * 0.5);  // Scale between 1.5 and 3.0
            
            // Update viewport
            double centerX = point.getReal();
            double centerY = point.getImag();
            double width = regionSize;
            double height = regionSize * SCREEN_HEIGHT / SCREEN_WIDTH;
            
            viewport.x_min = centerX - width/2;
            viewport.x_max = centerX + width/2;
            viewport.y_min = centerY - height/2;
            viewport.y_max = centerY + height/2;
            
            std::cout << "Navigated to interesting region at (" << centerX << ", " 
                      << centerY << ") with confidence " << confidence 
                      << " and zoom level " << viewport.getZoomLevel() << std::endl;
        } else {
            std::cout << "No interesting regions found." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error using GenAI to find interesting regions: " << e.what() << std::endl;
    }
}

int main(int, char*[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL_Init failed: " << SDL_GetError() << "\n";
        return 1;
    }
    if (TTF_Init() == -1) {
        std::cerr << "TTF_Init failed: " << TTF_GetError() << "\n";
        SDL_Quit();
        return 1;
    }
    
    // Get the highest resolution available
    if (!getHighestResolution(SCREEN_WIDTH, SCREEN_HEIGHT)) {
        std::cerr << "Warning: Using default resolution of " << SCREEN_WIDTH << "x" << SCREEN_HEIGHT << "\n";
    }
    
    std::cout << "Using " << NUM_THREADS << " threads for Mandelbrot calculation" << std::endl;
    
    // Create a fullscreen window by default
    SDL_Window* window = SDL_CreateWindow("Mandelbrot Explorer",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 
        SCREEN_WIDTH, SCREEN_HEIGHT, 
        SDL_WINDOW_FULLSCREEN_DESKTOP);  // Start in fullscreen by default
    
    if (!window) {
        std::cerr << "SDL_CreateWindow failed: " << SDL_GetError() << "\n";
        TTF_Quit();
        SDL_Quit();
        return 1;
    }
    
    // Get the actual window size (in case we're using SDL_WINDOW_FULLSCREEN_DESKTOP)
    SDL_GetWindowSize(window, &SCREEN_WIDTH, &SCREEN_HEIGHT);
    std::cout << "Actual window size: " << SCREEN_WIDTH << "x" << SCREEN_HEIGHT << "\n";
    
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        std::cerr << "SDL_CreateRenderer failed: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        TTF_Quit();
        SDL_Quit();
        return 1;
    }
    
    TTF_Font* font = nullptr;
    const char* fontPaths[] = {
        "/System/Library/Fonts/Times.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",  // Added Linux font path
        "/Library/Fonts/Times New Roman.ttf",
        nullptr
    };
    
    for (int i = 0; fontPaths[i] != nullptr && !font; ++i) {
        font = TTF_OpenFont(fontPaths[i], 24);
        if (font) {
            std::cout << "Loaded font: " << fontPaths[i] << "\n";
        } else {
            std::cerr << "Failed to load " << fontPaths[i] << ": " << TTF_GetError() << "\n";
        }
    }
    
    if (!font) {
        std::cerr << "All font loading attempts failed.\n";
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        TTF_Quit();
        SDL_Quit();
        return 1;
    }
    
    // Adjust viewport for the aspect ratio of the screen
    Viewport viewport;
    double screenAspectRatio = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
    double fractalAspectRatio = (viewport.x_max - viewport.x_min) / (viewport.y_max - viewport.y_min);
    
    // Adjust viewport to match screen aspect ratio
    if (screenAspectRatio > fractalAspectRatio) {
        // Screen is wider than fractal, adjust fractal width
        double fractalHeight = viewport.y_max - viewport.y_min;
        double newWidth = fractalHeight * screenAspectRatio;
        double widthDiff = newWidth - (viewport.x_max - viewport.x_min);
        viewport.x_min -= widthDiff / 2;
        viewport.x_max += widthDiff / 2;
    } else if (screenAspectRatio < fractalAspectRatio) {
        // Screen is taller than fractal, adjust fractal height
        double fractalWidth = viewport.x_max - viewport.x_min;
        double newHeight = fractalWidth / screenAspectRatio;
        double heightDiff = newHeight - (viewport.y_max - viewport.y_min);
        viewport.y_min -= heightDiff / 2;
        viewport.y_max += heightDiff / 2;
    }
    
    bool running = true;
    SDL_Event event;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    bool showPopup = false;
    Uint32 popupStartTime = 0;
    std::string popupMessage = "";
    
    // Selection rectangle variables
    bool isSelecting = false;
    int selectionStartX = 0, selectionStartY = 0;
    int selectionEndX = 0, selectionEndY = 0;
    
    // Smooth zoom variables
    float zoomVelocity = 0.0f;  // Accumulated scroll velocity
    const float ZOOM_SENSITIVITY = 0.05f;  // Scroll input multiplier
    const float ZOOM_DAMPING = 0.9f;  // Smoothing factor (0.0 = instant, 1.0 = no change)
    const float MIN_ZOOM_FACTOR = 0.99f;  // Minimum frame-by-frame zoom in
    const float MAX_ZOOM_FACTOR = 1.01f;  // Maximum frame-by-frame zoom out
    
    // Progress tracking for rendering
    std::atomic<int> renderProgress(0);
    
    // Show initial loading message
    {
        std::lock_guard<std::mutex> lock(rendererMutex);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
    }
    renderPopup(renderer, font, "Initializing Mandelbrot Explorer with " + std::to_string(NUM_THREADS) + " threads...");
    {
        std::lock_guard<std::mutex> lock(rendererMutex);
        SDL_RenderPresent(renderer);
    }
    
    // Pre-render initial Mandelbrot set
    SDL_Texture* mandelbrotTexture = nullptr;
    
    // Get dynamic iterations based on zoom level
    int maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
    
    // Use a future to render in a separate thread while showing progress
    std::future<SDL_Texture*> renderFuture = std::async(std::launch::async, [&]() {
        return renderMandelbrotToTexture(renderer, viewport, maxIterations, &renderProgress);
    });
    
    // Show progress while rendering
    while (renderFuture.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready) {
        {
            std::lock_guard<std::mutex> lock(rendererMutex);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);
        }
        
        // Show rendering progress
        renderProgressBar(renderer, renderProgress, SCREEN_HEIGHT, 
                         SCREEN_WIDTH / 2, 30, 
                         SCREEN_WIDTH / 4, SCREEN_HEIGHT / 2 - 15);
        
        // Display status message
        std::string statusMsg = "Rendering Mandelbrot Set... " + 
                              std::to_string(renderProgress * 100 / SCREEN_HEIGHT) + "%";
        renderPopup(renderer, font, statusMsg);
        
        {
            std::lock_guard<std::mutex> lock(rendererMutex);
            SDL_RenderPresent(renderer);
        }
        
        SDL_Delay(10);
    }
    
    // Get the rendered texture
    mandelbrotTexture = renderFuture.get();
    
    // Show initial welcome message
    popupMessage = "Multithreaded Mandelbrot Explorer at " + std::to_string(SCREEN_WIDTH) + "x" + std::to_string(SCREEN_HEIGHT) + 
                  " using " + std::to_string(NUM_THREADS) + " threads";
    popupStartTime = SDL_GetTicks();
    showPopup = true;
    
    // Flag to indicate if rendering is in progress
    bool isRendering = false;
    
    // Timestamps for rendering control
    Uint32 lastRenderTime = 0;
    const Uint32 MIN_RENDER_INTERVAL = 200; // Minimum milliseconds between renders
    
    // Enable GenAI exploration mode
    bool genAIExplorationMode = false;
    
    while (running) {
        bool needsRedraw = false;
        
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            } else if (event.type == SDL_MOUSEBUTTONDOWN) {
                if (event.button.button == SDL_BUTTON_LEFT && !isSelecting && !isRendering) {
                    // Start selection rectangle
                    isSelecting = true;
                    selectionStartX = selectionEndX = event.button.x;
                    selectionStartY = selectionEndY = event.button.y;
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                if (isSelecting) {
                    // Update selection rectangle end point
                    selectionEndX = event.motion.x;
                    selectionEndY = event.motion.y;
                }
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (event.button.button == SDL_BUTTON_LEFT && isSelecting && !isRendering) {
                    // Finish selection and zoom to rectangle
                    selectionEndX = event.motion.x;
                    selectionEndY = event.motion.y;
                    
                    // Calculate width and height
                    int width = abs(selectionEndX - selectionStartX);
                    int height = abs(selectionEndY - selectionStartY);
                    
                    // Only zoom if the rectangle has a reasonable size
                    if (width > 10 && height > 10) {
                        // Get coordinates of actual selection rectangle with aspect ratio preservation
                        int x1 = selectionStartX;
                        int y1 = selectionStartY;
                        int x2 = selectionEndX;
                        int y2 = selectionEndY;
                        
                        // Calculate dimensions for aspect ratio lock
                        double aspectRatio = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
                        
                        if (width > height * aspectRatio) {
                            // Width is proportionally larger, adjust height
                            height = width / aspectRatio;
                            
                            // Recalculate y2 based on drag direction
                            if (y2 > y1) {
                                y2 = y1 + height;
                            } else {
                                y2 = y1 - height;
                            }
                        } else {
                            // Height is proportionally larger, adjust width
                            width = height * aspectRatio;
                            
                            // Recalculate x2 based on drag direction
                            if (x2 > x1) {
                                x2 = x1 + width;
                            } else {
                                x2 = x1 - width;
                            }
                        }
                        
                        // Zoom to the aspect-corrected rectangle
                        viewport.zoomToRect(x1, y1, x2, y2);
                        needsRedraw = true;
                        
                        // Update max iterations based on new zoom level
                        maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                        
                        popupMessage = "Zoomed to selection, Zoom Level: " + std::to_string(viewport.getZoomLevel()) + 
                                       ", Iterations: " + std::to_string(maxIterations);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                    }
                    
                    isSelecting = false;
                }
            } else if (event.type == SDL_MOUSEWHEEL && !isRendering) {
                // Accumulate zoom velocity based on scroll direction
                zoomVelocity += event.wheel.y * ZOOM_SENSITIVITY;
                int x, y;
                SDL_GetMouseState(&x, &y);
                double cx = viewport.x_min + (x / (double)SCREEN_WIDTH) * (viewport.x_max - viewport.x_min);
                double cy = viewport.y_min + (y / (double)SCREEN_HEIGHT) * (viewport.y_max - viewport.y_min);
                popupMessage = "Zooming " + std::string(event.wheel.y > 0 ? "in" : "out") + 
                               " at (" + std::to_string(cx) + ", " + std::to_string(cy) + 
                               "), Zoom Level: " + std::to_string(viewport.getZoomLevel());
                popupStartTime = SDL_GetTicks();
                showPopup = true;
            } else if (event.type == SDL_KEYDOWN && !isRendering) {
                bool shiftPressed = (SDL_GetModState() & KMOD_SHIFT) != 0;
                double panFactor = shiftPressed ? 0.01 : 0.1;
                double dx = (viewport.x_max - viewport.x_min) * panFactor;
                double dy = (viewport.y_max - viewport.y_min) * panFactor;
                
                switch (event.key.keysym.sym) {
                    case SDLK_LEFT: {
                        viewport.pan(-dx, 0);
                        needsRedraw = true;
                        popupMessage = "Panned left by " + std::to_string(dx);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_RIGHT: {
                        viewport.pan(dx, 0);
                        needsRedraw = true;
                        popupMessage = "Panned right by " + std::to_string(dx);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_UP: {
                        viewport.pan(0, -dy);
                        needsRedraw = true;
                        popupMessage = "Panned up by " + std::to_string(dy);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_DOWN: {
                        viewport.pan(0, dy);
                        needsRedraw = true;
                        popupMessage = "Panned down by " + std::to_string(dy);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_PLUS:
                    case SDLK_EQUALS: {
                        int x_plus, y_plus;
                        SDL_GetMouseState(&x_plus, &y_plus);
                        double cx_plus = viewport.x_min + (x_plus / (double)SCREEN_WIDTH) * (viewport.x_max - viewport.x_min);
                        double cy_plus = viewport.y_min + (y_plus / (double)SCREEN_HEIGHT) * (viewport.y_max - viewport.y_min);
                        viewport.zoom(cx_plus, cy_plus, 0.9);
                        
                        // Update max iterations based on new zoom level
                        maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                        
                        needsRedraw = true;
                        popupMessage = "Zoomed in to (" + std::to_string(cx_plus) + ", " + 
                                       std::to_string(cy_plus) + "), Zoom Level: " + 
                                       std::to_string(viewport.getZoomLevel()) + 
                                       ", Iterations: " + std::to_string(maxIterations);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_MINUS: {
                        int x_minus, y_minus;
                        SDL_GetMouseState(&x_minus, &y_minus);
                        double cx_minus = viewport.x_min + (x_minus / (double)SCREEN_WIDTH) * (viewport.x_max - viewport.x_min);
                        double cy_minus = viewport.y_min + (y_minus / (double)SCREEN_HEIGHT) * (viewport.y_max - viewport.y_min);
                        viewport.zoom(cx_minus, cy_minus, 1.1);
                        
                        // Update max iterations based on new zoom level
                        maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                        
                        needsRedraw = true;
                        popupMessage = "Zoomed out to (" + std::to_string(cx_minus) + ", " + 
                                       std::to_string(cy_minus) + "), Zoom Level: " + 
                                       std::to_string(viewport.getZoomLevel()) + 
                                       ", Iterations: " + std::to_string(maxIterations);
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_r: {
                        viewport.x_min = -2.0; viewport.x_max = 1.0;
                        viewport.y_min = -1.5; viewport.y_max = 1.5;
                        
                        // Adjust viewport to match screen aspect ratio
                        double screenAspectRatio = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
                        double fractalAspectRatio = (viewport.x_max - viewport.x_min) / (viewport.y_max - viewport.y_min);
                        
                        if (screenAspectRatio > fractalAspectRatio) {
                            // Screen is wider than fractal, adjust fractal width
                            double fractalHeight = viewport.y_max - viewport.y_min;
                            double newWidth = fractalHeight * screenAspectRatio;
                            double widthDiff = newWidth - (viewport.x_max - viewport.x_min);
                            viewport.x_min -= widthDiff / 2;
                            viewport.x_max += widthDiff / 2;
                        } else if (screenAspectRatio < fractalAspectRatio) {
                            // Screen is taller than fractal, adjust fractal height
                            double fractalWidth = viewport.x_max - viewport.x_min;
                            double newHeight = fractalWidth / screenAspectRatio;
                            double heightDiff = newHeight - (viewport.y_max - viewport.y_min);
                            viewport.y_min -= heightDiff / 2;
                            viewport.y_max += heightDiff / 2;
                        }
                        
                        // Reset iterations to default
                        maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                        
                        needsRedraw = true;
                        popupMessage = "Reset view";
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_g: {
                        // Handle GenAI navigation
                        try {
                            popupMessage = "Using GenAI to find interesting regions...";
                            popupStartTime = SDL_GetTicks();
                            showPopup = true;
                            
                            // Render the message before starting computation
                            {
                                std::lock_guard<std::mutex> lock(rendererMutex);
                                SDL_RenderClear(renderer);
                                SDL_RenderCopy(renderer, mandelbrotTexture, NULL, NULL);
                            }
                            renderPopup(renderer, font, popupMessage);
                            {
                                std::lock_guard<std::mutex> lock(rendererMutex);
                                SDL_RenderPresent(renderer);
                            }
                            
                            // Use ML to find interesting regions and navigate there
                            navigateToInterestingRegion(viewport, MLModelType::KNN);
                            
                            // Update max iterations based on new zoom level
                            maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                            
                            popupMessage = "Navigated to AI-suggested region with " + 
                                          std::to_string(maxIterations) + " iterations";
                            popupStartTime = SDL_GetTicks();
                            showPopup = true;
                            needsRedraw = true;
                            
                            // Toggle exploration mode
                            genAIExplorationMode = !genAIExplorationMode;
                            if (genAIExplorationMode) {
                                popupMessage += " - GenAI Exploration Mode ENABLED";
                            }
                        } catch (const std::exception& e) {
                            popupMessage = "GenAI Error: " + std::string(e.what());
                            popupStartTime = SDL_GetTicks();
                            showPopup = true;
                            std::cerr << "GenAI Navigation failed: " << e.what() << "\n";
                        }
                        break;
                    }
                    case SDLK_f: {
                        // Toggle fullscreen
                        Uint32 flags = SDL_GetWindowFlags(window);
                        if (flags & SDL_WINDOW_FULLSCREEN_DESKTOP) {
                            SDL_SetWindowFullscreen(window, 0);
                            popupMessage = "Windowed mode";
                        } else {
                            SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP);
                            popupMessage = "Fullscreen mode";
                        }
                        SDL_GetWindowSize(window, &SCREEN_WIDTH, &SCREEN_HEIGHT);
                        
                        // Adjust the viewport aspect ratio to match the new window size
                        double screenAspectRatio = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
                        double fractalAspectRatio = (viewport.x_max - viewport.x_min) / (viewport.y_max - viewport.y_min);
                        
                        if (screenAspectRatio > fractalAspectRatio) {
                            // Screen is wider than fractal, adjust fractal width
                            double fractalHeight = viewport.y_max - viewport.y_min;
                            double newWidth = fractalHeight * screenAspectRatio;
                            double widthDiff = newWidth - (viewport.x_max - viewport.x_min);
                            viewport.x_min -= widthDiff / 2;
                            viewport.x_max += widthDiff / 2;
                        } else if (screenAspectRatio < fractalAspectRatio) {
                            // Screen is taller than fractal, adjust fractal height
                            double fractalWidth = viewport.x_max - viewport.x_min;
                            double newHeight = fractalWidth / screenAspectRatio;
                            double heightDiff = newHeight - (viewport.y_max - viewport.y_min);
                            viewport.y_min -= heightDiff / 2;
                            viewport.y_max += heightDiff / 2;
                        }
                        
                        needsRedraw = true;
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    case SDLK_ESCAPE: {
                        // If we're currently selecting, cancel the selection
                        if (isSelecting) {
                            isSelecting = false;
                        } else {
                            // Exit fullscreen or quit
                            Uint32 flags = SDL_GetWindowFlags(window);
                            if (flags & SDL_WINDOW_FULLSCREEN_DESKTOP) {
                                SDL_SetWindowFullscreen(window, 0);
                                popupMessage = "Windowed mode";
                                popupStartTime = SDL_GetTicks();
                                showPopup = true;
                            } else {
                                running = false;
                            }
                        }
                        break;
                    }
                    // Handle number keys 1-5 for direct color scheme selection
                    case SDLK_1:
                    case SDLK_2:
                    case SDLK_3:
                    case SDLK_4:
                    case SDLK_5: {
                        int scheme = event.key.keysym.sym - SDLK_1;
                        std::string schemeNames[] = {
                            "Rainbow",
                            "Blue-Yellow",
                            "Purple-Green",
                            "RGB Cycling",
                            "Fire"
                        };
                        
                        // Force redraw with this color scheme
                        needsRedraw = true;
                        popupMessage = "Color scheme: " + schemeNames[scheme];
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    // Add thread control keys
                    case SDLK_t: {
                        if (shiftPressed && NUM_THREADS < 64) {
                            // Increase thread count
                            int newThreads = NUM_THREADS + 1;
                            const_cast<int&>(NUM_THREADS) = newThreads;
                            popupMessage = "Increased threads to " + std::to_string(NUM_THREADS);
                            popupStartTime = SDL_GetTicks();
                            showPopup = true;
                            needsRedraw = true;
                        } else if (!shiftPressed && NUM_THREADS > 1) {
                            // Decrease thread count
                            int newThreads = NUM_THREADS - 1;
                            const_cast<int&>(NUM_THREADS) = newThreads;
                            popupMessage = "Decreased threads to " + std::to_string(NUM_THREADS);
                            popupStartTime = SDL_GetTicks();
                            showPopup = true;
                            needsRedraw = true;
                        }
                        break;
                    }
                    // Add key to toggle GenAI exploration mode
                    case SDLK_a: {
                        genAIExplorationMode = !genAIExplorationMode;
                        popupMessage = genAIExplorationMode ? 
                                      "GenAI Exploration Mode ENABLED" : 
                                      "GenAI Exploration Mode DISABLED";
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                    // Add key to toggle iteration display
                    case SDLK_i: {
                        popupMessage = "Current iterations: " + std::to_string(maxIterations) + 
                                      " at zoom level " + std::to_string(viewport.getZoomLevel());
                        popupStartTime = SDL_GetTicks();
                        showPopup = true;
                        break;
                    }
                }
            }
        }
        
        // Apply smooth zoom if there's velocity and not currently rendering
        if (std::abs(zoomVelocity) > 0.01 && !isRendering) {
            int x, y;
            SDL_GetMouseState(&x, &y);
            double cx = viewport.x_min + (x / (double)SCREEN_WIDTH) * (viewport.x_max - viewport.x_min);
            double cy = viewport.y_min + (y / (double)SCREEN_HEIGHT) * (viewport.y_max - viewport.y_min);
            
            // Calculate zoom factor based on velocity
            float zoomFactor = 1.0f + zoomVelocity;
            zoomFactor = std::max(MIN_ZOOM_FACTOR, std::min(MAX_ZOOM_FACTOR, zoomFactor));  // Clamp factor
            viewport.zoom(cx, cy, zoomFactor);
            
            // Only trigger a redraw after sufficient zooming to avoid continuous re-rendering
            static const float REDRAW_THRESHOLD = 0.03f;
            static float accumulatedZoom = 0.0f;
            accumulatedZoom += std::abs(zoomVelocity);
            
            // Update iterations based on zoom occasionally
            static float iterationUpdateThreshold = 0.1f;
            static float iterAccumulated = 0.0f;
            iterAccumulated += std::abs(zoomVelocity);
            
            if (iterAccumulated > iterationUpdateThreshold) {
                maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                iterAccumulated = 0.0f;
            }
            
            if (accumulatedZoom > REDRAW_THRESHOLD) {
                needsRedraw = true;
                accumulatedZoom = 0.0f;
            }
            
            // Dampen the velocity for smooth decay
            zoomVelocity *= ZOOM_DAMPING;
            
            // Update popup with current zoom info
            popupMessage = "Zooming " + std::string(zoomVelocity > 0 ? "in" : "out") + 
                           " at (" + std::to_string(cx) + ", " + std::to_string(cy) + 
                           "), Zoom Level: " + std::to_string(viewport.getZoomLevel()) +
                           ", Iterations: " + std::to_string(maxIterations);
            popupStartTime = SDL_GetTicks();
            showPopup = true;
        }
        
        // GenAI automatic exploration mode logic
        if (genAIExplorationMode && !isRendering && SDL_GetTicks() - lastRenderTime > 5000) {
            try {
                // Find a new interesting region
                navigateToInterestingRegion(viewport, MLModelType::KNN);
                maxIterations = getMaxIterationsForZoom(viewport.getZoomLevel());
                needsRedraw = true;
                popupMessage = "GenAI exploration: Found new region at zoom level " + 
                              std::to_string(viewport.getZoomLevel());
                popupStartTime = SDL_GetTicks();
                showPopup = true;
            } catch (const std::exception& e) {
                std::cerr << "Error in GenAI exploration: " << e.what() << std::endl;
                // Turn off exploration mode if it keeps failing
                genAIExplorationMode = false;
            }
        }
        
        // Redraw if necessary and not already rendering
        Uint32 currentTime = SDL_GetTicks();
        if (needsRedraw && !isRendering && (currentTime - lastRenderTime > MIN_RENDER_INTERVAL)) {
            isRendering = true;
            renderProgress = 0;
            needsRedraw = false;  // Reset the flag to prevent multiple render triggers
            lastRenderTime = currentTime;  // Update the last render time
            
            // Show initial rendering message
            {
                std::lock_guard<std::mutex> lock(rendererMutex);
                SDL_RenderClear(renderer);
            }
            renderPopup(renderer, font, "Rendering with " + std::to_string(NUM_THREADS) + 
                                        " threads, " + std::to_string(maxIterations) + " iterations...");
            {
                std::lock_guard<std::mutex> lock(rendererMutex);
                SDL_RenderPresent(renderer);
            }
            
            // Clean up old texture if it exists
            if (mandelbrotTexture) {
                std::lock_guard<std::mutex> lock(rendererMutex);
                SDL_DestroyTexture(mandelbrotTexture);
                mandelbrotTexture = nullptr;
            }
            
            // Use async instead of a detached thread for rendering
            renderFuture = std::async(std::launch::async, [&]() {
                return renderMandelbrotToTexture(renderer, viewport, maxIterations, &renderProgress);
            });
            
            std::cout << "Started new render at time: " << currentTime 
                      << " with iterations: " << maxIterations << std::endl;
        }
        
        // Check if rendering is complete
        if (isRendering && renderFuture.valid() && 
            renderFuture.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            try {
                // Get the new texture
                mandelbrotTexture = renderFuture.get();
                isRendering = false;
            } catch (const std::exception& e) {
                std::cerr << "Error during rendering: " << e.what() << std::endl;
                isRendering = false;
            }
        }
        
        // Clear the screen
        {
            std::lock_guard<std::mutex> lock(rendererMutex);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);
        }
        
        // If we're rendering, show progress
        if (isRendering) {
            // Show progress bar
            renderProgressBar(renderer, renderProgress, SCREEN_HEIGHT,
                             SCREEN_WIDTH / 2, 30,
                             SCREEN_WIDTH / 4, SCREEN_HEIGHT / 2 - 15);
            
            // Display status message
            int percent = renderProgress * 100 / SCREEN_HEIGHT;
            std::string statusMsg = "Rendering Mandelbrot Set... " + std::to_string(percent) + "% complete";
            renderPopup(renderer, font, statusMsg);
            
            // Display additional info
            std::string zoomInfo = "Zoom Level: " + std::to_string(viewport.getZoomLevel()) +
                                 ", Using " + std::to_string(NUM_THREADS) + " threads" +
                                 ", Iterations: " + std::to_string(maxIterations);
            
            SDL_Color textColor = {200, 200, 200, 255};
            SDL_Surface* textSurface = TTF_RenderText_Blended(font, zoomInfo.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = nullptr;
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                }
                int textW = 0, textH = 0;
                SDL_QueryTexture(textTexture, NULL, NULL, &textW, &textH);
                
                SDL_Rect textRect = {(SCREEN_WIDTH - textW) / 2, SCREEN_HEIGHT / 2 + 20, textW, textH};
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
                }
                
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    SDL_DestroyTexture(textTexture);
                }
                SDL_FreeSurface(textSurface);
            }
        } else if (mandelbrotTexture) {
            // Render the current fractal if available
            {
                std::lock_guard<std::mutex> lock(rendererMutex);
                SDL_RenderCopy(renderer, mandelbrotTexture, NULL, NULL);
            }
            
            // Draw selection rectangle if active
            if (isSelecting) {
                renderSelectionRect(renderer, selectionStartX, selectionStartY, selectionEndX, selectionEndY);
            }
            
            // Show popup if necessary
            Uint32 currentTime = SDL_GetTicks();
            if (showPopup && currentTime - popupStartTime < POPUP_DURATION) {
                renderPopup(renderer, font, popupMessage);
            } else {
                showPopup = false;
            }
            
            // Draw help information in the corner
            std::string helpText = "F: Fullscreen | ESC: Exit | R: Reset | G: AI | +/-: Zoom | T: Threading | A: Auto-Explore";
            SDL_Color textColor = {200, 200, 200, 255};
            SDL_Surface* textSurface = TTF_RenderText_Blended(font, helpText.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = nullptr;
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                }
                int textW = 0, textH = 0;
                SDL_QueryTexture(textTexture, NULL, NULL, &textW, &textH);
                
                // Dark background for text legibility
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
                    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 180);
                    SDL_Rect bgRect = {10, SCREEN_HEIGHT - textH - 15, textW + 10, textH + 10};
                    SDL_RenderFillRect(renderer, &bgRect);
                    
                    // Render text
                    SDL_Rect textRect = {15, SCREEN_HEIGHT - textH - 10, textW, textH};
                    SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
                    
                    SDL_DestroyTexture(textTexture);
                }
                SDL_FreeSurface(textSurface);
            }
            
            // Draw threading and performance info
            std::string threadInfo = "Using " + std::to_string(NUM_THREADS) + " threads | " + 
                                    "Zoom Level: " + std::to_string(viewport.getZoomLevel()) + " | " +
                                    "Iterations: " + std::to_string(maxIterations);
            
            // Add GenAI mode indicator if enabled
            if (genAIExplorationMode) {
                threadInfo += " | GenAI: ENABLED";
            }
            
            textSurface = TTF_RenderText_Blended(font, threadInfo.c_str(), textColor);
            if (textSurface) {
                SDL_Texture* textTexture = nullptr;
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                }
                int textW = 0, textH = 0;
                SDL_QueryTexture(textTexture, NULL, NULL, &textW, &textH);
                
                // Dark background for text legibility
                {
                    std::lock_guard<std::mutex> lock(rendererMutex);
                    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
                    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 180);
                    SDL_Rect bgRect = {SCREEN_WIDTH - textW - 20, 10, textW + 10, textH + 10};
                    SDL_RenderFillRect(renderer, &bgRect);
                    
                    // Render text
                    SDL_Rect textRect = {SCREEN_WIDTH - textW - 15, 15, textW, textH};
                    SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
                    
                    SDL_DestroyTexture(textTexture);
                }
                SDL_FreeSurface(textSurface);
            }
        }
        
        {
            std::lock_guard<std::mutex> lock(rendererMutex);
            SDL_RenderPresent(renderer);
        }
        SDL_Delay(10);
    }
    
    // Clean up
    if (mandelbrotTexture) {
        std::lock_guard<std::mutex> lock(rendererMutex);
        SDL_DestroyTexture(mandelbrotTexture);
    }
    
    TTF_CloseFont(font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}