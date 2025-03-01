#include "../Properties/Complex Number Handling (z, c).h"
#include "../Properties/Iteration Parameters.h"

bool checkDivergence(const Complex& z, double threshold) {
    // A point is considered to diverge if |z| > threshold
    // For optimization, we compare |z|² > threshold² to avoid square root calculation
    double thresholdSquared = threshold * threshold;
    return z.magnitudeSquared() > thresholdSquared;
}