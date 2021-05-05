#include "util.h"

#include <algorithm>
#include <cmath>

/**
 * Writes a PPM file.
 */
void dump_image_ldr(std::ostream& out, int width, int height, const float* pixels, float exposure) {
    out << "P6\n" << width << " " << height << "\n255\n";

    // gamma correction according to wikipedia
    float invgamma = 0.45f; // inverse of gamma=2.2f

    // turn off gamma correction by default
    bool gammacorrection = false;
    invgamma = (gammacorrection) ? invgamma : 1.0f;

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                float channel = exposure * pixels[3 * (width * j + i) + k];
                unsigned char encoding = std::clamp(std::pow(channel, invgamma), 0.f, 1.f) * 255.f;
                out << encoding;
            }
        }
    }
}

/**
 * Writes a PFM file.
 */
void dump_image_hdr(std::ostream& out, int width, int height, const float* pixels, float exposure) {
    out << "PF\n" << width << " " << height << "\n-1.0\n";

    for (int j = height - 1; j > -1; j--) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                float channel = exposure * pixels[3 * (width * j + i) + k];
                out.write((char*)&channel, sizeof(float));
            }
        }
    }
}

/**
 * Writes either a PPM (default) or PFM file.
 */
void dump_image(std::ostream& out, int width, int height, const float* pixels, bool hdr) {
    /* TODO: Possible Strategies to determine exposure:
        - user defined by exposure parameter in camera
        - choose exposure empirically
        - map max_val to 1
        - map 95th percentile to 1
    */

    float exposure = 1250.f; //

    /*
     float max_val = 0.f;
     for (int i = 0; i < width * height * 3; i++){
         max_val = std::max(max_val, pixels[i]);
     }
     float exposure = 1.f / max_val * 1.76;
     */

    if (hdr) {
        // store in binary .pfm format
        dump_image_hdr(out, width, height, pixels, exposure);
    } else {
        // store in binary .ppm format
        dump_image_ldr(out, width, height, pixels, exposure);
    }
}
