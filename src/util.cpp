#include "util.h"

#include <algorithm>

void dump_image_ldr(std::ostream& out, int width, int height, const float* pixels) {
    out << "P6\n" << width << " " << height << "\n255\n";

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                unsigned char channel = std::clamp(pixels[3 * (width * j + i) + k], 0.f, 1.f) * 255.f;
                out << channel;
            }
        }
    }
}

void dump_image_hdr(std::ostream& out, int width, int height, const float* pixels) {
    out << "PF\n" << width << " " << height << "\n-1.0\n";

    for (int j = height - 1; j > -1; j--) {
        for (int i = 0; i < width; i++) {
            out.write((char*)&pixels[3 * (width * j + i)], 3 * sizeof(float));
        }
    }
}

void dump_image(std::ostream& out, int width, int height, const float* pixels, bool hdr) {
    if (hdr) {
        // store in binary .pfm format
        dump_image_hdr(out, width, height, pixels);
    } else {
        // store in binary .ppm format
        dump_image_ldr(out, width, height, pixels);
    }
}
