#include <inttypes.h>
#include <stdio.h>

#include <iostream>
#include <memory>

static void dump_image(int width, int height, const float* pixels) {
    printf("P3\n%d %d\n255\n", width, height);

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                uint8_t channel = pixels[3 * (width * j + i) + k] * 255;
                printf("%03" PRIu8 " ", channel);
            }
        }

        printf("\n");
    }
}

int main(void) {
    // Height of the resulting image in pixels
    int height = 720;
    int width = 1280;

    auto pixels = std::make_unique<float[]>(height * width * 3);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            pixels[3 * (width * i + j)] = static_cast<float>(j) / width;
            pixels[3 * (width * i + j) + 1] = static_cast<float>(i) / height;
            pixels[3 * (width * i + j) + 2] = 0;
        }
    }

    dump_image(width, height, pixels.get());

    return 0;
}
