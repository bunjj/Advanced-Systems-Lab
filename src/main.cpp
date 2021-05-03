#include <assert.h>
#include <dbg.h>
#include <inttypes.h>
#include <stdio.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>

#include "impl_ref/geometry.h"
#include "impl_ref/impl.hpp"
#include "impl_ref/scene.hpp"
#include "instrument.h"
#include "timing.h"
#include "util.h"

flops_t flops_counter;

// Image Files and Output Validation {{{
/**
 * Reads a PPM file into the given buffer.
 * Adapted from here: http://www.cplusplus.com/forum/general/208835/
 */
static void read_ppm(std::string filename, unsigned char* pixels_in) {
    FILE* fp = fopen(filename.c_str(), "rb");

    assert(fp);

    // read header
    char pSix[10];
    fscanf(fp, "%s", pSix);

    // check if it is a PPM file
    if (strncmp(pSix, "P6", 10) != 0) {
        std::cerr << "Input file is not PPM!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // read the rest of header
    int width, height;
    int maximum;
    fscanf(fp, "%d\n %d\n", &width, &height);
    fscanf(fp, "%d\n", &maximum);

    // unformatted read of binary pixel data
    fread(pixels_in, sizeof(unsigned char), width * height * 3, fp);

    if (ferror(fp)) {
        std::cerr << "error while reading file" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (feof(fp)) {
        std::cerr << "EOF reached earlier than expected" << std::endl;
        exit(EXIT_FAILURE);
    }

    // close file
    fclose(fp);
}

/**
 * Returns the fraction of pixels that are significantly different from the reference image.
 * Significantly different means the difference in RGB values (summed) is greater than rgb_tol.
 */
static float compare_pixels(int rgb_tol, unsigned char* pixels_reference, unsigned char* pixels_out, size_t size) {
    int num_different = 0;

    for (size_t i = 0; i < size; i++) {
        unsigned char ro = pixels_out[3 * i];
        unsigned char go = pixels_out[3 * i + 1];
        unsigned char bo = pixels_out[3 * i + 2];

        unsigned char rr = pixels_reference[3 * i];
        unsigned char gr = pixels_reference[3 * i + 1];
        unsigned char br = pixels_reference[3 * i + 2];

        int difference = abs(ro - rr) + abs(go - gr) + abs(bo - br);
        if (difference > rgb_tol) {
            num_different++;
        }
    }

    // return fraction of pixels that are different
    return (float)num_different / size;
}

/**
 * Asserts that no more than overall_tol (= fraction) pixels are significantly different from the reference image.
 * Significantly different means that the difference in RGB values (range 0..255) (summed) is greater than rgb_tol.
 */
static void validate_output(
    int rgb_tol, float overall_tol, std::string ref_filename, std::string out_filename, int height, int width) {
    int size = height * width;

    // read the two images
    auto pixels_ref = std::make_unique<unsigned char[]>(size * 3);
    read_ppm(ref_filename, pixels_ref.get());
    auto pixels_out = std::make_unique<unsigned char[]>(size * 3);
    read_ppm(out_filename, pixels_out.get());

    float fraction_different = compare_pixels(rgb_tol, pixels_ref.get(), pixels_out.get(), height * width);
    if (fraction_different > overall_tol) {
        std::cerr << "OUTPUT VALIDATION FAILED: " << std::setprecision(4) << fraction_different * 100 << "% different"
                  << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "output validation OK: " << std::setprecision(4) << fraction_different * 100 << "% different"
                  << std::endl;
    }
}

// }}}

void run(int width, int height, std::string output) {
    auto pixels = std::make_unique<float[]>(height * width * 3);

    ins_rst();
    impl::ref::render_init();
    ins_dump("Setup");

    ins_rst();

    timing_start();
    impl::ref::render(width, height, pixels.get());
    timing_t timing = timing_stop();

    ins_dump(NULL);

    fprintf(stderr, "Width: %d\n", width);
    fprintf(stderr, "Height: %d\n", height);
    fprintf(stderr, "Flops: %" PRIu64 "\n", ins_total());
    fprintf(stderr, "Cycles: %" PRIu64 "\n", timing.cycles);
    fprintf(stderr, "Microseconds: %" PRIu64 "\n", timing.usec);
    fprintf(stderr, "Seconds: %.2f\n", timing.usec * 1.f / 1e6);

    if (!output.empty()) {
        std::ofstream o;
        o.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try {
            o.open(output);
        } catch (std::system_error& e) {
            std::cerr << "Failed to open output file '" << output << "': " << strerror(errno) << std::endl;
            throw;
        }

        dump_image(o, width, height, pixels.get());
        o.close();
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input> <width> <height> [<output>] [<reference>]\n", argv[0]);
        return EXIT_FAILURE;
    }

    std::string input = argv[1];
    std::string width_str = argv[2];
    std::string height_str = argv[3];
    std::string output = argc > 4 ? argv[4] : "";
    std::string reference = argc > 5 ? argv[5] : "";

    impl::ref::load_scene(input);

    int width = std::stoi(width_str);
    int height = std::stoi(height_str);

    run(width, height, output);

    // compare against reference image
    if (!output.empty() && !reference.empty()) {
        validate_output(5, 0.01, reference, output, height, width);
    }

    return 0;
}
