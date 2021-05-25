#include <assert.h>
#include <dbg.h>
#include <inttypes.h>
#include <stdio.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>


#include "instrument.h"
#include "timing.h"
#include "util.h"

#include "impl_ref/impl.hpp"
#include "impl_ref/scene.hpp"

#include "impl_opt0/impl.hpp"
#include "impl_opt1/impl.hpp"
#include "impl_opt3/impl.hpp"
#include "impl_opt4/impl.hpp"
#include "impl_opt5/impl.hpp"

flops_t flops_counter;

typedef void (*fp_render_init)(std::string);
typedef void (*fp_render)(int, int, float*);

fp_render_init fun_render_init;
fp_render fun_render;

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

void run(int width, int height, std::string input, std::string output) {
    auto pixels = std::make_unique<float[]>(height * width * 3);

    ins_rst();
    fun_render_init(input);
    ins_dump("Setup");

    ins_rst();

#ifndef DO_INSTRUMENT
    timing_start();
#endif
    fun_render(width, height, pixels.get());
#ifndef DO_INSTRUMENT
    timing_t timing = timing_stop();
#endif

    ins_dump(NULL);

    fprintf(stderr, "Width: %d\n", width);
    fprintf(stderr, "Height: %d\n", height);
#ifdef DO_INSTRUMENT
    fprintf(stderr, "Flops: %" PRIu64 "\n", ins_total());
#else
    fprintf(stderr, "Cycles: %" PRIu64 "\n", timing.cycles);
    fprintf(stderr, "Microseconds: %" PRIu64 "\n", timing.usec);
    fprintf(stderr, "Seconds: %.2f\n", timing.usec * 1.f / 1e6);
#endif

    if (!output.empty()) {
        std::ofstream o;
        o.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        bool hdr;
        if (output.find(".") != std::string::npos) {
            std::string fileformat = output.substr(output.find_last_of("."));
            if (fileformat.compare(".pfm") == 0) {
                hdr = true;
            } else if (fileformat.compare(".ppm") == 0) {
                hdr = false;
            } else {
                std::cerr << "Unsupported file format '" << fileformat << "', defaulting to '.pfm'" << std::endl;
                hdr = true;
                output.append(".pfm");
            }
        } else {
            std::cerr << "File format not specified, defaulting to '.pfm'" << std::endl;
            hdr = true;
            output.append(".pfm");
        }

        try {
            o.open(output);
        } catch (std::system_error& e) {
            std::cerr << "Failed to open output file '" << output << "': " << strerror(errno) << std::endl;
            throw;
        }

        dump_image(o, width, height, pixels.get(), hdr);
        o.close();
    }
}

/**
 * Sets the two function pointers fun_render_init and fun_render to the one for the requested implementation.
 *
 * When adding a new implementation add the corresponding case here.
 */
void set_render_fp(const std::string& impl) {
    if (impl == "ref") {
        fun_render_init = &impl::ref::render_init;
        fun_render = &impl::ref::render;
    } else if (impl == "opt0") {
        fun_render_init = &impl::opt0::render_init;
        fun_render = &impl::opt0::render;
    } else if (impl == "opt1") {
        fun_render_init = &impl::opt1::render_init;
        fun_render = &impl::opt1::render;
    } else if (impl == "opt3") {
        fun_render_init = &impl::opt3::render_init;
        fun_render = &impl::opt3::render;
    } else if (impl == "opt4") {
        fun_render_init = &impl::opt4::render_init;
        fun_render = &impl::opt4::render;
    } else if (impl == "opt5") {
        fun_render_init = &impl::opt5::render_init;
        fun_render = &impl::opt5::render;
    } else {
        throw std::runtime_error("Unknown implementation '" + impl + "'");
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <impl> <input> <width> <height> [<output>] [<reference>]\n", argv[0]);
        return EXIT_FAILURE;
    }

    std::string impl = argv[1];
    std::string input = argv[2];
    std::string width_str = argv[3];
    std::string height_str = argv[4];
    std::string output = argc > 5 ? argv[5] : "";
    std::string reference = argc > 6 ? argv[6] : "";

    set_render_fp(impl);

    /*
     * Load the reference scene. All other implementations will use this to derive their own scenes.
     */
    impl::ref::load_scene(input);

    int width = std::stoi(width_str);
    int height = std::stoi(height_str);

    run(width, height, input, output);

    // compare against reference image
    if (!output.empty() && !reference.empty()) {
        validate_output(5, 0.01, reference, output, height, width);
    }

    return 0;
}
