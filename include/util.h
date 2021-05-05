#pragma once

#include <iostream>

void dump_image_ldr(std::ostream& out, int width, int height, const float* pixels, float exposure);
void dump_image_hdr(std::ostream& out, int width, int height, const float* pixels, float exposure);
void dump_image(std::ostream& out, int width, int height, const float* pixels, bool hdr = false);
