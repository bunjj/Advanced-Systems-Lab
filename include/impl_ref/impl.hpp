#pragma once

#include <iostream>

namespace impl::ref {
    void render_init(std::string input);
    void render(int width, int height, float* pixels);
} // namespace impl::ref
