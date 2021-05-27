#pragma once

#include <iostream>

namespace impl::vec2 {
    void render_init(std::string input);
    void render(int width, int height, float* pixels);
} // namespace impl::vec2
