#include <iostream>

#include <nlohmann/json.hpp>

#include "instrument.h"
#include "geometry.h"

flops_t flops_counter;

using json = nlohmann::json;

int main(void) {
    ins_rst();
    float x = FMA(2, 9, 11);

    std::cout << x << std::endl;

    ins_dump();
}
