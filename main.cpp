#include <iostream>

#include "instrument.h"

flops_t flops_counter;

int main(void) {
    ins_rst();
    float x = FMA(2, 9, 11);

    std::cout << x << std::endl;
}
