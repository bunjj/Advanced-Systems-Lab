#pragma once

#include <stdint.h>

typedef struct {
    /**
     * Runtime in cycles.
     */
    uint64_t cycles;
    /**
     * Runtime in microseconds.
     */
    uint64_t usec;
} timing_t;

void timing_start();
timing_t timing_stop();
