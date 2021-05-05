#include "timing.h"

#include <sys/time.h>
#include <x86intrin.h>

static timing_t start;

static timing_t get_timing_data() {
    struct timeval tval;
    gettimeofday(&tval, NULL);

    timing_t data;
    data.cycles = _rdtsc();
    data.usec = tval.tv_sec * 1e6 + tval.tv_usec;

    return data;
}

void timing_start() {
    start = get_timing_data();
}

timing_t timing_stop() {
    timing_t end = get_timing_data();

    timing_t res;
    res.cycles = end.cycles - start.cycles;
    res.usec = end.usec - start.usec;

    return res;
}
