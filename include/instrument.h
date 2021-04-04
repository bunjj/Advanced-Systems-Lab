#pragma once

#include <stdint.h>
#include <string.h>

#define NOP do {} while(0)

typedef struct {
    uint64_t add;
    uint64_t mul;
    uint64_t fma;
    uint64_t div;
} flops_t;

extern flops_t flops_counter;

#define DO_INSTRUMENT

#ifndef DO_INSTRUMENT

#define ins_rst() NOP
#define INS_INC1(name, offset) NOP

#else
static inline void ins_rst(void) {
    memset(&flops_counter, 0, sizeof(flops_counter));
}

#define INS_INC1(name, offset) (flops_counter.name += (offset))
#endif

#define INS_INC(name) INS_INC1(name, 1)
#define INS_ADD INS_INC(add)
#define INS_MUL INS_INC(mul)
#define INS_FMA INS_INC(fma)
#define INS_DIV INS_INC(div)

#define FADD(x, y) (INS_ADD, ((x) + (y)))
#define FMUL(x, y) (INS_MUL, ((x) * (y)))
#define FMA(x, y, z) (INS_FMA, (((x) * (y)) + (z)))
#define FDIV(x, y) (INS_DIV, ((x) / (y)))
// TODO add macros for vector operations
