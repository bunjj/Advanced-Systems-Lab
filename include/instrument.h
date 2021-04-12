#pragma once

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define NOP ((void)0)

typedef struct {
    uint64_t add;
    uint64_t mul;
    uint64_t fma;
    uint64_t div;
    uint64_t sqrt;
    uint64_t abs;
    uint64_t cmp;
    uint64_t pow;
    uint64_t sphere;
    uint64_t plane;
    uint64_t box;
    uint64_t torus;
    uint64_t cone;
    uint64_t octa;
} flops_t;

extern flops_t flops_counter;

#ifndef DO_INSTRUMENT

#define ins_dump(title) NOP
#define ins_rst() NOP
#define INS_INC1(name, offset) NOP

#else

static inline void ins_dump(const char* title) {
    if (title) {
        fprintf(stderr, "%s\n", title);
    }
    fprintf(stderr, "====================\n");
    fprintf(stderr, "=== FLOPS COUNTER ==\n");
    fprintf(stderr, "ADD   : %12" PRIu64 "\n", flops_counter.add);
    fprintf(stderr, "MUL   : %12" PRIu64 "\n", flops_counter.mul);
    fprintf(stderr, "FMA   : %12" PRIu64 "\n", flops_counter.fma);
    fprintf(stderr, "DIV   : %12" PRIu64 "\n", flops_counter.div);
    fprintf(stderr, "SQRT  : %12" PRIu64 "\n", flops_counter.sqrt);
    fprintf(stderr, "ABS   : %12" PRIu64 "\n", flops_counter.abs);
    fprintf(stderr, "CMP   : %12" PRIu64 "\n", flops_counter.cmp);
    fprintf(stderr, "POW   : %12" PRIu64 "\n", flops_counter.pow);
    fprintf(stderr, "SPHERE: %12" PRIu64 "\n", flops_counter.sphere);
    fprintf(stderr, "PLANE : %12" PRIu64 "\n", flops_counter.plane);
    fprintf(stderr, "BOX   : %12" PRIu64 "\n", flops_counter.box);
    fprintf(stderr, "TORUS : %12" PRIu64 "\n", flops_counter.torus);
    fprintf(stderr, "CONE  : %12" PRIu64 "\n", flops_counter.cone);
    fprintf(stderr, "OCTA  : %12" PRIu64 "\n", flops_counter.octa);
    fprintf(stderr, "====================\n");
}

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
#define INS_SQRT INS_INC(sqrt)
#define INS_ABS INS_INC(abs)
#define INS_CMP INS_INC(cmp)
#define INS_POW INS_INC(pow)

#define FADD(x, y) (INS_ADD, ((x) + (y)))
#define FMUL(x, y) (INS_MUL, ((x) * (y)))
#define FMA(x, y, z) (INS_FMA, (((x) * (y)) + (z)))
#define FDIV(x, y) (INS_DIV, ((x) / (y)))
#define FSQRT(x) (INS_SQRT, sqrtf((x)))
#define FABS(x) (INS_ABS, fabsf((x)))

// TODO add macros for vector operations
