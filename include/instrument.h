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
    uint64_t max; //< For max and min
    uint64_t pow;
    uint64_t tan;
    uint64_t sphere;
    uint64_t plane;
    uint64_t box;
    uint64_t torus;
    uint64_t cone;
    uint64_t octa;
    uint64_t sphere_r;
    uint64_t box_r;
    uint64_t torus_r;
    uint64_t cone_r;
    uint64_t octa_r;
    uint64_t sphere_n;
    uint64_t plane_n;
    uint64_t box_n;
    uint64_t torus_n;
    uint64_t cone_n;
    uint64_t octa_n;
} flops_t;

extern flops_t flops_counter;

#ifndef DO_INSTRUMENT

#define ins_dump(title) NOP
#define ins_total() UINT64_C(0)
#define ins_rst() NOP
#define INS_INC1(name, offset) NOP

#else

static inline uint64_t ins_total() {
    return flops_counter.add + flops_counter.mul + 2 * flops_counter.fma + flops_counter.div + flops_counter.sqrt +
           flops_counter.abs + flops_counter.cmp + flops_counter.max + flops_counter.pow + flops_counter.tan;
}

static inline void ins_dump(const char* title) {
    fprintf(stderr, "====================\n");
    if (title) {
        fprintf(stderr, "%s\n", title);
    }
    fprintf(stderr, "=== FLOPS COUNTER ==\n");
    fprintf(stderr, "ADD      : %12" PRIu64 "\n", flops_counter.add);
    fprintf(stderr, "MUL      : %12" PRIu64 "\n", flops_counter.mul);
    fprintf(stderr, "FMA      : %12" PRIu64 "\n", flops_counter.fma);
    fprintf(stderr, "DIV      : %12" PRIu64 "\n", flops_counter.div);
    fprintf(stderr, "SQRT     : %12" PRIu64 "\n", flops_counter.sqrt);
    fprintf(stderr, "ABS      : %12" PRIu64 "\n", flops_counter.abs);
    fprintf(stderr, "CMP      : %12" PRIu64 "\n", flops_counter.cmp);
    fprintf(stderr, "MAX      : %12" PRIu64 "\n", flops_counter.max);
    fprintf(stderr, "POW      : %12" PRIu64 "\n", flops_counter.pow);
    fprintf(stderr, "TAN      : %12" PRIu64 "\n", flops_counter.tan);
    fprintf(stderr, "SPHERE   : %12" PRIu64 "\n", flops_counter.sphere);
    fprintf(stderr, "SPHERE_R : %12" PRIu64 "\n", flops_counter.sphere_r);
    fprintf(stderr, "PLANE    : %12" PRIu64 "\n", flops_counter.plane);
    fprintf(stderr, "BOX      : %12" PRIu64 "\n", flops_counter.box);
    fprintf(stderr, "BOX_R    : %12" PRIu64 "\n", flops_counter.box_r);
    fprintf(stderr, "TORUS    : %12" PRIu64 "\n", flops_counter.torus);
    fprintf(stderr, "TORUS_R  : %12" PRIu64 "\n", flops_counter.torus_r);
    fprintf(stderr, "CONE     : %12" PRIu64 "\n", flops_counter.cone);
    fprintf(stderr, "CONE_R   : %12" PRIu64 "\n", flops_counter.cone_r);
    fprintf(stderr, "OCTA     : %12" PRIu64 "\n", flops_counter.octa);
    fprintf(stderr, "OCTA_R   : %12" PRIu64 "\n", flops_counter.octa_r);
    fprintf(stderr, "SPHERE_N : %12" PRIu64 "\n", flops_counter.sphere_n);
    fprintf(stderr, "PLANE_N  : %12" PRIu64 "\n", flops_counter.plane_n);
    fprintf(stderr, "BOX_N    : %12" PRIu64 "\n", flops_counter.box_n);
    fprintf(stderr, "TORUS_N  : %12" PRIu64 "\n", flops_counter.torus_n);
    fprintf(stderr, "CONE_N   : %12" PRIu64 "\n", flops_counter.cone_n);
    fprintf(stderr, "OCTA_N   : %12" PRIu64 "\n", flops_counter.octa_n);
    fprintf(stderr, "=======================\n");
    fprintf(stderr, "TOTAL    : %12" PRIu64 "\n", ins_total());
    fprintf(stderr, "=======================\n");
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
#define INS_TAN INS_INC(tan)

#define FADD(x, y) (INS_ADD, ((x) + (y)))
#define FMUL(x, y) (INS_MUL, ((x) * (y)))
#define FMA(x, y, z) (INS_FMA, (((x) * (y)) + (z)))
#define FDIV(x, y) (INS_DIV, ((x) / (y)))
#define FSQRT(x) (INS_SQRT, sqrtf((x)))
#define FABS(x) (INS_ABS, fabsf((x)))

// TODO add macros for vector operations
