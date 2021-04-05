#pragma once

#include <math.h>
#include <iostream>
#include <iomanip>

#include "instrument.h"

extern float M_PI_F;
#define TO_RAD(angle) ((angle) / 180.0f * M_PI_F)

// Vec {{{

struct vec {
    float x;
    float y;
    float z;
};

std::ostream& operator<<(std::ostream& out, const vec& v);

#define VEC_OP(v1, v2, OP) \
    vec { (v1).x OP(v2).x, (v1).y OP(v2).y, (v1).z OP(v2).z }

static inline vec vec_add(vec v1, vec v2) { 
    INS_INC1(add, 3);
    return VEC_OP(v1, v2, +);
}

static inline vec vec_sub(vec v1, vec v2) { 
    INS_INC1(add, 3);
    return VEC_OP(v1, v2, -);
}

static inline vec vec_scale(vec v1, float factor) {
    INS_INC1(mul, 3);
    return VEC_OP(v1, (vec{factor, factor, factor}), *);
}

static inline float vec_dot(vec v1, vec v2) {
    INS_INC1(add, 2);
    INS_INC1(mul, 3);
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

static inline float vec_dot2(vec v) {
    return vec_dot(v, v);
}

static inline float vec_length(vec v) { 
    return FSQRT(vec_dot2(v));
}

static inline vec vec_normalize(vec v) {
    float len = vec_length(v);
    INS_INC1(div, 3);
    return vec{v.x / len, v.y / len, v.z / len};
}

static inline vec vec_abs(vec v) { 
    INS_INC1(abs, 3);
    return {fabsf(v.x), fabsf(v.y), fabsf(v.z)}; 
}

static inline vec vec_max(vec v, float val) {
    INS_INC1(max, 3);
    return {std::max(v.x, val), std::max(v.y, val), std::max(v.z, val)};
}
// }}}

// Vec2 {{{
struct vec2 {
    float x;
    float y;
};

#define VEC2_OP(v1, v2, OP) \
    vec2{ (v1).x OP(v2).x, (v1).y OP(v2).y }

static inline vec2 vec2_sub(vec2 v1, vec2 v2) {
    INS_INC1(add, 2);
    return VEC2_OP(v1, v2, -);
}

static inline vec2 vec2_add(vec2 v1, vec2 v2) {
    INS_INC1(add, 2);
    return VEC2_OP(v1, v2, +);
}

static inline vec2 vec2_scale(vec2 v, float f) {
    INS_INC1(mul, 2);
    return VEC2_OP(v, (vec2{f, f}), *);
}

static inline float vec2_dot(vec2 v1, vec2 v2) {
    INS_INC(add);
    INS_INC1(mul, 2);
    return v1.x * v2.x + v1.y * v2.y;
}

static inline float vec2_dot2(vec2 v) {
    return vec2_dot(v, v);
}

static inline float vec2_length(vec2 v) {
    return FSQRT(vec2_dot2(v));
}
// }}}

// Vec4 {{{

struct vec4 {
    float x;
    float y;
    float z;
    float t;
};

static inline vec4 vec4_init(float* values) {
    return {values[0], values[1], values[2], values[3]};
}

static inline float vec4_dot(vec4 v1, vec4 v2) {
    INS_INC1(add, 3);
    INS_INC1(mul, 4);
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.t * v2.t;
}

static inline vec4 vec4_from_point(vec p) { return {p.x, p.y, p.z, 1}; }

static inline vec4 vec4_from_dir(vec p) { return {p.x, p.y, p.z, 0}; }

static inline vec vec4_to_vec(vec4 p) { return {p.x, p.y, p.z}; }

static inline vec4 get_base_vec4(int idx) {
    switch (idx) {
        case 0:
            return {1, 0, 0, 0};
        case 1:
            return {0, 1, 0, 0};
        case 2:
            return {0, 0, 1, 0};
        case 3:
            return {0, 0, 0, 1};
        default:
            return {0, 0, 0, 0};
    }
}

// }}}

// M44 {{{

/**
 * Defines a 4x4 matrix.
 *
 * Used for transformations
 */
struct m44 {
    /**
     * val[i][j] is the i-th row, j-th column
     *
     * This means it is stored in row-major format.
     */
    float val[4][4];

    m44(){};

    m44(vec4 e1, vec4 e2, vec4 e3, vec4 e4) {
        val[0][0] = e1.x;
        val[1][0] = e1.y;
        val[2][0] = e1.z;
        val[3][0] = e1.t;
        val[0][1] = e2.x;
        val[1][1] = e2.y;
        val[2][1] = e2.z;
        val[3][1] = e2.t;
        val[0][2] = e3.x;
        val[1][2] = e3.y;
        val[2][2] = e3.z;
        val[3][2] = e3.t;
        val[0][3] = e4.x;
        val[1][3] = e4.y;
        val[2][3] = e4.z;
        val[3][3] = e4.t;
    };
};

static const m44 identity =
    m44(get_base_vec4(0), get_base_vec4(1), get_base_vec4(2), get_base_vec4(3));

m44 get_transf_matrix(vec pos, vec rot);
m44 get_rot_matrix(vec rot);
m44 m44_inv(m44 m);
vec4 m44_mul_vec(m44 m, vec4 v);
std::ostream& operator<<(std::ostream& out, const m44& m);
// }}}
