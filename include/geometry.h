#pragma once

#include <math.h>
#include <iostream>
#include <iomanip>

#include "instrument.h"

static float M_PI_F = M_PI;
#define TO_RAD(angle) ((angle) / 180.0f * M_PI_F)

// Vec {{{

struct vec {
    float x;
    float y;
    float z;
};

std::ostream& operator<<(std::ostream& out, const vec& v) {
    out << "{" << v.x << ", " << v.y << ", " << v.z << "}";
    return out;
}

#define VEC_OP(v1, v2, OP) \
    vec { (v1).x OP(v2).x, (v1).y OP(v2).y, (v1).z OP(v2).z }

inline vec vec_add(vec v1, vec v2) { 
    INS_INC1(add, 3);
    return VEC_OP(v1, v2, +);
}

inline vec vec_sub(vec v1, vec v2) { 
    INS_INC1(add, 3);
    return VEC_OP(v1, v2, -);
}

inline vec vec_scale(vec v1, float factor) {
    INS_INC1(mul, 3);
    return VEC_OP(v1, (vec{factor, factor, factor}), *);
}

inline float vec_dot(vec v1, vec v2) {
    INS_INC1(add, 2);
    INS_INC1(mul, 3);
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline float vec_dot2(vec v) {
    return vec_dot(v, v);
}

inline float vec_length(vec v) { 
    return FSQRT(vec_dot2(v));
}

inline vec vec_normalize(vec v) {
    float len = vec_length(v);
    INS_INC1(div, 3);
    return vec{v.x / len, v.y / len, v.z / len};
}

inline vec vec_abs(vec v) { 
    INS_INC1(abs, 3);
    return {fabsf(v.x), fabsf(v.y), fabsf(v.z)}; 
}

inline vec vec_max(vec v, float val) {
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

inline vec2 vec2_sub(vec2 v1, vec2 v2) {
    INS_INC1(add, 2);
    return VEC2_OP(v1, v2, -);
}

inline vec2 vec2_add(vec2 v1, vec2 v2) {
    INS_INC1(add, 2);
    return VEC2_OP(v1, v2, +);
}

inline vec2 vec2_scale(vec2 v, float f) {
    INS_INC1(mul, 2);
    return VEC2_OP(v, (vec2{f, f}), *);
}

inline float vec2_dot(vec2 v1, vec2 v2) {
    INS_INC(add);
    INS_INC1(mul, 2);
    return v1.x * v2.x + v1.y * v2.y;
}

inline float vec2_dot2(vec2 v) {
    return vec2_dot(v, v);
}

inline float vec2_length(vec2 v) {
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

inline vec4 vec4_init(float* values) {
    return {values[0], values[1], values[2], values[3]};
}

inline float vec4_dot(vec4 v1, vec4 v2) {
    INS_INC1(add, 3);
    INS_INC1(mul, 4);
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.t * v2.t;
}

inline vec4 vec4_from_point(vec p) { return {p.x, p.y, p.z, 1}; }

inline vec4 vec4_from_dir(vec p) { return {p.x, p.y, p.z, 0}; }

inline vec vec4_to_vec(vec4 p) { return {p.x, p.y, p.z}; }

inline vec4 get_base_vec4(int idx) {
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

std::ostream& operator<<(std::ostream& out, const m44& m) {
    int width = 8;
    out << std::fixed << std::setprecision(2) << std::right << std::endl;
    // clang-format off
    out << "┌ " << std::setw(width) << m.val[0][0] << " " << std::setw(width) << m.val[0][1] << " " << std::setw(width) << m.val[0][2] << " " << std::setw(width) << m.val[0][3] << " ┐" << std::endl;
    out << "│ " << std::setw(width) << m.val[1][0] << " " << std::setw(width) << m.val[1][1] << " " << std::setw(width) << m.val[1][2] << " " << std::setw(width) << m.val[1][3] << " │" << std::endl;
    out << "│ " << std::setw(width) << m.val[2][0] << " " << std::setw(width) << m.val[2][1] << " " << std::setw(width) << m.val[2][2] << " " << std::setw(width) << m.val[2][3] << " │" << std::endl;
    out << "└ " << std::setw(width) << m.val[3][0] << " " << std::setw(width) << m.val[3][1] << " " << std::setw(width) << m.val[3][2] << " " << std::setw(width) << m.val[3][3] << " ┘" << std::endl;
    // clang-format on
    return out;
}

/**
 * Calculates m * v
 */
vec4 m44_mul_vec(m44 m, vec4 v) {
    float x = vec4_dot(vec4_init(m.val[0]), v);
    float y = vec4_dot(vec4_init(m.val[1]), v);
    float z = vec4_dot(vec4_init(m.val[2]), v);
    float t = vec4_dot(vec4_init(m.val[3]), v);

    return {x, y, z, t};
}

m44 m44_inv(m44 m) {
    float a11 = m.val[0][0];
    float a12 = m.val[0][1];
    float a13 = m.val[0][2];
    float a14 = m.val[0][3];
    float a21 = m.val[1][0];
    float a22 = m.val[1][1];
    float a23 = m.val[1][2];
    float a24 = m.val[1][3];
    float a31 = m.val[2][0];
    float a32 = m.val[2][1];
    float a33 = m.val[2][2];
    float a34 = m.val[2][3];
    float a41 = m.val[3][0];
    float a42 = m.val[3][1];
    float a43 = m.val[3][2];
    float a44 = m.val[3][3];

    INS_INC1(add, 18);
    INS_INC1(mul, 36);

    float A2323 = a33 * a44 - a34 * a43;
    float A1323 = a32 * a44 - a34 * a42;
    float A1223 = a32 * a43 - a33 * a42;
    float A0323 = a31 * a44 - a34 * a41;
    float A0223 = a31 * a43 - a33 * a41;
    float A0123 = a31 * a42 - a32 * a41;
    float A2313 = a23 * a44 - a24 * a43;
    float A1313 = a22 * a44 - a24 * a42;
    float A1213 = a22 * a43 - a23 * a42;
    float A2312 = a23 * a34 - a24 * a33;
    float A1312 = a22 * a34 - a24 * a32;
    float A1212 = a22 * a33 - a23 * a32;
    float A0313 = a21 * a44 - a24 * a41;
    float A0213 = a21 * a43 - a23 * a41;
    float A0312 = a21 * a34 - a24 * a31;
    float A0212 = a21 * a33 - a23 * a31;
    float A0113 = a21 * a42 - a22 * a41;
    float A0112 = a21 * a32 - a22 * a31;

    INS_INC1(add, 7);
    INS_INC1(mul, 16);

    float det = a11 * (a22 * A2323 - a23 * A1323 + a24 * A1223) -
                a12 * (a21 * A2323 - a23 * A0323 + a24 * A0223) +
                a13 * (a21 * A1323 - a22 * A0323 + a24 * A0123) -
                a14 * (a21 * A1223 - a22 * A0223 + a23 * A0123);
    det = 1 / det;

    m44 inv;

    INS_INC1(add, 40);
    INS_INC1(mul, 64);

    inv.val[0][0] = det * (a22 * A2323 - a23 * A1323 + a24 * A1223);
    inv.val[0][1] = det * -(a12 * A2323 - a13 * A1323 + a14 * A1223);
    inv.val[0][2] = det * (a12 * A2313 - a13 * A1313 + a14 * A1213);
    inv.val[0][3] = det * -(a12 * A2312 - a13 * A1312 + a14 * A1212);
    inv.val[1][0] = det * -(a21 * A2323 - a23 * A0323 + a24 * A0223);
    inv.val[1][1] = det * (a11 * A2323 - a13 * A0323 + a14 * A0223);
    inv.val[1][2] = det * -(a11 * A2313 - a13 * A0313 + a14 * A0213);
    inv.val[1][3] = det * (a11 * A2312 - a13 * A0312 + a14 * A0212);
    inv.val[2][0] = det * (a21 * A1323 - a22 * A0323 + a24 * A0123);
    inv.val[2][1] = det * -(a11 * A1323 - a12 * A0323 + a14 * A0123);
    inv.val[2][2] = det * (a11 * A1313 - a12 * A0313 + a14 * A0113);
    inv.val[2][3] = det * -(a11 * A1312 - a12 * A0312 + a14 * A0112);
    inv.val[3][0] = det * -(a21 * A1223 - a22 * A0223 + a23 * A0123);
    inv.val[3][1] = det * (a11 * A1223 - a12 * A0223 + a13 * A0123);
    inv.val[3][2] = det * -(a11 * A1213 - a12 * A0213 + a13 * A0113);
    inv.val[3][3] = det * (a11 * A1212 - a12 * A0212 + a13 * A0112);
    return inv;
}

m44 get_rot_matrix(vec rot) {
    float c = TO_RAD(rot.x);
    float b = TO_RAD(rot.y);
    float a = TO_RAD(rot.z);

    float ca = cosf(a);
    float cb = cosf(b);
    float cc = cosf(c);
    float sa = sinf(a);
    float sb = sinf(b);
    float sc = sinf(c);

    m44 m = identity;

    INS_INC1(add, 5);
    INS_INC1(mul, 16);

    m.val[0][0] = ca * cb;
    m.val[1][0] = sa * cb;
    m.val[2][0] = -sb;

    m.val[0][1] = ca * sb * sc - sa * cc;
    m.val[1][1] = sa * sb * sc + ca * cc;
    m.val[2][1] = cb * sc;

    m.val[0][2] = ca * sb * cc + sa * sc;
    m.val[1][2] = sa * sb * cc - ca * sc;
    m.val[2][2] = cb * cc;

    return m;
}

/**
 * Creates a object-space to world-space transformation matrix
 *
 * pos is the origin of the object space and the direction of the axis is given
 * as rotation along the three axis in degrees.
 */
m44 get_transf_matrix(vec pos, vec rot) {
    m44 camera_matrix = get_rot_matrix(rot);

    camera_matrix.val[0][3] = pos.x;
    camera_matrix.val[1][3] = pos.y;
    camera_matrix.val[2][3] = pos.z;

    return camera_matrix;
}

// }}}
