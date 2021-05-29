#pragma once

#include <immintrin.h>
#include <math.h>

#include <iomanip>
#include <iostream>

#include "instrument.h"

namespace impl::vec3 {

    extern float M_PI_F;
#define TO_RAD(angle) ((angle) / 180.0f * M_PI_F)

    // 1 / sqrt(3)
#define SQRT3_INV .57735026918962576451f

    static inline float max(float x, float y) {
        INS_CMP;
        return x > y ? x : y;
    }

    static inline float min(float x, float y) {
        INS_CMP;
        return x < y ? x : y;
    }

    static inline float clamp(float x, float lo, float hi) {
        INS_CMP;
        if (x < lo) {
            return lo;
        }

        INS_CMP;
        if (x > hi) {
            return hi;
        }

        return x;
    }

    // Vec {{{

    struct vec {
        float x;
        float y;
        float z;
    };

    struct vec256 {
        __m256 x;
        __m256 y;
        __m256 z;
    };

    std::ostream& operator<<(std::ostream& out, const vec& v);

#define VEC_OP(v1, v2, OP)                                \
    vec {                                                 \
        (v1).x OP(v2).x, (v1).y OP(v2).y, (v1).z OP(v2).z \
    }

    static inline vec vec_add(vec v1, vec v2) {
        INS_INC1(add, 3);
        return VEC_OP(v1, v2, +);
    }

    static inline vec vec_sub(vec v1, vec v2) {
        INS_INC1(add, 3);
        return VEC_OP(v1, v2, -);
    }

    static inline vec vec_mul(vec v1, vec v2) {
        INS_INC1(mul, 3);
        return VEC_OP(v1, v2, *);
    }

    static inline vec vec_div(vec v1, vec v2) {
        INS_INC1(div, 3);
        return VEC_OP(v1, v2, /);
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
        return {max(v.x, val), max(v.y, val), max(v.z, val)};
    }
    // }}}

    static inline vec vec_init(float* values) {
        return {values[0], values[1], values[2]};
    }
    // Vec2 {{{
    struct vec2 {
        float x;
        float y;
    };

#define VEC2_OP(v1, v2, OP)              \
    vec2 {                               \
        (v1).x OP(v2).x, (v1).y OP(v2).y \
    }

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

    static inline vec4 vec4_from_point(vec p) {
        return {p.x, p.y, p.z, 1};
    }

    static inline vec4 vec4_from_dir(vec p) {
        return {p.x, p.y, p.z, 0};
    }

    static inline vec vec4_to_vec(vec4 p) {
        return {p.x, p.y, p.z};
    }

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

    static inline vec get_base_vec(int idx) {
        switch (idx) {
            case 0:
                return {1, 0, 0};
            case 1:
                return {0, 1, 0};
            case 2:
                return {0, 0, 1};
            default:
                return {
                    0,
                    0,
                    0,
                };
        }
    }

    /**    Define a 3x3 matrix
     *
     *    Used for rotations
     */

    struct m33 {
        /** similarly defined as m33
         */
        float val[3][3];
        m33(){};

        m33(vec e1, vec e2, vec e3) {
            val[0][0] = e1.x;
            val[1][0] = e1.y;
            val[2][0] = e1.z;
            val[0][1] = e2.x;
            val[1][1] = e2.y;
            val[2][1] = e2.z;
            val[0][2] = e3.x;
            val[1][2] = e3.y;
            val[2][2] = e3.z;
        };
    };

    static const m33 identity_33 = m33(get_base_vec(0), get_base_vec(1), get_base_vec(2));
    // }}}

    m33 m33_inv(m33 m);
    m33 get_rot_matrix_33(vec rot);

    static inline vec m33_mul_vec(m33 m, vec v) {
        float x = vec_dot(vec_init(m.val[0]), v);
        float y = vec_dot(vec_init(m.val[1]), v);
        float z = vec_dot(vec_init(m.val[2]), v);

        return {x, y, z};
    }


    static inline vec transform_point(const m33 R, const vec t, const vec p) {
        return vec_add(m33_mul_vec(R, p), t);
    }

    static inline vec invtransform_point(const m33 Rt, vec t, const vec p) {
        return m33_mul_vec(Rt, vec_sub(p, t));
    }

    struct Ray {
        vec o; // origin
        vec d; // direction
    };

    static inline vec trace_ray(const Ray r, const float t) {
        return vec_add(r.o, vec_scale(r.d, t));
    }

    static inline vec256 trace_ray_vectorized(int idx, const float* o_x_p, const float* o_y_p, const float* o_z_p, const float* d_x_p, const float* d_y_p, const float* d_z_p, const float t) {
        __m256 o_x = _mm256_loadu_ps(o_x_p + idx);
        __m256 o_y = _mm256_loadu_ps(o_y_p + idx);
        __m256 o_z = _mm256_loadu_ps(o_z_p + idx);
        __m256 d_x = _mm256_loadu_ps(d_x_p + idx);
        __m256 d_y = _mm256_loadu_ps(d_y_p + idx);
        __m256 d_z = _mm256_loadu_ps(d_z_p + idx);

        __m256 t_v = _mm256_set1_ps(t);

        // vec_add(r.o, vec_scale(r.d, t));
        INS_INC1(fma, 24);
        __m256 res_x = _mm256_fmadd_ps(d_x, t_v, o_x);
        __m256 res_y = _mm256_fmadd_ps(d_y, t_v, o_y);
        __m256 res_z = _mm256_fmadd_ps(d_z, t_v, o_z);

        return {res_x, res_y, res_z};
    }

    static inline Ray transform_ray(const m33 R, const vec t, const Ray r) {
        Ray res;
        res.o = transform_point(R, t, r.o);
        res.d = m33_mul_vec(R, r.d);
        return res;
    }

    static inline Ray invtransform_ray(const m33 Rt, const vec t, const Ray r) {
        Ray res;
        res.o = invtransform_point(Rt, t, r.o);
        res.d = m33_mul_vec(Rt, r.d);
        return res;
    }

    // Vectorized operations {{{

    static inline __m256 vectorized_vec_dot(__m256 a_x, __m256 a_y, __m256 a_z, __m256 b_x, __m256 b_y, __m256 b_z) {
        INS_INC1(mul, 8);
        __m256 x = _mm256_mul_ps(a_x, b_x);
        INS_INC1(fma, 16);
        __m256 xy = _mm256_fmadd_ps(a_y, b_y, x);
        __m256 xyz = _mm256_fmadd_ps(a_z, b_z, xy);
        return xyz;
    }

    // }}}

} // namespace impl::vec3
