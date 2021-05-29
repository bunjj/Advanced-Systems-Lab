#pragma once

#include <algorithm>
#include <nlohmann/json.hpp>
#include <immintrin.h>

#include "impl_vec2/geometry.h"

using json = nlohmann::json;

namespace impl::vec2 {

    struct camera {
        float fov;
        vec pos;
        vec rotation;
    };

    struct light {
        vec pos;
        vec emission;
    };

    enum shape_type {
        SHAPE_SPHERE = 0,
        SHAPE_PLANE,
        SHAPE_BOX,
        SHAPE_TORUS,
        SHAPE_CONE,
        SHAPE_OCTA,
    };

    struct sphere {
        vec center;
        float radius;
        vec color;
        float reflection;
        float shininess;
    };

    struct plane {
        // Normal vector
        vec normal;
        // Point on the plane
        vec point;
        vec color;
        float reflection;
        float shininess;
    };

    struct box {
        vec bottom_left;
        vec extents;
        m33 rot;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
        m33 inv_rot;
        float r;
    };

    struct torus {
        vec center;
        float r1;
        float r2;
        float r;
        m33 rot;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
        m33 inv_rot;
    };

    struct cone {
        vec center;
        float r1;
        float r2;
        float height;
        m33 rot;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
        m33 inv_rot;
        float r;
    };

    struct octa {
        vec center;
        float s;
        m33 rot;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
        m33 inv_rot;
    };

    // vectorized data layout
    struct sphere_vectors {
        float* center_x;
        float* center_y;
        float* center_z;
        float* radius;
    };

    struct box_vectors {
        float* bottom_left_x;
        float* bottom_left_y;
        float* bottom_left_z;
        float* extents_x;
        float* extents_y;
        float* extents_z;
        float* r;

        // arrays for each entry of the matrix
        float* inv_rot[3][3];
    };

    struct torus_vectors {
        float* center_x;
        float* center_y;
        float* center_z;
        float* r1;
        float* r2;
        float* r;

        // arrays for each entry of the matrix
        float* inv_rot[3][3];
    };

    struct cone_vectors {
        float* center_x;
        float* center_y;
        float* center_z;
        float* r1;
        float* r2;
        float* height;
        float* r;

        // arrays for each entry of the matrix
        float* inv_rot[3][3];
    };

    struct octa_vectors {
        float* center_x;
        float* center_y;
        float* center_z;
        float* s;

        // arrays for each entry of the matrix
        float* inv_rot[3][3];
    };

    struct scene {
        light* lights;
        int num_lights;
        camera cam;

        int num_spheres;
        int num_planes;
        int num_boxes;
        int num_tori;
        int num_cones;
        int num_octahedra;

        sphere* spheres;
        plane* planes;
        box* boxes;
        torus* tori;
        cone* cones;
        octa* octahedra;

        sphere_vectors sphere_vecs;
        box_vectors box_vecs;
        torus_vectors torus_vecs;
        cone_vectors cone_vecs;
        octa_vectors octa_vecs;
    };

    extern struct scene scene;

    void load_scene(std::string& input);
    void from_ref_scene();

    // distance and normal functions

    // Sphere
    static inline float sphere_distance(const sphere sp, const vec from) {
        INS_INC(sphere);
        INS_ADD;
        return vec_length(vec_sub(sp.center, from)) - sp.radius;
    }

    /**
     * Vectorized sphere distance function (without early termination).
     */
    static inline void sphere_distance_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* radius, const vec from) {
        INS_INC1(sphere, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r = _mm256_loadu_ps(radius + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);

        // vec t = vec_sub(sp.center, from);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(c_x, from_x);
        __m256 t_y = _mm256_sub_ps(c_y, from_y);
        __m256 t_z = _mm256_sub_ps(c_z, from_z);

        // float t_len = vec_len(t);
        __m256 t_dot = vectorized_vec_dot(t_x, t_y, t_z, t_x, t_y, t_z);
        INS_INC1(sqrt, 8);
        __m256 t_len = _mm256_sqrt_ps(t_dot);

        // float res = t_len - sp.radius;
        INS_INC1(add, 8);
        __m256 dist = _mm256_sub_ps(t_len, r);
        _mm256_storeu_ps(res, dist);
    }

    static inline float sphere_distance_short(const sphere sp, const vec from, const float current_min) {
        INS_INC(sphere);
        INS_ADD;
        float upper_bound = current_min + sp.radius;
        float squared_distance = vec_dot2(vec_sub(sp.center, from));
        INS_MUL;
        INS_CMP;
        if (squared_distance >= upper_bound * upper_bound) {
            return current_min;
        }
        INS_ADD;
        return FSQRT(squared_distance) - sp.radius;
    }

    /**
     * Computes the sphere distance function with early termination.
     * Returns zero if early termination is possible, and a nonzero value otherwise.
     */
    static inline int sphere_distance_short_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* radius, const vec from, const float current_min) {
        INS_INC1(sphere, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r = _mm256_loadu_ps(radius + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);
        __m256 curr_min = _mm256_set1_ps(current_min);

        // vec t = vec_sub(sp.center, from);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(c_x, from_x);
        __m256 t_y = _mm256_sub_ps(c_y, from_y);
        __m256 t_z = _mm256_sub_ps(c_z, from_z);

        // float t_len = vec_len(t);
        __m256 t_dot = vectorized_vec_dot(t_x, t_y, t_z, t_x, t_y, t_z);

        // short circuit termination mask
        INS_INC1(add, 8);
        __m256 upper_bound = _mm256_add_ps(curr_min, r);
        INS_INC1(mul, 8);
        __m256 upper_bound_square = _mm256_mul_ps(upper_bound, upper_bound);
        INS_INC1(cmp, 8);
        __m256 mask = _mm256_cmp_ps(t_dot, upper_bound_square, _CMP_LT_OQ);
        int mask_int = _mm256_movemask_ps(mask);

        if (!mask_int) {
            return 0;
        }

        INS_INC1(sqrt, 8);
        __m256 t_len = _mm256_sqrt_ps(t_dot);

        // float res = t_len - sp.radius;
        INS_INC1(add, 8);
        __m256 dist = _mm256_sub_ps(t_len, r);
        _mm256_storeu_ps(res, dist);

        return 1;
    }

    // Box
    static inline float box_distance(const box b, const vec from) {
        INS_INC(box);
        vec pos = m33_mul_vec(b.inv_rot, vec_sub(from, b.bottom_left));
        vec q = vec_sub(vec_abs(pos), b.extents);
        return FADD(vec_length(vec_max(q, 0)), min(0.0f, max(max(q.x, q.y), q.z)));
    }

    /**
     * This is a more explicit version of the box_distance() above.
     * This makes it easier to vectorize.
     */
    static inline float box_distance_elaborate(const box b, const vec from) {
        vec t = vec_sub(from, b.bottom_left);
        vec pos = m33_mul_vec(b.inv_rot, t);
        vec pos_abs = vec_abs(pos);
        vec q = vec_sub(pos_abs, b.extents);
        vec max_q_0 = vec_max(q, 0);
        float left = vec_length(max_q_0);
        float max_qx_qy = max(q.x, q.y);
        float max_qx_qy_qz = max(max_qx_qy, q.z);
        float right = min(0.0f, max_qx_qy_qz);
        return left + right;
    }

    /**
     * Vectorized box distance function (without early termination).
     */
    static inline void box_distance_vectorized(int idx, float* res, float* bottom_left_x, float* bottom_left_y, float* bottom_left_z, float* extents_x, float* extents_y, float* extents_z, float* inv_rot[3][3], const vec from) {
        INS_INC1(box, 8);

        // load everything
        __m256 bl_x = _mm256_loadu_ps(bottom_left_x + idx);
        __m256 bl_y = _mm256_loadu_ps(bottom_left_y + idx);
        __m256 bl_z = _mm256_loadu_ps(bottom_left_z + idx);

        __m256 ext_x = _mm256_loadu_ps(extents_x + idx);
        __m256 ext_y = _mm256_loadu_ps(extents_y + idx);
        __m256 ext_z = _mm256_loadu_ps(extents_z + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);

        // vec t = vec_sub(from, b.bottom_left);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(bl_x, from_x);
        __m256 t_y = _mm256_sub_ps(bl_y, from_y);
        __m256 t_z = _mm256_sub_ps(bl_z, from_z);

        // vec pos = m33_mul_vec(b.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // vec pos_abs = vec_abs(pos);
        // can compute absolute value by setting the sign bits to 0
        INS_INC1(abs, 24);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 pos_abs_x = _mm256_and_ps(abs_mask, pos_x);
        __m256 pos_abs_y = _mm256_and_ps(abs_mask, pos_y);
        __m256 pos_abs_z = _mm256_and_ps(abs_mask, pos_z);

        // vec q = vec_sub(pos_abs, b.extents);
        INS_INC1(add, 24);
        __m256 q_x = _mm256_sub_ps(pos_abs_x, ext_x);
        __m256 q_y = _mm256_sub_ps(pos_abs_y, ext_y);
        __m256 q_z = _mm256_sub_ps(pos_abs_z, ext_z);

        // vec max_q_0 = vec_max(q, 0);
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 24);
        __m256 max_q_0_x = _mm256_max_ps(q_x, zero);
        __m256 max_q_0_y = _mm256_max_ps(q_y, zero);
        __m256 max_q_0_z = _mm256_max_ps(q_z, zero);

        // float left = vec_length(max_q_0);
        __m256 left_square = vectorized_vec_dot(max_q_0_x, max_q_0_y, max_q_0_z, max_q_0_x, max_q_0_y, max_q_0_z);

        INS_INC1(sqrt, 8);
        __m256 left = _mm256_sqrt_ps(left_square);

        // float max_qx_qy = max(q.x, q.y);
        INS_INC1(cmp, 8);
        __m256 max_qx_qy = _mm256_max_ps(q_x, q_y);

        // float max_qx_qy_qz = max(max_qx_qy, q.z);
        INS_INC1(cmp, 8);
        __m256 max_qx_qy_qz = _mm256_max_ps(max_qx_qy, q_z);

        // float right = min(0.0f, max_qx_qy_qz);
        INS_INC1(cmp, 8);
        __m256 right = _mm256_min_ps(zero, max_qx_qy_qz);

        // return left + right;
        INS_INC1(add, 8);
        __m256 dist = _mm256_add_ps(left, right);

        _mm256_storeu_ps(res, dist);
    }

    static inline float box_distance_short(const box b, const vec from, const float current_min) {
        INS_INC(box);
        vec pos = m33_mul_vec(b.inv_rot, vec_sub(from, b.bottom_left));

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = b.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

        INS_INC(box_r);

        vec q = vec_sub(vec_abs(pos), b.extents);
        float extent_values = min(0.0f, max(max(q.x, q.y), q.z));
        float intermediate_squared_dist = vec_dot2(vec_max(q, 0));
        INS_ADD;
        float upper_bound = extent_values + current_min;
        INS_MUL;
        INS_CMP;
        if (intermediate_squared_dist >= upper_bound * upper_bound) {
            return current_min;
        }
        return FADD(FSQRT(intermediate_squared_dist), extent_values);
    }

    /**
     * Computes the box distance function with early termination.
     * Returns zero if early termination is possible, and a nonzero value otherwise.
     */
    static inline int box_distance_short_vectorized(int idx, float* res, float* bottom_left_x, float* bottom_left_y, float* bottom_left_z, float* extents_x, float* extents_y, float* extents_z, float* rad, float* inv_rot[3][3], const vec from, float current_min) {
        INS_INC1(box, 8);

        // load vectors
        __m256 bl_x = _mm256_loadu_ps(bottom_left_x + idx);
        __m256 bl_y = _mm256_loadu_ps(bottom_left_y + idx);
        __m256 bl_z = _mm256_loadu_ps(bottom_left_z + idx);

        __m256 r = _mm256_loadu_ps(rad + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);
        __m256 curr_min = _mm256_set1_ps(current_min);

        // vec t = vec_sub(from, b.bottom_left);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(bl_x, from_x);
        __m256 t_y = _mm256_sub_ps(bl_y, from_y);
        __m256 t_z = _mm256_sub_ps(bl_z, from_z);

        // vec pos = m33_mul_vec(b.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // enclosing sphere check
        __m256 pos_square = vectorized_vec_dot(pos_x, pos_y, pos_z, pos_x, pos_y, pos_z);
        INS_INC1(add, 8);
        __m256 duf_bound = _mm256_add_ps(r, curr_min);
        INS_INC1(mul, 8);
        __m256 duf_bound_square = _mm256_mul_ps(duf_bound, duf_bound);
        INS_INC1(cmp, 8);
        __m256 duf_mask = _mm256_cmp_ps(pos_square, duf_bound_square, _CMP_LT_OQ);
        int duf_mask_int = _mm256_movemask_ps(duf_mask);
        if (!duf_mask_int) {
            return 0;
        }

        INS_INC1(box_r, 8);

        // load remaining vectors
        __m256 ext_x = _mm256_loadu_ps(extents_x + idx);
        __m256 ext_y = _mm256_loadu_ps(extents_y + idx);
        __m256 ext_z = _mm256_loadu_ps(extents_z + idx);

        // vec pos_abs = vec_abs(pos);
        // can compute absolute value by setting the sign bits to 0
        INS_INC1(abs, 24);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 pos_abs_x = _mm256_and_ps(abs_mask, pos_x);
        __m256 pos_abs_y = _mm256_and_ps(abs_mask, pos_y);
        __m256 pos_abs_z = _mm256_and_ps(abs_mask, pos_z);

        // vec q = vec_sub(pos_abs, b.extents);
        INS_INC1(add, 24);
        __m256 q_x = _mm256_sub_ps(pos_abs_x, ext_x);
        __m256 q_y = _mm256_sub_ps(pos_abs_y, ext_y);
        __m256 q_z = _mm256_sub_ps(pos_abs_z, ext_z);

        // vec max_q_0 = vec_max(q, 0);
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 24);
        __m256 max_q_0_x = _mm256_max_ps(q_x, zero);
        __m256 max_q_0_y = _mm256_max_ps(q_y, zero);
        __m256 max_q_0_z = _mm256_max_ps(q_z, zero);

        // float left = vec_length(max_q_0);
        __m256 left_square = vectorized_vec_dot(max_q_0_x, max_q_0_y, max_q_0_z, max_q_0_x, max_q_0_y, max_q_0_z);

        // float max_qx_qy = max(q.x, q.y);
        INS_INC1(cmp, 8);
        __m256 max_qx_qy = _mm256_max_ps(q_x, q_y);

        // float max_qx_qy_qz = max(max_qx_qy, q.z);
        INS_INC1(cmp, 8);
        __m256 max_qx_qy_qz = _mm256_max_ps(max_qx_qy, q_z);

        // float right = min(0.0f, max_qx_qy_qz);
        INS_INC1(cmp, 8);
        __m256 right = _mm256_min_ps(zero, max_qx_qy_qz);

        // short circuit termination mask
        INS_INC1(add, 8);
        __m256 upper_bound = _mm256_add_ps(curr_min, right);
        INS_INC1(mul, 8);
        __m256 upper_bound_square = _mm256_mul_ps(upper_bound, upper_bound);
        INS_INC1(cmp, 8);
        __m256 mask = _mm256_cmp_ps(left_square, upper_bound_square, _CMP_LT_OQ);
        int mask_int = _mm256_movemask_ps(mask);

        if (!mask_int) {
            return 0;
        }

        INS_INC1(sqrt, 8);
        __m256 left = _mm256_sqrt_ps(left_square);

        // return left + right;
        INS_INC1(add, 8);
        __m256 dist = _mm256_add_ps(left, right);

        _mm256_store_ps(res, dist);

        return 1;
    }

    // Plane
    static inline float plane_distance(const plane p, const vec from) {
        INS_INC(plane);
        return vec_dot(p.normal, vec_sub(from, p.point));
    }

    // Torus
    static inline float torus_distance(const torus t, const vec from) {
        INS_INC(torus);
        vec pos = m33_mul_vec(t.inv_rot, vec_sub(from, t.center));
        vec2 posxz = {pos.x, pos.z};
        INS_ADD;
        vec2 q = {vec2_length(posxz) - t.r1, pos.y};

        INS_ADD;
        return vec2_length(q) - t.r2;
    }

    /**
     * This is a more explicit version of the torus_distance() above.
     * This makes it easier to vectorize.
     */
    static inline float torus_distance_elaborate(const torus to, const vec from) {
        vec t = vec_sub(from, to.center);
        vec pos = m33_mul_vec(to.inv_rot, t);
        vec2 posxz = {pos.x, pos.z};
        float posxz_len = vec2_length(posxz);
        float q1 = posxz_len - to.r1;
        vec2 q = {q1, pos.y};
        float q_len = vec2_length(q);
        return q_len - to.r2;
    }

    /**
     * Vectorized torus distance function (without early termination).
     */
    static inline void torus_distance_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* rad1, float* rad2, float* inv_rot[3][3], const vec from) {
        INS_INC1(torus, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r1 = _mm256_loadu_ps(rad1 + idx);
        __m256 r2 = _mm256_loadu_ps(rad2 + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);

        // vec t = vec_sub(from, to.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(to.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // vec2 posxz = {pos.x, pos.z};
        // float posxz_len = vec2_length(posxz);
        INS_INC1(mul, 8);
        __m256 pos_x_square = _mm256_mul_ps(pos_x, pos_x);
        INS_INC1(fma, 8);
        __m256 pos_xz_square = _mm256_fmadd_ps(pos_z, pos_z, pos_x_square);

        INS_INC1(sqrt, 8);
        __m256 pos_xz_len = _mm256_sqrt_ps(pos_xz_square);

        // float q1 = posxz_len - to.r1;
        INS_INC1(add, 8);
        __m256 q1 = _mm256_sub_ps(pos_xz_len, r1);

        // vec2 q = {q1, pos.y};
        // float q_len = vec2_length(q);
        INS_INC1(mul, 8);
        __m256 q1_square = _mm256_mul_ps(q1, q1);
        INS_INC1(fma, 8);
        __m256 q_square = _mm256_fmadd_ps(pos_y, pos_y, q1_square);

        INS_INC1(sqrt, 8);
        __m256 q_len = _mm256_sqrt_ps(q_square);

        // return q_len - to.r2;
        INS_INC1(add, 8);
        __m256 dist = _mm256_sub_ps(q_len, r2);

        _mm256_storeu_ps(res, dist);
    }

    static inline float torus_distance_short(const torus t, const vec from, const float current_min) {
        INS_INC(torus);
        vec pos = m33_mul_vec(t.inv_rot, vec_sub(from, t.center));

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = t.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;
        INS_INC(torus_r);
        vec2 posxz = {pos.x, pos.z};
        INS_ADD;
        vec2 q = {vec2_length(posxz) - t.r1, pos.y};
        float q_squared = vec2_dot2(q);

        INS_ADD;
        float upper_bound = current_min + t.r2;
        INS_MUL;
        INS_CMP;
        if (q_squared >= upper_bound * upper_bound) {
            return current_min;
        }

        INS_ADD;
        return FSQRT(q_squared) - t.r2;
    }

    /**
     * Computes the torus distance function with early termination.
     * Returns zero if early termination is possible, and a nonzero value otherwise.
     */
    static inline int torus_distance_short_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* rad1, float* rad2, float* rad, float* inv_rot[3][3], const vec from, float current_min) {
        INS_INC1(torus, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r = _mm256_loadu_ps(rad + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);
        __m256 curr_min = _mm256_set1_ps(current_min);

        // vec t = vec_sub(from, to.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(to.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // enclosing sphere check
        INS_INC1(mul, 8);
        __m256 pos_x_square = _mm256_mul_ps(pos_x, pos_x);
        INS_INC1(fma, 16);
        __m256 pos_xz_square = _mm256_fmadd_ps(pos_z, pos_z, pos_x_square);
        __m256 pos_square = _mm256_fmadd_ps(pos_y, pos_y, pos_xz_square);
        INS_INC1(add, 8);
        __m256 duf_bound = _mm256_add_ps(r, curr_min);
        INS_INC1(mul, 8);
        __m256 duf_bound_square = _mm256_mul_ps(duf_bound, duf_bound);
        INS_INC1(cmp, 8);
        __m256 duf_mask = _mm256_cmp_ps(pos_square, duf_bound_square, _CMP_LT_OQ);
        int duf_mask_int = _mm256_movemask_ps(duf_mask);
        if (!duf_mask_int) {
            return 0;
        }

        INS_INC1(torus_r, 8);

        // vec2 posxz = {pos.x, pos.z};
        // float posxz_len = vec2_length(posxz);
        INS_INC1(sqrt, 8);
        __m256 pos_xz_len = _mm256_sqrt_ps(pos_xz_square);

        // load remaining vectors
        __m256 r1 = _mm256_loadu_ps(rad1 + idx);
        __m256 r2 = _mm256_loadu_ps(rad2 + idx);

        // float q1 = posxz_len - to.r1;
        INS_INC1(add, 8);
        __m256 q1 = _mm256_sub_ps(pos_xz_len, r1);

        // vec2 q = {q1, pos.y};
        // float q_len = vec2_length(q);
        INS_INC1(mul, 8);
        __m256 q1_square = _mm256_mul_ps(q1, q1);
        INS_INC1(fma, 8);
        __m256 q_square = _mm256_fmadd_ps(pos_y, pos_y, q1_square);

        // short circuit termination mask
        INS_INC1(add, 8);
        __m256 upper_bound = _mm256_add_ps(curr_min, r2);
        INS_INC1(mul, 8);
        __m256 upper_bound_square = _mm256_mul_ps(upper_bound, upper_bound);
        INS_INC1(cmp, 8);
        __m256 mask = _mm256_cmp_ps(q_square, upper_bound_square, _CMP_LT_OQ);
        int mask_int = _mm256_movemask_ps(mask);

        if (!mask_int) {
            return 0;
        }

        INS_INC1(sqrt, 8);
        __m256 q_len = _mm256_sqrt_ps(q_square);

        // return q_len - to.r2;
        INS_INC1(add, 8);
        __m256 dist = _mm256_sub_ps(q_len, r2);

        _mm256_store_ps(res, dist);

        return 1;
    }

    // Cone
    static inline float cone_distance(const cone c, const vec from) {
        INS_INC(cone);
        vec pos = m33_mul_vec(c.inv_rot, vec_sub(from, c.center));

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        vec2 q = {vec2_length({pos.x, pos.z}), pos.y};
        vec2 k1 = {r2, h};
        INS_ADD;
        INS_MUL;
        vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 2);
        INS_ABS;
        INS_CMP;
        vec2 ca = {q.x - min(q.x, (q.y < 0 ? r1 : r2)), fabsf(q.y) - h};
        INS_DIV;
        vec2 cb =
            vec2_add(vec2_sub(q, k1), vec2_scale(k2, clamp(vec2_dot(vec2_sub(k1, q), k2) / vec2_dot2(k2), 0.0f, 1.0f)));
        INS_INC1(cmp, 2);
        float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;

        INS_MUL;
        INS_SQRT;
        return s * sqrtf(min(vec2_dot2(ca), vec2_dot2(cb)));
    }

    /**
     * This is a more explicit version of the cone_distance() above.
     * This makes it easier to vectorize.
     */
    static inline float cone_distance_elaborate(const cone c, const vec from) {
        vec t = vec_sub(from, c.center);
        vec pos = m33_mul_vec(c.inv_rot, t);

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        float xz_len = vec2_length({pos.x, pos.z});
        vec2 q = {xz_len, pos.y};
        vec2 k1 = {r2, h};
        vec2 k2 = {r2 - r1, 2 * h};

        float r1_or_r2 = q.y < 0 ? r1 : r2;
        float min_qx_r1_or_r2 = min(q.x, r1_or_r2);
        float ca1 = q.x - min_qx_r1_or_r2;
        float qy_abs = fabsf(q.y);
        float ca2 = qy_abs - h;
        vec2 ca = {ca1, ca2};

        vec2 left = vec2_sub(q, k1);

        float k2_dot2 = vec2_dot2(k2);
        vec2 k1_minus_q = vec2_sub(k1, q);
        float k1_minus_q_dot_k2 = vec2_dot(k1_minus_q, k2);
        float to_clamp = k1_minus_q_dot_k2 / k2_dot2;
        float clamped = clamp(to_clamp, 0.0f, 1.0f);
        vec2 right = vec2_scale(k2, clamped);

        vec2 cb = vec2_add(left, right);
        float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;

        float ca_dot2 = vec2_dot2(ca);
        float cb_dot2 = vec2_dot2(cb);
        float min_square = min(ca_dot2, cb_dot2);
        float min = sqrtf(min_square);

        return s * min;
    }

    /**
     * Vectorized cone distance function (without early termination).
     */
    static inline void cone_distance_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* rad1, float* rad2, float* height, float* inv_rot[3][3], const vec from) {
        INS_INC1(cone, 8);

        // float r1 = c.r1;
        // float r2 = c.r2;
        // float h = c.height;

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r1 = _mm256_loadu_ps(rad1 + idx);
        __m256 r2 = _mm256_loadu_ps(rad2 + idx);
        __m256 h = _mm256_loadu_ps(height + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);

        // vec t = vec_sub(from, c.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(c.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // float xz_len = vec2_length({pos.x, pos.z});
        INS_INC1(mul, 8);
        __m256 pos_x_square = _mm256_mul_ps(pos_x, pos_x);
        INS_INC1(fma, 8);
        __m256 pos_xz_square = _mm256_fmadd_ps(pos_z, pos_z, pos_x_square);
        INS_INC1(sqrt, 8);
        __m256 xz_len = _mm256_sqrt_ps(pos_xz_square);

        // vec2 q = {xz_len, pos.y};
        // vec2 k1 = {r2, h};

        // float r1_or_r2 = q.y < 0 ? r1 : r2;
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 8);
        __m256 qy_lt_0_mask = _mm256_cmp_ps(pos_y, zero, _CMP_LT_OQ);
        __m256 r1_or_r2 = _mm256_blendv_ps(r2, r1, qy_lt_0_mask);

        // float min_qx_r1_or_r2 = min(q.x, r1_or_r2);
        INS_INC1(cmp, 8);
        __m256 min_qx_r1_or_r2 = _mm256_min_ps(xz_len, r1_or_r2);

        // float ca1 = q.x - min_qx_r1_or_r2;
        INS_INC1(add, 8);
        __m256 ca1 = _mm256_sub_ps(xz_len, min_qx_r1_or_r2);

        // float qy_abs = fabsf(q.y);
        INS_INC1(abs, 8);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 qy_abs = _mm256_and_ps(abs_mask, pos_y);

        // float ca2 = qy_abs - h;
        INS_INC1(add, 8);
        __m256 ca2 = _mm256_sub_ps(qy_abs, h);

        // vec2 ca = {ca1, ca2};

        // vec2 left = vec2_sub(q, k1);
        INS_INC1(add, 16);
        __m256 left1 = _mm256_sub_ps(xz_len, r2);
        __m256 left2 = _mm256_sub_ps(pos_y, h);

        // vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 8);
        __m256 r2_minus_r1 = _mm256_sub_ps(r2, r1);
        __m256 two = _mm256_set1_ps(2.f);
        INS_INC1(mul, 8);
        __m256 h_twice = _mm256_mul_ps(h, two);

        // float k2_dot2 = vec2_dot2(k2);
        INS_INC1(mul, 8);
        __m256 r2_minus_r1_square = _mm256_mul_ps(r2_minus_r1, r2_minus_r1);
        INS_INC1(fma, 8);
        __m256 k2_dot2 = _mm256_fmadd_ps(h_twice, h_twice, r2_minus_r1_square);

        // vec2 k1_minus_q = vec2_sub(k1, q);
        INS_INC1(add, 16);
        __m256 k1_minus_q_1 = _mm256_sub_ps(r2, xz_len);
        __m256 k1_minus_q_2 = _mm256_sub_ps(h, pos_y);

        // float k1_minus_q_dot_k2 = vec2_dot(k1_minus_q, k2);
        INS_INC1(mul, 8);
        __m256 k1_minus_q_dot_k2_1 = _mm256_mul_ps(k1_minus_q_1, r2_minus_r1);
        INS_INC1(fma, 8);
        __m256 k1_minus_q_dot_k2 = _mm256_fmadd_ps(k1_minus_q_2, h_twice, k1_minus_q_dot_k2_1);

        // float to_clamp = k1_minus_q_dot_k2 / k2_dot2;
        INS_INC1(div, 8);
        __m256 to_clamp = _mm256_div_ps(k1_minus_q_dot_k2, k2_dot2);

        // float clamped = clamp(to_clamp, 0.0f, 1.0f);
        __m256 one = _mm256_set1_ps(1.f);
        INS_INC1(cmp, 16);
        __m256 clamped_upper = _mm256_min_ps(to_clamp, one);
        __m256 clamped = _mm256_max_ps(clamped_upper, zero);

        // vec2 right = vec2_scale(k2, clamped);
        INS_INC1(mul, 16);
        __m256 right1 = _mm256_mul_ps(r2_minus_r1, clamped);
        __m256 right2 = _mm256_mul_ps(h_twice, clamped);

        // vec2 cb = vec2_add(left, right);
        INS_INC1(add, 16);
        __m256 cb1 = _mm256_add_ps(left1, right1);
        __m256 cb2 = _mm256_add_ps(left2, right2);

        // float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;
        INS_INC1(cmp, 16);
        __m256 cb1_lt_0_mask = _mm256_cmp_ps(cb1, zero, _CMP_LT_OQ);
        __m256 ca2_lt_0_mask = _mm256_cmp_ps(ca2, zero, _CMP_LT_OQ);
        __m256 cond_mask = _mm256_and_ps(cb1_lt_0_mask, ca2_lt_0_mask);
        __m256 minusone = _mm256_set1_ps(-1.f);
        __m256 s = _mm256_blendv_ps(one, minusone, cond_mask);

        // float ca_dot2 = vec2_dot2(ca);
        INS_INC1(mul, 8);
        __m256 ca1_square = _mm256_mul_ps(ca1, ca1);
        INS_INC1(fma, 8);
        __m256 ca_dot2 = _mm256_fmadd_ps(ca2, ca2, ca1_square);

        // float cb_dot2 = vec2_dot2(cb);
        INS_INC1(mul, 8);
        __m256 cb1_square = _mm256_mul_ps(cb1, cb1);
        INS_INC1(fma, 8);
        __m256 cb_dot2 = _mm256_fmadd_ps(cb2, cb2, cb1_square);

        // float min_square = min(ca_dot2, cb_dot2);
        INS_INC1(cmp, 8);
        __m256 min_square = _mm256_min_ps(ca_dot2, cb_dot2);

        // float min = sqrtf(min_square);
        INS_INC1(sqrt, 8);
        __m256 min = _mm256_sqrt_ps(min_square);

        // return s * min;
        INS_INC1(mul, 8);
        __m256 dist = _mm256_mul_ps(s, min);

        _mm256_storeu_ps(res, dist);
    }

    static inline float cone_distance_short(const cone c, const vec from, const float current_min) {
        INS_INC(cone);
        vec pos = m33_mul_vec(c.inv_rot, vec_sub(from, c.center));

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = c.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

        INS_INC(cone_r);

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        vec2 q = {vec2_length({pos.x, pos.z}), pos.y};
        vec2 k1 = {r2, h};
        INS_ADD;
        INS_MUL;
        vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 2);
        INS_ABS;
        INS_CMP;
        vec2 ca = {q.x - min(q.x, (q.y < 0 ? r1 : r2)), fabsf(q.y) - h};
        INS_DIV;
        vec2 cb =
            vec2_add(vec2_sub(q, k1), vec2_scale(k2, clamp(vec2_dot(vec2_sub(k1, q), k2) / vec2_dot2(k2), 0.0f, 1.0f)));
        INS_INC1(cmp, 2);
        float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;

        float squared_min = min(vec2_dot2(ca), vec2_dot2(cb));
        INS_MUL;
        INS_CMP;
        if (squared_min >= current_min * current_min) {
            return current_min;
        }

        INS_MUL;
        INS_SQRT;
        return s * sqrtf(squared_min);
    }

    /**
     * Computes the cone distance function with early termination.
     * Returns zero if early termination is possible, and a nonzero value otherwise.
     */
    static inline int cone_distance_short_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* rad1, float* rad2, float* height, float* rad, float* inv_rot[3][3], const vec from, float current_min) {
        INS_INC1(cone, 8);

        // float r1 = c.r1;
        // float r2 = c.r2;
        // float h = c.height;

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 r = _mm256_loadu_ps(rad + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);
        __m256 curr_min = _mm256_set1_ps(current_min);

        // vec t = vec_sub(from, c.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(c.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // enclosing sphere check
        INS_INC1(mul, 8);
        __m256 pos_x_square = _mm256_mul_ps(pos_x, pos_x);
        INS_INC1(fma, 16);
        __m256 pos_xz_square = _mm256_fmadd_ps(pos_z, pos_z, pos_x_square);
        __m256 pos_square = _mm256_fmadd_ps(pos_y, pos_y, pos_xz_square);
        INS_INC1(add, 8);
        __m256 duf_bound = _mm256_add_ps(r, curr_min);
        INS_INC1(mul, 8);
        __m256 duf_bound_square = _mm256_mul_ps(duf_bound, duf_bound);
        INS_INC1(cmp, 8);
        __m256 duf_mask = _mm256_cmp_ps(pos_square, duf_bound_square, _CMP_LT_OQ);
        int duf_mask_int = _mm256_movemask_ps(duf_mask);
        if (!duf_mask_int) {
            return 0;
        }

        INS_INC1(cone_r, 8);

        // load remaining vectors
        __m256 r1 = _mm256_loadu_ps(rad1 + idx);
        __m256 r2 = _mm256_loadu_ps(rad2 + idx);
        __m256 h = _mm256_loadu_ps(height + idx);

        // float xz_len = vec2_length({pos.x, pos.z});
        INS_INC1(sqrt, 8);
        __m256 xz_len = _mm256_sqrt_ps(pos_xz_square);

        // vec2 q = {xz_len, pos.y};
        // vec2 k1 = {r2, h};

        // float r1_or_r2 = q.y < 0 ? r1 : r2;
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 8);
        __m256 qy_lt_0_mask = _mm256_cmp_ps(pos_y, zero, _CMP_LT_OQ);
        __m256 r1_or_r2 = _mm256_blendv_ps(r2, r1, qy_lt_0_mask);

        // float min_qx_r1_or_r2 = min(q.x, r1_or_r2);
        INS_INC1(cmp, 8);
        __m256 min_qx_r1_or_r2 = _mm256_min_ps(xz_len, r1_or_r2);

        // float ca1 = q.x - min_qx_r1_or_r2;
        INS_INC1(add, 8);
        __m256 ca1 = _mm256_sub_ps(xz_len, min_qx_r1_or_r2);

        // float qy_abs = fabsf(q.y);
        INS_INC1(abs, 8);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 qy_abs = _mm256_and_ps(abs_mask, pos_y);

        // float ca2 = qy_abs - h;
        INS_INC1(add, 8);
        __m256 ca2 = _mm256_sub_ps(qy_abs, h);

        // vec2 ca = {ca1, ca2};

        // vec2 left = vec2_sub(q, k1);
        INS_INC1(add, 16);
        __m256 left1 = _mm256_sub_ps(xz_len, r2);
        __m256 left2 = _mm256_sub_ps(pos_y, h);

        // vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 8);
        __m256 r2_minus_r1 = _mm256_sub_ps(r2, r1);
        __m256 two = _mm256_set1_ps(2.f);
        INS_INC1(mul, 8);
        __m256 h_twice = _mm256_mul_ps(h, two);

        // float k2_dot2 = vec2_dot2(k2);
        INS_INC1(mul, 8);
        __m256 r2_minus_r1_square = _mm256_mul_ps(r2_minus_r1, r2_minus_r1);
        INS_INC1(fma, 8);
        __m256 k2_dot2 = _mm256_fmadd_ps(h_twice, h_twice, r2_minus_r1_square);

        // vec2 k1_minus_q = vec2_sub(k1, q);
        INS_INC1(add, 16);
        __m256 k1_minus_q_1 = _mm256_sub_ps(r2, xz_len);
        __m256 k1_minus_q_2 = _mm256_sub_ps(h, pos_y);

        // float k1_minus_q_dot_k2 = vec2_dot(k1_minus_q, k2);
        INS_INC1(mul, 8);
        __m256 k1_minus_q_dot_k2_1 = _mm256_mul_ps(k1_minus_q_1, r2_minus_r1);
        INS_INC1(fma, 8);
        __m256 k1_minus_q_dot_k2 = _mm256_fmadd_ps(k1_minus_q_2, h_twice, k1_minus_q_dot_k2_1);

        // float to_clamp = k1_minus_q_dot_k2 / k2_dot2;
        INS_INC1(div, 8);
        __m256 to_clamp = _mm256_div_ps(k1_minus_q_dot_k2, k2_dot2);

        // float clamped = clamp(to_clamp, 0.0f, 1.0f);
        __m256 one = _mm256_set1_ps(1.f);
        INS_INC1(cmp, 16);
        __m256 clamped_upper = _mm256_min_ps(to_clamp, one);
        __m256 clamped = _mm256_max_ps(clamped_upper, zero);

        // vec2 right = vec2_scale(k2, clamped);
        INS_INC1(mul, 16);
        __m256 right1 = _mm256_mul_ps(r2_minus_r1, clamped);
        __m256 right2 = _mm256_mul_ps(h_twice, clamped);

        // vec2 cb = vec2_add(left, right);
        INS_INC1(add, 16);
        __m256 cb1 = _mm256_add_ps(left1, right1);
        __m256 cb2 = _mm256_add_ps(left2, right2);

        // float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;
        INS_INC1(cmp, 16);
        __m256 cb1_lt_0_mask = _mm256_cmp_ps(cb1, zero, _CMP_LT_OQ);
        __m256 ca2_lt_0_mask = _mm256_cmp_ps(ca2, zero, _CMP_LT_OQ);
        __m256 cond_mask = _mm256_and_ps(cb1_lt_0_mask, ca2_lt_0_mask);
        __m256 minusone = _mm256_set1_ps(-1.f);
        __m256 s = _mm256_blendv_ps(one, minusone, cond_mask);

        // float ca_dot2 = vec2_dot2(ca);
        INS_INC1(mul, 8);
        __m256 ca1_square = _mm256_mul_ps(ca1, ca1);
        INS_INC1(fma, 8);
        __m256 ca_dot2 = _mm256_fmadd_ps(ca2, ca2, ca1_square);

        // float cb_dot2 = vec2_dot2(cb);
        INS_INC1(mul, 8);
        __m256 cb1_square = _mm256_mul_ps(cb1, cb1);
        INS_INC1(fma, 8);
        __m256 cb_dot2 = _mm256_fmadd_ps(cb2, cb2, cb1_square);

        // float min_square = min(ca_dot2, cb_dot2);
        INS_INC1(cmp, 8);
        __m256 min_square = _mm256_min_ps(ca_dot2, cb_dot2);

        // short circuit termination mask
        INS_INC1(mul, 8);
        __m256 upper_bound_square = _mm256_mul_ps(curr_min, curr_min);
        INS_INC1(cmp, 8);
        __m256 mask = _mm256_cmp_ps(min_square, upper_bound_square, _CMP_LT_OQ);
        int mask_int = _mm256_movemask_ps(mask);

        if (!mask_int) {
            return 0;
        }

        // float min = sqrtf(min_square);
        INS_INC1(sqrt, 8);
        __m256 min = _mm256_sqrt_ps(min_square);

        // return s * min;
        INS_INC1(mul, 8);
        __m256 dist = _mm256_mul_ps(s, min);

        _mm256_storeu_ps(res, dist);

        return 1;
    }

    // Octahedron
    static inline float octahedron_distance(const octa o, const vec from) {
        INS_INC(octa);
        vec pos = m33_mul_vec(o.inv_rot, vec_sub(from, o.center));
        pos = vec_abs(pos);

        float s = o.s;

        INS_INC1(add, 3);
        float m = pos.x + pos.y + pos.z - s;
        vec q;

        if (3 * pos.x < m) {
            INS_INC1(mul, 1);
            INS_INC1(cmp, 1);
            q = pos;
        } else if (3 * pos.y < m) {
            INS_INC1(mul, 2);
            INS_INC1(cmp, 2);
            q = {pos.y, pos.x, pos.z};
        } else if (3 * pos.z < m) {
            INS_INC1(mul, 3);
            INS_INC1(cmp, 3);
            q = {pos.z, pos.x, pos.y};
        } else {
            INS_INC1(mul, 4);
            INS_INC1(cmp, 3);
            return m * 0.57735027;
        }

        INS_MUL;
        INS_INC1(add, 2);
        float k = clamp(0.5f * (q.z - q.y + s), 0.0f, s);

        INS_INC1(add, 3);
        return vec_length({q.x, q.y - s + k, q.z - k});
    }

    /**
     * This is a more explicit version of the octahedron_distance() above.
     * This makes it easier to vectorize.
     */
    static inline float octahedron_distance_elaborate(const octa o, const vec from) {
        vec t = vec_sub(from, o.center);
        vec pos = m33_mul_vec(o.inv_rot, t);
        vec pos_abs = vec_abs(pos);

        float s = o.s;

        float m = pos_abs.x + pos_abs.y + pos_abs.z - s;
        vec q;

        float pos_abs_x_times3 = 3 * pos_abs.x;
        float pos_abs_y_times3 = 3 * pos_abs.y;
        float pos_abs_z_times3 = 3 * pos_abs.z;

        int cond1 = pos_abs_x_times3 < m;
        int cond2 = pos_abs_y_times3 < m;
        int cond3 = pos_abs_z_times3 < m;

        if (cond1) {
            q = pos_abs;
        } else if (cond2) {
            q = {pos_abs.y, pos_abs.x, pos_abs.z};
        } else {
            q = {pos_abs.z, pos_abs.x, pos_abs.y};
        }

        float qz_minus_qy = q.z - q.y;
        float qz_minus_qy_plus_s = qz_minus_qy + s;
        float to_clamp = 0.5f * qz_minus_qy_plus_s;
        float k = clamp(to_clamp, 0.0f, s);

        vec dist_vec = {q.x, q.y - s + k, q.z - k};
        float dist_vec_len = vec_length(dist_vec);

        float m_scaled = m * 0.57735027;
        return (!cond1 && !cond2 && !cond3) ? m_scaled : dist_vec_len;
    }

    static inline void octahedron_distance_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* s_, float* inv_rot[3][3], const vec from) {
        INS_INC1(octa, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 s = _mm256_loadu_ps(s_ + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);

        // vec t = vec_sub(from, o.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(o.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // vec pos_abs = vec_abs(pos);
        INS_INC1(abs, 24);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 pos_abs_x = _mm256_and_ps(abs_mask, pos_x);
        __m256 pos_abs_y = _mm256_and_ps(abs_mask, pos_y);
        __m256 pos_abs_z = _mm256_and_ps(abs_mask, pos_z);

        // float s = o.s;

        // float m = pos_abs.x + pos_abs.y + pos_abs.z - s;
        INS_INC1(add, 24);
        __m256 pos_abs_xy = _mm256_add_ps(pos_abs_x, pos_abs_y);
        __m256 pos_abs_xyz = _mm256_add_ps(pos_abs_xy, pos_abs_z);
        __m256 m = _mm256_sub_ps(pos_abs_xyz, s);

        // vec q;

        // float pos_abs_x_times3 = 3 * pos_abs.x;
        // float pos_abs_y_times3 = 3 * pos_abs.y;
        // float pos_abs_z_times3 = 3 * pos_abs.z;
        __m256 three = _mm256_set1_ps(3.f);
        INS_INC1(mul, 24);
        __m256 pos_abs_x_times3 = _mm256_mul_ps(three, pos_abs_x);
        __m256 pos_abs_y_times3 = _mm256_mul_ps(three, pos_abs_y);
        __m256 pos_abs_z_times3 = _mm256_mul_ps(three, pos_abs_z);

        // int cond1 = pos_abs_x_times3 < m;
        // int cond2 = pos_abs_y_times3 < m;
        // int cond3 = pos_abs_z_times3 < m;
        INS_INC1(cmp, 24);
        __m256 cond1 = _mm256_cmp_ps(pos_abs_x_times3, m, _CMP_LT_OQ);
        __m256 cond2 = _mm256_cmp_ps(pos_abs_y_times3, m, _CMP_LT_OQ);
        __m256 cond3 = _mm256_cmp_ps(pos_abs_z_times3, m, _CMP_LT_OQ);

        // if (cond1) {
        //     q = pos_abs;
        // } else if (cond2) {
        //     q = {pos_abs.y, pos_abs.x, pos_abs.z};
        // } else {
        //     q = {pos_abs.z, pos_abs.x, pos_abs.y};
        // }

        __m256 q_branch_23_x = _mm256_blendv_ps(pos_abs_z, pos_abs_y, cond2);
        __m256 q_branch_23_y = pos_abs_x;
        __m256 q_branch_23_z = _mm256_blendv_ps(pos_abs_y, pos_abs_z, cond2);

        __m256 q_x = _mm256_blendv_ps(q_branch_23_x, pos_abs_x, cond1);
        __m256 q_y = _mm256_blendv_ps(q_branch_23_y, pos_abs_y, cond1);
        __m256 q_z = _mm256_blendv_ps(q_branch_23_z, pos_abs_z, cond1);

        // float qz_minus_qy = q.z - q.y;
        INS_INC1(add, 8);
        __m256 qz_minus_qy = _mm256_sub_ps(q_z, q_y);

        // float qz_minus_qy_plus_s = qz_minus_qy + s;
        INS_INC1(add, 8);
        __m256 qz_minus_qy_plus_s = _mm256_add_ps(qz_minus_qy, s);

        // float to_clamp = 0.5f * qz_minus_qy_plus_s;
        __m256 pointfive = _mm256_set1_ps(0.5f);
        INS_INC1(mul, 8);
        __m256 to_clamp = _mm256_mul_ps(pointfive, qz_minus_qy_plus_s);

        // float k = clamp(to_clamp, 0.0f, s);
        INS_INC1(cmp, 8);
        __m256 clamped_upper = _mm256_min_ps(to_clamp, s);
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 8);
        __m256 k = _mm256_max_ps(zero, clamped_upper);

        // vec dist_vec = {q.x, q.y - s + k, q.z - k};
        __m256 dist_vec_x = q_x;
        INS_INC1(add, 24);
        __m256 qy_minus_s = _mm256_sub_ps(q_y, s);
        __m256 dist_vec_y = _mm256_add_ps(qy_minus_s, k);
        __m256 dist_vec_z = _mm256_sub_ps(q_z, k);

        // float dist_vec_len = vec_length(dist_vec);
        __m256 dist_vec_dot2 = vectorized_vec_dot(dist_vec_x, dist_vec_y, dist_vec_z, dist_vec_x, dist_vec_y, dist_vec_z);
        INS_INC1(sqrt, 8);
        __m256 dist_vec_len = _mm256_sqrt_ps(dist_vec_dot2);

        // float m_scaled = m * 0.57735027;
        __m256 constant = _mm256_set1_ps(0.57735027f);
        INS_INC1(mul, 8);
        __m256 m_scaled = _mm256_mul_ps(m, constant);

        // return (!cond1 && !cond2 && !cond3) ? m_scaled : dist_vec_len;
        __m256 cond1_or_cond2 = _mm256_or_ps(cond1, cond2);
        __m256 cond1_or_cond2_or_cond3 = _mm256_or_ps(cond1_or_cond2, cond3);
        __m256 dist = _mm256_blendv_ps(m_scaled, dist_vec_len, cond1_or_cond2_or_cond3);

        _mm256_store_ps(res, dist);
    }

    static inline float octahedron_distance_short(const octa o, const vec from, const float current_min) {
        INS_INC(octa);
        vec pos = m33_mul_vec(o.inv_rot, vec_sub(from, o.center));

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = o.s + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

        INS_INC(octa_r);

        pos = vec_abs(pos);

        float s = o.s;

        INS_INC1(add, 3);
        float m = pos.x + pos.y + pos.z - s;
        vec q;

        if (3 * pos.x < m) {
            INS_INC1(mul, 1);
            INS_INC1(cmp, 1);
            q = pos;
        } else if (3 * pos.y < m) {
            INS_INC1(mul, 2);
            INS_INC1(cmp, 2);
            q = {pos.y, pos.x, pos.z};
        } else if (3 * pos.z < m) {
            INS_INC1(mul, 3);
            INS_INC1(cmp, 3);
            q = {pos.z, pos.x, pos.y};
        } else {
            INS_INC1(mul, 4);
            INS_INC1(cmp, 3);
            return m * 0.57735027;
        }

        INS_MUL;
        INS_INC1(add, 2);
        float k = clamp(0.5f * (q.z - q.y + s), 0.0f, s);

        INS_INC1(add, 3);
        float squared_distance = vec_dot2({q.x, q.y - s + k, q.z - k});

        INS_INC1(mul, 1);
        INS_CMP;
        if (squared_distance >= current_min * current_min) {
            return current_min;
        }
        return FSQRT(squared_distance);
    }

    /**
     * Computes the octahedron distance function with early termination.
     * Returns zero if early termination is possible, and a nonzero value otherwise.
     */
    static inline int octahedron_distance_short_vectorized(int idx, float* res, float* center_x, float* center_y, float* center_z, float* s_, float* inv_rot[3][3], const vec from, float current_min) {
        INS_INC1(octa, 8);

        __m256 c_x = _mm256_loadu_ps(center_x + idx);
        __m256 c_y = _mm256_loadu_ps(center_y + idx);
        __m256 c_z = _mm256_loadu_ps(center_z + idx);
        __m256 s = _mm256_loadu_ps(s_ + idx);

        __m256 inv_rot_00 = _mm256_loadu_ps(inv_rot[0][0] + idx);
        __m256 inv_rot_01 = _mm256_loadu_ps(inv_rot[0][1] + idx);
        __m256 inv_rot_02 = _mm256_loadu_ps(inv_rot[0][2] + idx);
        __m256 inv_rot_10 = _mm256_loadu_ps(inv_rot[1][0] + idx);
        __m256 inv_rot_11 = _mm256_loadu_ps(inv_rot[1][1] + idx);
        __m256 inv_rot_12 = _mm256_loadu_ps(inv_rot[1][2] + idx);
        __m256 inv_rot_20 = _mm256_loadu_ps(inv_rot[2][0] + idx);
        __m256 inv_rot_21 = _mm256_loadu_ps(inv_rot[2][1] + idx);
        __m256 inv_rot_22 = _mm256_loadu_ps(inv_rot[2][2] + idx);

        __m256 from_x = _mm256_set1_ps(from.x);
        __m256 from_y = _mm256_set1_ps(from.y);
        __m256 from_z = _mm256_set1_ps(from.z);
        __m256 curr_min = _mm256_set1_ps(current_min);

        // vec t = vec_sub(from, o.center);
        INS_INC1(add, 24);
        __m256 t_x = _mm256_sub_ps(from_x, c_x);
        __m256 t_y = _mm256_sub_ps(from_y, c_y);
        __m256 t_z = _mm256_sub_ps(from_z, c_z);

        // vec pos = m33_mul_vec(o.inv_rot, t);
        __m256 pos_x = vectorized_vec_dot(inv_rot_00, inv_rot_01, inv_rot_02, t_x, t_y, t_z);
        __m256 pos_y = vectorized_vec_dot(inv_rot_10, inv_rot_11, inv_rot_12, t_x, t_y, t_z);
        __m256 pos_z = vectorized_vec_dot(inv_rot_20, inv_rot_21, inv_rot_22, t_x, t_y, t_z);

        // enclosing sphere check
        __m256 pos_square = vectorized_vec_dot(pos_x, pos_y, pos_z, pos_x, pos_y, pos_z);
        INS_INC1(add, 8);
        __m256 duf_bound = _mm256_add_ps(s, curr_min);
        INS_INC1(mul, 8);
        __m256 duf_bound_square = _mm256_mul_ps(duf_bound, duf_bound);
        INS_INC1(cmp, 8);
        __m256 duf_mask = _mm256_cmp_ps(pos_square, duf_bound_square, _CMP_LT_OQ);
        int duf_mask_int = _mm256_movemask_ps(duf_mask);
        if (!duf_mask_int) {
            return 0;
        }

        INS_INC1(octa_r, 8);

        // vec pos_abs = vec_abs(pos);
        INS_INC1(abs, 24);
        __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        __m256 pos_abs_x = _mm256_and_ps(abs_mask, pos_x);
        __m256 pos_abs_y = _mm256_and_ps(abs_mask, pos_y);
        __m256 pos_abs_z = _mm256_and_ps(abs_mask, pos_z);

        // float s = o.s;

        // float m = pos_abs.x + pos_abs.y + pos_abs.z - s;
        INS_INC1(add, 24);
        __m256 pos_abs_xy = _mm256_add_ps(pos_abs_x, pos_abs_y);
        __m256 pos_abs_xyz = _mm256_add_ps(pos_abs_xy, pos_abs_z);
        __m256 m = _mm256_sub_ps(pos_abs_xyz, s);

        // vec q;

        // float pos_abs_x_times3 = 3 * pos_abs.x;
        // float pos_abs_y_times3 = 3 * pos_abs.y;
        // float pos_abs_z_times3 = 3 * pos_abs.z;
        __m256 three = _mm256_set1_ps(3.f);
        INS_INC1(mul, 24);
        __m256 pos_abs_x_times3 = _mm256_mul_ps(three, pos_abs_x);
        __m256 pos_abs_y_times3 = _mm256_mul_ps(three, pos_abs_y);
        __m256 pos_abs_z_times3 = _mm256_mul_ps(three, pos_abs_z);

        // int cond1 = pos_abs_x_times3 < m;
        // int cond2 = pos_abs_y_times3 < m;
        // int cond3 = pos_abs_z_times3 < m;
        INS_INC1(cmp, 24);
        __m256 cond1 = _mm256_cmp_ps(pos_abs_x_times3, m, _CMP_LT_OQ);
        __m256 cond2 = _mm256_cmp_ps(pos_abs_y_times3, m, _CMP_LT_OQ);
        __m256 cond3 = _mm256_cmp_ps(pos_abs_z_times3, m, _CMP_LT_OQ);

        // if (cond1) {
        //     q = pos_abs;
        // } else if (cond2) {
        //     q = {pos_abs.y, pos_abs.x, pos_abs.z};
        // } else {
        //     q = {pos_abs.z, pos_abs.x, pos_abs.y};
        // }

        __m256 q_branch_23_x = _mm256_blendv_ps(pos_abs_z, pos_abs_y, cond2);
        __m256 q_branch_23_y = pos_abs_x;
        __m256 q_branch_23_z = _mm256_blendv_ps(pos_abs_y, pos_abs_z, cond2);

        __m256 q_x = _mm256_blendv_ps(q_branch_23_x, pos_abs_x, cond1);
        __m256 q_y = _mm256_blendv_ps(q_branch_23_y, pos_abs_y, cond1);
        __m256 q_z = _mm256_blendv_ps(q_branch_23_z, pos_abs_z, cond1);

        // float qz_minus_qy = q.z - q.y;
        INS_INC1(add, 8);
        __m256 qz_minus_qy = _mm256_sub_ps(q_z, q_y);

        // float qz_minus_qy_plus_s = qz_minus_qy + s;
        INS_INC1(add, 8);
        __m256 qz_minus_qy_plus_s = _mm256_add_ps(qz_minus_qy, s);

        // float to_clamp = 0.5f * qz_minus_qy_plus_s;
        __m256 pointfive = _mm256_set1_ps(0.5f);
        INS_INC1(mul, 8);
        __m256 to_clamp = _mm256_mul_ps(pointfive, qz_minus_qy_plus_s);

        // float k = clamp(to_clamp, 0.0f, s);
        INS_INC1(cmp, 8);
        __m256 clamped_upper = _mm256_min_ps(to_clamp, s);
        __m256 zero = _mm256_setzero_ps();
        INS_INC1(cmp, 8);
        __m256 k = _mm256_max_ps(zero, clamped_upper);

        // vec dist_vec = {q.x, q.y - s + k, q.z - k};
        __m256 dist_vec_x = q_x;
        INS_INC1(add, 24);
        __m256 qy_minus_s = _mm256_sub_ps(q_y, s);
        __m256 dist_vec_y = _mm256_add_ps(qy_minus_s, k);
        __m256 dist_vec_z = _mm256_sub_ps(q_z, k);

        // float dist_vec_len = vec_length(dist_vec);
        __m256 dist_vec_dot2 = vectorized_vec_dot(dist_vec_x, dist_vec_y, dist_vec_z, dist_vec_x, dist_vec_y, dist_vec_z);

        // float m_scaled = m * 0.57735027;
        __m256 constant = _mm256_set1_ps(0.57735027f);
        INS_INC1(mul, 8);
        __m256 m_scaled = _mm256_mul_ps(m, constant);

        // for early termination, we want all intermediate results to be the square of the final result
        INS_INC1(mul, 8);
        __m256 m_scaled_square = _mm256_mul_ps(m_scaled, m_scaled);

        // return (!cond1 && !cond2 && !cond3) ? m_scaled : dist_vec_len;
        __m256 cond1_or_cond2 = _mm256_or_ps(cond1, cond2);
        __m256 cond1_or_cond2_or_cond3 = _mm256_or_ps(cond1_or_cond2, cond3);
        __m256 dist_square = _mm256_blendv_ps(m_scaled_square, dist_vec_dot2, cond1_or_cond2_or_cond3);

        // short circuit termination mask
        INS_INC1(mul, 8);
        __m256 upper_bound_square = _mm256_mul_ps(curr_min, curr_min);
        INS_INC1(cmp, 8);
        __m256 mask = _mm256_cmp_ps(dist_square, upper_bound_square, _CMP_LT_OQ);
        int mask_int = _mm256_movemask_ps(mask);

        if (!mask_int) {
            return 0;
        }

        INS_INC1(sqrt, 8);
        __m256 dist = _mm256_sqrt_ps(dist_square);

        _mm256_storeu_ps(res, dist);

        return 1;
    }

    static inline vec sphere_normal(const sphere sp, const vec pos) {
        INS_INC(sphere_n);
        return vec_normalize(vec_sub(pos, sp.center));
    }

    static inline vec box_normal(const box b, const vec from) {
        INS_INC(box_n);
        vec pos = m33_mul_vec(b.inv_rot, vec_sub(from, b.bottom_left));

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec sign = vec_div(pos, abs);
        vec q = vec_sub(abs, b.extents);

        // argmax(q.x, q.y, q.z)
        vec n_obj = {0, 0, 0};
        if (q.x > q.y && q.x > q.z && q.x > 0) {
            INS_INC1(cmp, 3);
            n_obj = {1, 0, 0};
        } else if (q.y > q.z && q.y > 0) {
            INS_INC1(cmp, 5);
            n_obj = {0, 1, 0};
        } else if (q.z > 0) {
            INS_INC1(cmp, 6);
            n_obj = {0, 0, 1};
        } else {
            INS_INC1(cmp, 6);
        }

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        return m33_mul_vec(b.rot, n_obj);
    }

    static inline vec plane_normal(const plane p, const vec) {
        INS_INC(plane_n);
        return p.normal;
    }

    static inline vec torus_normal(const torus t, const vec from) {
        INS_INC(torus_n);
        vec pos = m33_mul_vec(t.inv_rot, vec_sub(from, t.center));

        vec2 posxz = {pos.x, pos.z};

        INS_ADD;
        INS_DIV;
        posxz = vec2_scale(posxz, 1 - t.r1 / vec2_length(posxz));
        vec q = {posxz.x, pos.y, posxz.y};
        vec n_obj = vec_normalize(q);

        return m33_mul_vec(t.rot, n_obj);
    }

    static inline vec cone_normal(const cone c, const vec from) {
        INS_INC(cone_n);
        vec pos = m33_mul_vec(c.inv_rot, vec_sub(from, c.center));

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        // transform into rotation invariant subspace around y-axis
        vec2 posxz = {pos.x, pos.z};
        vec2 q = {vec2_length(posxz), pos.y};
        vec2 k1 = {r2, h};
        INS_ADD;
        INS_MUL;
        vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 2);
        INS_ABS;
        INS_CMP;
        vec2 ca = {q.x - min(q.x, (q.y < 0 ? r1 : r2)), fabsf(q.y) - h};
        INS_DIV;
        vec2 cb =
            vec2_add(vec2_sub(q, k1), vec2_scale(k2, clamp(vec2_dot(vec2_sub(k1, q), k2) / vec2_dot2(k2), 0.0f, 1.0f)));

        // invert transform from rotation invariant subspace
        vec n_obj = {0, 0, 0};
        INS_CMP;
        if (vec2_dot2(ca) < vec2_dot2(cb)) {
            INS_DIV;
            posxz = vec2_scale(posxz, ca.x / vec2_length(posxz));
            n_obj = {posxz.x, ca.y, posxz.y};
        } else {
            INS_DIV;
            posxz = vec2_scale(posxz, cb.x / vec2_length(posxz));
            n_obj = {posxz.x, cb.y, posxz.y};
        }
        n_obj = vec_normalize(n_obj);

        return m33_mul_vec(c.rot, n_obj);
    }

    static inline vec octahedron_normal(const octa o, const vec from) {
        INS_INC(octa_n);
        vec pos = m33_mul_vec(o.inv_rot, vec_sub(from, o.center));

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec sign = vec_div(pos, abs);

        vec n_obj = {1, 1, 1};
        n_obj = vec_normalize(n_obj);

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        return m33_mul_vec(o.rot, n_obj);
    }

} // namespace impl::vec2
