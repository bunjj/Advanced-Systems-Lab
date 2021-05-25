#pragma once

#include <algorithm>
#include <nlohmann/json.hpp>

#include "impl_opt5/geometry.h"

using json = nlohmann::json;

namespace impl::opt5 {

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
        float k2d2inv; // 1 / vec_dot2(k2)
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
    };

    extern struct scene scene;

    void load_scene(std::string& input);
    void from_ref_scene();

    // distance and normal functions

    // Sphere
    static inline float sphere_distance_short(const sphere sp, const vec p_world, const float current_min) {
        INS_INC(sphere);
        INS_ADD;
        float upper_bound = current_min + sp.radius;
        float squared_distance = vec_dot2(vec_sub(sp.center, p_world));
        INS_MUL;
        INS_CMP;
        if (squared_distance >= upper_bound * upper_bound) {
            return current_min;
        }
        INS_ADD;
        return FSQRT(squared_distance) - sp.radius;
    }

    // Box
    static inline float box_distance_short(const box b, const vec pos, const float current_min) {
        INS_INC(box);

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = b.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

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
        INS_ADD;
        return FSQRT(intermediate_squared_dist) + extent_values;
    }

    // Plane
    static inline float plane_distance(const plane p, const vec p_world) {
        INS_INC(plane);
        return vec_dot(p.normal, vec_sub(p_world, p.point));
    }

    // Torus
    static inline float torus_distance_short(const torus t, const vec pos, const float current_min) {
        INS_INC(torus);

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = t.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;
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

    // Cone
    static inline float cone_distance_short(const cone c, const vec pos, const float current_min) {
        INS_INC(cone);

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = c.r + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        vec2 q = {vec2_length({pos.x, pos.z}), pos.y};
        vec2 k1 = {r2, h};
        // Computing k2 again instead of storing it in the shape seems to perform better. Maybe because the compiler can
        // optimize better?
        INS_ADD;
        INS_MUL;
        vec2 k2 = {r2 - r1, 2 * h};
        INS_INC1(add, 2);
        INS_ABS;
        INS_CMP;
        vec2 ca = {q.x - min(q.x, (q.y < 0 ? r1 : r2)), fabsf(q.y) - h};
        INS_MUL;
        vec2 cb = vec2_add(vec2_sub(q, k1), vec2_scale(k2, clamp(vec2_dot(vec2_sub(k1, q), k2) * c.k2d2inv, 0.0f, 1.0f)));

        float squared_min = min(vec2_dot2(ca), vec2_dot2(cb));
        INS_MUL;
        INS_CMP;
        if (squared_min >= current_min * current_min) {
            return current_min;
        }

        INS_INC1(cmp, 2);
        float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;

        INS_MUL;
        INS_SQRT;
        return s * sqrtf(squared_min);
    }

    // Octahedron
    static inline float octahedron_distance_short(const octa o, const vec pos, const float current_min) {
        INS_INC(octa);

        //computation of DUF
        INS_ADD;
        INS_MUL;
        INS_CMP;
        float pos_squared = vec_dot2(pos);
        float duf_bound = o.s + current_min;
        if( pos_squared >= duf_bound * duf_bound) return current_min;

        vec abs = vec_abs(pos);

        INS_INC1(add, 3);
        float m = abs.x + abs.y + abs.z - o.s;
        float qx;
        float qy;
        float qz;

        if (3 * abs.x < m) {
            INS_INC1(mul, 1);
            INS_INC1(cmp, 1);
            qx = abs.x;
            qy = abs.y;
            qz = abs.z;
        } else if (3 * abs.y < m) {
            INS_INC1(mul, 2);
            INS_INC1(cmp, 2);
            qx = abs.y;
            qy = abs.x;
            qz = abs.z;
        } else if (3 * abs.z < m) {
            INS_INC1(mul, 3);
            INS_INC1(cmp, 3);
            qx = abs.z;
            qy = abs.x;
            qz = abs.y;
        } else {
            INS_INC1(mul, 4);
            INS_INC1(cmp, 3);
            return m * SQRT3_INV;
        }

        INS_MUL;
        INS_INC1(add, 2);
        float k = clamp(0.5f * (qz - qy + o.s), 0.0f, o.s);

        INS_INC1(add, 3);
        float squared_distance = vec_dot2({qx, qy - o.s + k, qz - k});

        INS_INC1(mul, 1);
        INS_CMP;
        if (squared_distance >= current_min * current_min) {
            return current_min;
        }
        return FSQRT(squared_distance);
    }

    static inline vec sphere_normal(const sphere sp, const vec p_world) {
        INS_INC(sphere_n);
        return vec_normalize(vec_sub(p_world, sp.center));
    }

    static inline vec box_normal(const box b, const vec pos) {
        INS_INC(box_n);

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec q = vec_sub(abs, b.extents);

        vec n_obj = {0.f , 0.f, 0.f};

        // We only need to find the normal of the side we touched
        if (q.x > q.y && q.x > q.z && q.x > 0) {
            INS_INC1(cmp, 4);
            n_obj = {pos.x >= 0.f ? 1.f : -1.f, 0.f, 0.f};
        } else if (q.y > q.z && q.y > 0) {
            INS_INC1(cmp, 6);
            n_obj = {0.f, pos.y >= 0.f ? 1.f : -1.f, 0.f};
        } else if (q.z > 0) {
            INS_INC1(cmp, 7);
            n_obj = {0.f, 0.f, pos.z >= 0.f ? 1.f : -1.f};
        } else {
            INS_INC1(cmp, 6);
        }

        return m33_mul_vec(b.rot, n_obj);
    }

    static inline vec plane_normal(const plane p, const vec) {
        INS_INC(plane_n);
        return p.normal;
    }

    static inline vec torus_normal(const torus t, const vec pos) {
        INS_INC(torus_n);

        vec2 posxz = {pos.x, pos.z};

        INS_ADD;
        INS_DIV;
        posxz = vec2_scale(posxz, 1 - t.r1 / vec2_length(posxz));
        vec q = {posxz.x, pos.y, posxz.y};
        vec n_obj = vec_normalize(q);

        return m33_mul_vec(t.rot, n_obj);
    }

    static inline vec cone_normal(const cone c, const vec pos) {
        INS_INC(cone_n);

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
        INS_MUL;
        vec2 cb = vec2_add(vec2_sub(q, k1), vec2_scale(k2, clamp(vec2_dot(vec2_sub(k1, q), k2) * c.k2d2inv, 0.0f, 1.0f)));

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

    static inline vec octahedron_normal(const octa o, const vec pos) {
        INS_INC(octa_n);

        INS_INC1(cmp, 3);
        float x = pos.x >= 0.f ? SQRT3_INV : -SQRT3_INV;
        float y = pos.y >= 0.f ? SQRT3_INV : -SQRT3_INV;
        float z = pos.z >= 0.f ? SQRT3_INV : -SQRT3_INV;
        return m33_mul_vec(o.rot, {x, y, z});
    }

} // namespace impl::opt5
