#pragma once

#include <algorithm>
#include <nlohmann/json.hpp>

#include "impl_opt1/geometry.h"

using json = nlohmann::json;

namespace impl::opt1 {

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
        m44 matrix;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct torus {
        vec center;
        float r1;
        float r2;
        m44 matrix;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct cone {
        vec center;
        float r1;
        float r2;
        float height;
        m44 matrix;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct octa {
        vec center;
        float s;
        m44 matrix;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
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
    static inline float sphere_distance(const sphere sp, const vec from) {
        INS_INC(sphere);
        INS_ADD;
        return vec_length(vec_sub(sp.center, from)) - sp.radius;
    }

    // Box
    static inline float box_distance(const box b, const vec from) {
        INS_INC(box);
        vec pos = vec4_to_vec(m44_mul_vec(b.inv_matrix, vec4_from_point(from)));
        vec q = vec_sub(vec_abs(pos), b.extents);
        return FADD(vec_length(vec_max(q, 0)), min(0.0f, max(max(q.x, q.y), q.z)));
    }

    // Plane
    static inline float plane_distance(const plane p, const vec from) {
        INS_INC(plane);
        return vec_dot(p.normal, vec_sub(from, p.point));
    }

    // Torus
    static inline float torus_distance(const torus t, const vec from) {
        INS_INC(torus);
        vec pos = vec4_to_vec(m44_mul_vec(t.inv_matrix, vec4_from_point(from)));
        vec2 posxz = {pos.x, pos.z};
        INS_ADD;
        vec2 q = {vec2_length(posxz) - t.r1, pos.y};

        INS_ADD;
        return vec2_length(q) - t.r2;
    }

    // Cone
    static inline float cone_distance(const cone c, const vec from) {
        INS_INC(cone);
        vec pos = vec4_to_vec(m44_mul_vec(c.inv_matrix, vec4_from_point(from)));

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

    // Octahedron
    static inline float octahedron_distance(const octa o, const vec from) {
        INS_INC(octa);
        vec pos = vec4_to_vec(m44_mul_vec(o.inv_matrix, vec4_from_point(from)));
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

    static inline vec sphere_normal(const sphere sp, const vec pos) {
        INS_INC(sphere_n);
        return vec_normalize(vec_sub(pos, sp.center));
    }

    static inline vec box_normal(const box b, const vec from) {
        INS_INC(box_n);
        vec pos = vec4_to_vec(m44_mul_vec(b.inv_matrix, vec4_from_point(from)));

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

        vec n_world = m44_rotate_only(b.matrix, n_obj);
        return n_world;
    }

    static inline vec plane_normal(const plane p, const vec) {
        INS_INC(plane_n);
        return p.normal;
    }

    static inline vec torus_normal(const torus t, const vec from) {
        INS_INC(torus_n);
        vec pos = vec4_to_vec(m44_mul_vec(t.inv_matrix, vec4_from_point(from)));

        vec2 posxz = {pos.x, pos.z};

        INS_ADD;
        INS_DIV;
        posxz = vec2_scale(posxz, 1 - t.r1 / vec2_length(posxz));
        vec q = {posxz.x, pos.y, posxz.y};
        vec n_obj = vec_normalize(q);

        vec n_world = m44_rotate_only(t.matrix, n_obj);
        return n_world;
    }

    static inline vec cone_normal(const cone c, const vec from) {
        INS_INC(cone_n);
        vec pos = vec4_to_vec(m44_mul_vec(c.inv_matrix, vec4_from_point(from)));

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

        vec n_world = m44_rotate_only(c.matrix, n_obj);
        return n_world;
    }

    static inline vec octahedron_normal(const octa o, const vec from) {
        INS_INC(octa_n);
        vec pos = vec4_to_vec(m44_mul_vec(o.inv_matrix, vec4_from_point(from)));

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec sign = vec_div(pos, abs);

        vec n_obj = {1, 1, 1};
        n_obj = vec_normalize(n_obj);

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        vec n_world = m44_rotate_only(o.matrix, n_obj);
        return n_world;
    }

} // namespace impl::opt1
