#include "dbg.h"
#include <inttypes.h>
#include <stdio.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>
#include "nlohmann/json.hpp"

using json = nlohmann::json;
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

float vec_norm(vec v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

float vec_length(vec v) { return sqrtf(vec_norm(v)); }

#define VEC_OP(v1, v2, OP) \
    vec { (v1).x OP(v2).x, (v1).y OP(v2).y, (v1).z OP(v2).z }

vec vec_add(vec v1, vec v2) { return VEC_OP(v1, v2, +); }

vec vec_sub(vec v1, vec v2) { return VEC_OP(v1, v2, -); }

vec vec_scale(vec v1, float factor) {
    return VEC_OP(v1, (vec{factor, factor, factor}), *);
}

vec vec_normalize(vec v) {
    float len = vec_length(v);
    return vec{v.x / len, v.y / len, v.z / len};
}

float vec_dot(vec v1, vec v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec vec_abs(vec v) { return {fabsf(v.x), fabsf(v.y), fabsf(v.z)}; }

vec vec_max(vec v, float val) {
    return {std::max(v.x, val), std::max(v.y, val), std::max(v.z, val)};
}

// }}}

// Vec2 {{{
struct vec2 {
    float x;
    float y;
};

float vec2_length(vec2 v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

#define VEC2_OP(v1, v2, OP) \
    vec2{ (v1).x OP(v2).x, (v1).y OP(v2).y }

vec2 vec2_sub(vec2 v1, vec2 v2) {
    return VEC2_OP(v1, v2, -);
}

vec2 vec2_add(vec2 v1, vec2 v2) {
    return VEC2_OP(v1, v2, +);
}

vec2 vec2_scale(vec2 v, float f) {
    return VEC2_OP(v, (vec2{f, f}), *);
}

float vec2_dot(vec2 v1, vec2 v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

float vec2_dot2(vec2 v) {
    return vec2_dot(v, v);
}
// }}}

// Vec4 {{{

struct vec4 {
    float x;
    float y;
    float z;
    float t;
};

vec4 vec4_init(float* values) {
    return {values[0], values[1], values[2], values[3]};
}

float vec4_dot(vec4 v1, vec4 v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.t * v2.t;
}

vec4 vec4_from_point(vec p) { return {p.x, p.y, p.z, 1}; }

vec4 vec4_from_dir(vec p) { return {p.x, p.y, p.z, 0}; }

vec vec4_to_vec(vec4 p) { return {p.x, p.y, p.z}; }

vec4 get_base_vec4(int idx) {
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

    float det = a11 * (a22 * A2323 - a23 * A1323 + a24 * A1223) -
                a12 * (a21 * A2323 - a23 * A0323 + a24 * A0223) +
                a13 * (a21 * A1323 - a22 * A0323 + a24 * A0123) -
                a14 * (a21 * A1223 - a22 * A0223 + a23 * A0123);
    det = 1 / det;

    m44 inv;

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

// Scene {{{
struct camera {
    float fov;
    vec pos;
    vec rotation;
};

struct light {
    vec pos;
    vec color;
    float intensity;
};

struct sphere {
    vec center;
    float radius;
};

struct plane {
    // Normal vector
    vec normal;
    // Point on the plane
    vec point;
};

struct box {
    vec bottom_left;
    vec extents;
};

struct torus {
    vec center;
    float r1;
    float r2;
};

struct cone {
    vec center;
    float r1;
    float r2;
    float height;
};

struct octa {
    vec center;
    float s;
};

typedef float (*distance_fun)(const struct shape s, const vec pos);

struct shape {
    distance_fun distance;
    char data[std::max({sizeof(sphere), sizeof(plane), sizeof(box), sizeof(torus), sizeof(cone), sizeof(octa)})];
    // The matrix for transforming any point in the object space into the world
    // space.
    m44 matrix;
    // Inverse of the above matrix. Transforms points in the world space into
    // the object space.
    m44 inv_matrix;
};

shape make_shape(const m44 matrix, distance_fun f, void* data,
                 size_t data_size) {
    shape shap;
    shap.distance = f;
    shap.matrix = matrix;
    shap.inv_matrix = m44_inv(matrix);
    memcpy(&shap.data, data, data_size);
    return shap;
}

// Sphere {{{
float sphere_distance(const shape s, const vec from) {
    sphere sp = *((sphere*)s.data);
    return vec_length(vec_sub(sp.center, from)) - sp.radius;
}

vec sphere_normal(sphere s, vec pos) {
    return vec_normalize(vec_sub(pos, s.center));
}

shape make_sphere(float x, float y, float z, float r) {
    sphere s = {.center = {x, y, z}, .radius = r};
    return make_shape(get_transf_matrix({x, y, z}, {0, 0, 0}), sphere_distance,
                      &s, sizeof(s));
}

shape load_sphere(json& j){
    float x,y,z,r;
    x = j["position"]["x"];
    y = j["position"]["y"];
    z = j["position"]["z"];
    r = j["params"]["radius"];
    return make_sphere(x,y,z,r);

}

// }}}

// Box {{{
float box_distance(const shape s, const vec from) {
    box b = *((box*)s.data);
    vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
    vec q = vec_sub(vec_abs(pos), b.extents);
    return vec_length(vec_max(q, 0)) +
           std::min(0.0f, std::max({q.x, q.y, q.z}));
}

shape make_box(vec bottom_left, vec extents, vec rot) {
    box s = {bottom_left, extents};
    return make_shape(get_transf_matrix(bottom_left, rot), box_distance, &s,
                      sizeof(s));
}

shape load_box(json& j){
    vec pos, extents, rot;
    pos.x = j["position"]["x"];
    pos.y = j["position"]["y"];
    pos.z = j["position"]["z"];
    extents.x = j["params"]["extents"]["x"];
    extents.y = j["params"]["extents"]["y"];
    extents.z = j["params"]["extents"]["z"];
    rot.x = j["rotation"]["x"];
    rot.y = j["rotation"]["y"];
    rot.z = j["rotation"]["z"];
    return make_box(pos, extents, rot);
}

// }}}

// Plane {{{

float plane_distance(const shape s, const vec from) {
    plane p = *((plane*)s.data);
    return vec_dot(p.normal, vec_sub(from, p.point));
}

shape make_plane(vec normal, vec point) {
    plane p = {.normal = vec_normalize(normal), .point = point};
    return make_shape(identity, plane_distance, &p, sizeof(p));
}

shape load_plane(json& j){
    vec normal, point;
    float displacement = j["params"]["displacement"];
    normal.x = j["params"]["normal"]["x"];
    normal.y = j["params"]["normal"]["y"];
    normal.z = j["params"]["normal"]["z"];
    point = vec_scale(normal, displacement);
    return make_plane(normal, point);  
    
}
// }}}

// Torus {{{
float torus_distance(const shape s, const vec from) {
    torus t = *((torus*)s.data);
    vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
    vec2 posxz = {pos.x, pos.z};
    vec2 q = {vec2_length(posxz) - t.r1, pos.y};

    return vec2_length(q) - t.r2;
}

shape make_torus(vec center, float r1, float r2, vec rot) {
    torus t = {center, r1, r2};
    return make_shape(get_transf_matrix(center, rot), torus_distance, &t,
                      sizeof(t));
}

shape load_torus(json& j){
    vec pos, rot;
    float r1, r2;
    pos.x = j["position"]["x"];
    pos.y = j["position"]["y"];
    pos.z = j["position"]["z"];
    rot.x = j["rotation"]["x"];
    rot.y = j["rotation"]["y"];
    rot.z = j["rotation"]["z"];

    r1 = j["params"]["r1"];
    r2 = j["params"]["r2"];

    return make_torus(pos, r1, r2, rot);
}

// }}}

// Cone {{{
float cone_distance(const shape shap, const vec from) {
    cone c = *((cone*)shap.data);
    vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));

    float r1 = c.r1;
    float r2 = c.r2;
    float h = c.height;

    vec2 q = {vec2_length({pos.x, pos.z}), pos.y};
    vec2 k1 = {r2, h};
    vec2 k2 = {r2 - r1, 2 * h};
    vec2 ca = {q.x - std::min(q.x, (q.y < 0? r1 : r2)), fabsf(q.y) - h};
    vec2 cb = vec2_add(vec2_sub(q, k1), vec2_scale(k2, std::clamp(vec2_dot(vec2_sub(k1, q), k2) / vec2_dot2(k2), 0.0f, 1.0f)));
    float s = (cb.x < 0 && ca.y < 0) ? -1 : 1;

    return s * std::sqrt(std::min(vec2_dot2(ca), vec2_dot2(cb)));
}

shape make_cone(vec center, float r1, float r2, float height, vec rot) {
    cone c = {center, r1, r2, height};
    return make_shape(get_transf_matrix(center, rot), cone_distance, &c, sizeof(c));
}

shape load_cone(json& j){
    //TODO: adjust cone definition according to mail!
    vec pos, rot;
    float r1,r2,height;
    pos.x = j["position"]["x"];
    pos.y = j["position"]["y"];
    pos.z = j["position"]["z"];
    rot.x = j["rotation"]["x"];
    rot.y = j["rotation"]["y"];
    rot.z = j["rotation"]["z"];

    r1 = j["params"][0];
    r2 = j["params"][1];
    height = j["params"][2];

    return make_cone(pos, r1,r2,height, rot);
}

// }}}

// Octahedron {{{
float octahedron_distance(const shape shap, const vec from) {
    octa o = *((octa*)shap.data);
    vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));
    pos = vec_abs(pos);

    float s = o.s;

    float m = pos.x + pos.y + pos.z - s;
    vec q;

    if (3 * pos.x < m) {
        q = pos;
    } else if (3 * pos.y < m) {
        q = {pos.y, pos.x, pos.z};
    } else if (3 * pos.z < m) {
        q = {pos.z, pos.x, pos.y};
    } else {
        return m * 0.57735027;
    }

    float k = std::clamp(0.5f * (q.z - q.y + s), 0.0f, s);

    return vec_length({q.x, q.y - s + k, q.z - k});
}

shape make_octahedron(vec center, float s, vec rot) {
    octa o = {center, s};
    return make_shape(get_transf_matrix(center, rot), octahedron_distance, &o, sizeof(o));
}

shape load_octa(json& j){
    vec pos, rot; float s;
    pos.x = j["position"]["x"];
    pos.y = j["position"]["y"];
    pos.z = j["position"]["z"];
    rot.x = j["rotation"]["x"];
    rot.y = j["rotation"]["y"];
    rot.z = j["rotation"]["z"];

    s = j["params"]["s"];
    return make_octahedron(pos, s, rot);

}

// }}}

#ifdef SCENE0

// clang-format off
static shape shapes[] = {
    make_plane({0, 1, 0}, {0, -3, 0}),
    make_box({-4, -1.5, 15}, {0.25, 0.5, 1}, {-1.5, -1.5, 12}),
    make_sphere(0, 0, 20, 3),         
    make_cone({1.5, -1.5, 12}, 1, 0.5, 1, {0, 0, 0}),
    make_torus({4, -1.5, 15}, 1, 0.5, {0, 0, 0}),
    make_octahedron({-1.5, -1.5, 12}, 1, {0, 0, 0}),
};
// clang-format on
static int num_shapes = sizeof(shapes) / sizeof(shape);

/*
 * Static light sources
 */
static light lights[] = {
    {{0, 100, 0}, {200.0f/255, 200.0f/255, 200.0f/255}, 150000},
};

static int num_lights = sizeof(lights) / sizeof(light);

static vec camera_pos = {0, 0, 0};
static vec camera_rot = {0, 0, 0};

// Field of view in degrees
static float fov = 30;

#else


// clang-format off
static shape shapes[] = {
    make_sphere(1, -5, 20, 7),         
    make_sphere(0, 1, 15, 2),
    make_sphere(0, 120, 200, 100),     
    make_plane({0, 0, -1}, {0, 0, 200}),
    make_plane({0, 1, 0}, {0, -20, 0}),
    make_box({40, 10, 80}, {10, 20, 10}, {0, -30, 0}),
};
// clang-format on
static int num_shapes = sizeof(shapes) / sizeof(shape);



/*
 * Static light sources
 */
static light lights[] = {
    {{0, 10, 0}, {1.0, 0.9, 0.7}, 3000},
    {{20, 10, 15}, {0, 0, 1.0}, 6000},
    {{-20, 10, 15}, {0.847, 0.2588, 0.2588}, 6000},
    {{-100, 40, 80}, {0, 0.9, 0.7}, 90000},
    {{100, 40, 80}, {1.0, 0, 0.7}, 90000},
};

static int num_lights = sizeof(lights) / sizeof(light);

static vec camera_pos = {0, 0, 0};
static vec camera_rot = {0, 0, 0};

// Field of view in degrees
static float fov = 80;
#endif
static shape* shape_load;
// }}}

// Sphere Tracing {{{

// max distance
static float D = 2048;
static float EPS = 0.001;

struct hit {
    bool is_hit;
    float distance;
    int steps;
    vec color;
};

static bool sphere_trace_shadow(vec point, vec light_dir, float max_distance) {
    // TODO if we start with t = 0 this function causes some pixels to be black
    // because it erroneously detects a collision with the original object
    float t = EPS;

    while (t < max_distance) {
        vec pos = vec_add(point, vec_scale(light_dir, t));

        float min_distance = INFINITY;

        for (int k = 0; k < num_shapes; k++) {
            float distance = shapes[k].distance(shapes[k], pos);

            if (distance < min_distance) {
                min_distance = distance;

                if (min_distance <= EPS * t) {
                    return true;
                }
            }
        }

        t += min_distance;
    }

    return false;
}

static hit sphere_trace(vec origin, vec dir) {
    float t = 0;

    int steps = 0;

    while (t < D) {
        vec pos = vec_add(origin, vec_scale(dir, t));

        float min_distance = INFINITY;
        int shape_idx = -1;

        for (int k = 0; k < num_shapes; k++) {
            float distance = shapes[k].distance(shapes[k], pos);

            if (distance < min_distance) {
                min_distance = distance;
                shape_idx = k;
            }
        }

        if (min_distance < EPS) {
            shape s = shapes[shape_idx];
            static const float delta = 10e-5;
            vec delta1 = {delta, 0, 0};
            vec delta2 = {0, delta, 0};
            vec delta3 = {0, 0, delta};
            // Some shapes can calculate this directly
            vec normal = vec_normalize({
                s.distance(s, vec_add(pos, delta1)) -
                    s.distance(s, vec_sub(pos, delta1)),
                s.distance(s, vec_add(pos, delta2)) -
                    s.distance(s, vec_sub(pos, delta2)),
                s.distance(s, vec_add(pos, delta3)) -
                    s.distance(s, vec_sub(pos, delta3)),
            });
            vec color{0, 0, 0};

            for (int i = 0; i < num_lights; i++) {
                vec light_point = vec_sub(lights[i].pos, pos);

                if (vec_dot(light_point, normal) > 0) {
                    vec light_dir = vec_normalize(light_point);
                    // Squared distance from light to point
                    float dist_sq = vec_norm(light_point);

                    if (!sphere_trace_shadow(pos, light_dir, sqrt(dist_sq))) {
                        color =
                            vec_add(color, vec_scale(lights[i].color,
                                                   vec_dot(light_dir, normal) *
                                                       lights[i].intensity /
                                                       (4 * M_PI_F * dist_sq)));
                    }
                }
            }

            return {true, t, steps, color};
        }

        t += min_distance;
        steps++;
    }

    return hit{false, t, steps, 0};
}

// }}}

static void dump_image(int width, int height, const float* pixels) {
    printf("P3\n%d %d\n255\n", width, height);

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                uint8_t channel =
                    std::min(1.0f, pixels[3 * (width * j + i) + k]) * 255;
                printf("%03" PRIu8 " ", channel);
            }
        }

        printf("\n");
    }
}


/**
 * 
 * Methods to read in the data from the json file
 * 
 */

static camera load_camera(json& j){
    camera cam;
    cam.fov = j["camera"]["fov"];
    cam.pos.x = j["camera"]["position"]["x"];
    cam.pos.y = j["camera"]["position"]["y"];
    cam.pos.z = j["camera"]["position"]["z"];
    cam.rotation.x = j["camera"]["rotation"]["x"];
    cam.rotation.y = j["camera"]["rotation"]["y"];
    cam.rotation.z = j["camera"]["rotation"]["z"];

    return cam;
}

static light load_light(json& j){
    //Loads a single light since all the scenes seem to only have one light.
    //But can simply be adjusted to load multiple lights instead
    light l;
    l.pos.x = j["pointlight"]["position"]["x"];
    l.pos.y = j["pointlight"]["position"]["y"];
    l.pos.z = j["pointlight"]["position"]["z"];
    l.color.x = j["pointlight"]["emission"]["x"];
    l.color.y = j["pointlight"]["emission"]["y"];
    l.color.z = j["pointlight"]["emission"]["z"];
    l.intensity = 150000;

    return l;
}


static void load_shapes(json& j){
    num_shapes = j["objects"].size();
    //TODO check out implementation of shapes!
    shape_load = (shape*) malloc(sizeof(shape) * num_shapes);

    for(int i = 0; i < num_shapes; i++){
        shape new_shape;
        json current_shape = j["objects"][i];
        std::string current = current_shape["kind"].get<std::string>();
        if(current == "sphere"){
           new_shape = load_sphere(current_shape); 
        }else if(current == "plane"){
            new_shape = load_plane(current_shape);
        }else if(current == "box"){
            new_shape = load_box(current_shape);
        }else if(current == "torus"){
            new_shape = load_torus(current_shape);
        }else if(current == "cone"){
            new_shape = load_cone(current_shape);
        }else if(current == "octahedron"){
            new_shape = load_octa(current_shape);
        }else{
            std::cout << "something went wrong" << std::endl;
            break;
        }
        shape_load[i] = new_shape;

    }

}

/**
 * Very simple sphere tracer
 *
 * All objects are in world coordinates.
 *
 * The camera can be positioned and rotated arbitrarily
 *
 * There are hard-coded spheres and lights.
 *
 */
int main(void) {
    // Height of the resulting image in pixels
    int height = 480;
    int width = 640;

    //Path to scene
    std::string path = "../../scenes/scene0.json";
    std::ifstream i(path);
    json j;
    i>> j;

    //load camera paremters
    camera cam = load_camera(j);
    //load light (as defined in scene, not with intensity parameter)
    light l = load_light(j);
    //don't have to save shapes since it uses static variable
    load_shapes(j);



    m44 camera_matrix = get_transf_matrix(camera_pos, camera_rot);

    float aspect_ratio = static_cast<float>(width) / height;

    float fov_factor = tanf(TO_RAD(fov / 2));

    auto pixels = std::make_unique<float[]>(height * width * 3);

    vec origin =
        vec4_to_vec(m44_mul_vec(camera_matrix, vec4_from_point({0, 0, 0})));

    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            /*
             * Position of the pixel in camera space.
             *
             * We assume that the camera is looking towards positive z and the
             * image plane is one unit away from the camera (z = -1 in this
             * case).
             */
            float x = (2 * (px + 0.5) / width - 1) * aspect_ratio * fov_factor;
            float y = (1 - 2 * (py + 0.5) / height) * fov_factor;
            float z = 1;

            // Direction in camera space.
            vec dir = vec_normalize({x, y, z});

            vec4 world_dir = m44_mul_vec(camera_matrix, vec4_from_dir(dir));

            auto h = sphere_trace(origin, vec4_to_vec(world_dir));

            vec color = h.is_hit ? h.color : vec{0, 0, 0};

            pixels[3 * (width * py + px)] = color.x;
            pixels[3 * (width * py + px) + 1] = color.y;
            pixels[3 * (width * py + px) + 2] = color.z;
        }
    }
    dump_image(width, height, pixels.get());

    return 0;
}
