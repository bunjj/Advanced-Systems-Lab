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
    typedef vec (*normal_fun)(const struct shape s, const vec pos);

    struct shape {
        distance_fun distance;
        normal_fun normal;
        char data[std::max({sizeof(sphere), sizeof(plane), sizeof(box), sizeof(torus), sizeof(cone), sizeof(octa)})];
        // The matrix for transforming any point in the object space into the world
        // space.
        m44 matrix;
        // Inverse of the above matrix. Transforms points in the world space into
        // the object space.
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
        enum shape_type type;
    };

    struct scene {
        shape* shapes;
        int num_shapes;
        light* lights;
        int num_lights;
        camera cam;
    };

    extern struct scene scene;

    void load_scene(std::string& input);

} // namespace impl::opt1
