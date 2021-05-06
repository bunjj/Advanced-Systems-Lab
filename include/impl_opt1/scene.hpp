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
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct box {
        vec bottom_left;
        vec extents;
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct torus {
        vec center;
        float r1;
        float r2;
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
        m44 inv_matrix;
        vec color;
        float reflection;
        float shininess;
    };

    struct octa {
        vec center;
        float s;
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

    // distance functions
    float sphere_distance(const sphere sp, const vec from);
    float box_distance(const box b, const vec from);
    float plane_distance(const plane p, const vec from);
    float torus_distance(const torus t, const vec from);
    float cone_distance(const cone c, const vec from);
    float octahedron_distance(const octa o, const vec from);

    // normal functions
    vec sphere_normal(sphere s, vec pos);
    vec box_normal(box s, vec pos);
    vec plane_normal(plane s, vec pos);
    vec torus_normal(torus s, vec pos);
    vec cone_normal(cone s, vec pos);
    vec octahedron_normal(octa s, vec pos);

} // namespace impl::opt1
