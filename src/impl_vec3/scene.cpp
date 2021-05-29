/**
 * Contains all functions related to the scene in the reference implementation.
 * All other implementations will use this to construct their scenes.
 */
#include "impl_vec3/scene.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "impl_ref/scene.hpp"

using json = nlohmann::json;

namespace impl::vec3 {

    static vec load_vec(json& j) {
        return {j["x"], j["y"], j["z"]};
    }

    static vec load_pos(json& j) {
        return load_vec(j["position"]);
    }

    static vec load_rot(json& j) {
        return load_vec(j["rotation"]);
    }

    /**
     * Global variable containing the scene for this implementation
     */
    struct scene scene;

    // Sphere {{{

    sphere make_sphere(vec center, float r, vec color, float reflection, float shininess) {
        sphere s;
        s.center = center;
        s.radius = r;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    sphere load_sphere(json& j) {
        float r;
        vec center = load_pos(j);
        r = j["params"]["radius"];
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];
        return make_sphere(center, r, color, reflection, shininess);
    }

    // }}}

    // Box {{{

    box make_box(vec bottom_left, vec extents, vec color, float reflection, float shininess, m33 rot, m33 rot_matrix) {
        box s;
        s.bottom_left = bottom_left;
        s.extents = extents;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        s.rot = rot;
        s.inv_rot = rot_matrix;
        s.r = vec_length(extents);
        return s;
    }

    box load_box(json& j) {
        vec pos = load_pos(j);
        vec extents = load_vec(j["params"]["extents"]);
        vec rot = load_rot(j);
        m33 rot_m = get_rot_matrix_33(rot);
        m33 inv_rot = m33_inv(rot_m);
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];
        return make_box(pos, extents, color, reflection, shininess, rot_m, inv_rot);
    }

    // }}}

    // Plane {{{

    plane make_plane(vec normal, vec point, vec color, float reflection, float shininess) {
        plane s;
        s.normal = vec_normalize(normal);
        s.point = point;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    plane load_plane(json& j) {
        float displacement = j["params"]["displacement"];
        vec normal = load_vec(j["params"]["normal"]);
        vec point = vec_scale(normal, displacement);
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];
        return make_plane(normal, point, color, reflection, shininess);
    }

    // }}}

    // Torus {{{

    torus make_torus(vec center, float r1, float r2, vec color, float reflection, float shininess,
        m33 rot, m33 inv_rot) {
        torus s;
        s.center = center;
        s.r1 = r1;
        s.r2 = r2;
        s.r = r1 + r2;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        s.rot = rot;
        s.inv_rot = inv_rot;
        return s;
    }

    torus load_torus(json& j) {
        float r1, r2;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        r1 = j["params"]["r1"];
        r2 = j["params"]["r2"];

        m33 rot_m = get_rot_matrix_33(rot);
        m33 inv_rot = m33_inv(rot_m);
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_torus(pos, r1, r2, color, reflection, shininess, rot_m, inv_rot);
    }

    // }}}

    // Cone {{{

    cone make_cone(vec center, float r1, float r2, float height, vec color, float reflection,
        float shininess, m33 rot, m33 inv_rot) {
        cone s;
        s.center = center;
        s.r1 = r1;
        s.r2 = r2;
        s.height = height;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        s.rot = rot;
        s.inv_rot = inv_rot;
        s.r = vec2_length({height, max(r1,r2)});
        vec2 k2 = {r2 - r1, 2 * height};
        s.k2d2inv = 1.f / vec2_dot2(k2);
        return s;
    }

    cone load_cone(json& j) {
        // TODO: adjust cone definition according to mail!
        float r1, r2, height;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        r1 = j["params"][0];
        r2 = j["params"][1];
        height = j["params"][2];

        m33 rot_m = get_rot_matrix_33(rot);
        m33 inv_rot = m33_inv(rot_m);

        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_cone(pos, r1, r2, height, color, reflection, shininess, rot_m, inv_rot);
    }

    // }}}

    // Octahedron {{{

    octa make_octahedron(
        vec center, float s_param, vec color, float reflection, float shininess, m33 rot, m33 inv_rot) {
        octa s;
        s.center = center;
        s.s = s_param;
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        s.rot = rot;
        s.inv_rot = inv_rot;
        return s;
    }

    octa load_octa(json& j) {
        float s;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        s = j["params"]["s"];

        m33 rot_m = get_rot_matrix_33(rot);
        m33 inv_rot = m33_inv(rot_m);

        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_octahedron(pos, s, color, reflection, shininess, rot_m, inv_rot);
    }

    // }}}

    // Load JSON {{{

    /**
     *
     * Methods to read in the data from the json file
     *
     */
    static camera load_camera(json& j) {
        camera cam;
        json c = j["camera"];
        cam.fov = j["camera"]["fov"];
        cam.pos = load_pos(c);
        cam.rotation = load_rot(c);
        return cam;
    }

    static light load_single_light(json& j) {
        light l;
        l.pos = load_pos(j);
        l.emission = load_vec(j["emission"]);
        return l;
    }

    static void load_light(json& j) {
        json light = j["pointlight"];

        scene.num_lights = light.is_array() ? light.size() : 1;

        scene.lights = (struct light*)malloc(sizeof(struct light) * scene.num_lights);

        if (light.is_array()) {
            int i = 0;
            for (auto& l : light) {
                scene.lights[i] = load_single_light(l);
                i++;
            }
        } else {
            scene.lights[0] = load_single_light(light);
        }
    }


    static void load_shapes(json& j) {
        int num_shapes = j["objects"].size();

        scene.num_spheres = 0;
        scene.num_planes = 0;
        scene.num_boxes = 0;
        scene.num_tori = 0;
        scene.num_cones = 0;
        scene.num_octahedra = 0;

        // first pass through scene to determine number of shapes of each kind
        for (int i = 0; i < num_shapes; i++) {
            json current_shape = j["objects"][i];
            std::string current = current_shape["kind"].get<std::string>();
            if (current == "sphere") {
                scene.num_spheres++;
            } else if (current == "plane") {
                scene.num_planes++;
            } else if (current == "box") {
                scene.num_boxes++;
            } else if (current == "torus") {
                scene.num_tori++;
            } else if (current == "cone") {
                scene.num_cones++;
            } else if (current == "octahedron") {
                scene.num_octahedra++;
            } else {
                throw std::runtime_error("Unknown shape " + current);
            }
        }

        // allocate memory for the shape arrays
        scene.spheres = (sphere*)malloc(sizeof(sphere) * scene.num_spheres);
        scene.planes = (plane*)malloc(sizeof(plane) * scene.num_planes);
        scene.boxes = (box*)malloc(sizeof(box) * scene.num_boxes);
        scene.tori = (torus*)malloc(sizeof(torus) * scene.num_tori);
        scene.cones = (cone*)malloc(sizeof(cone) * scene.num_cones);
        scene.octahedra = (octa*)malloc(sizeof(octa) * scene.num_octahedra);

        // allocate memory for the vectorized data layout
        // spheres
        scene.sphere_vecs.center_x = (float*)malloc(sizeof(float) * scene.num_spheres);
        scene.sphere_vecs.center_y = (float*)malloc(sizeof(float) * scene.num_spheres);
        scene.sphere_vecs.center_z = (float*)malloc(sizeof(float) * scene.num_spheres);
        scene.sphere_vecs.radius = (float*)malloc(sizeof(float) * scene.num_spheres);

        // boxes
        scene.box_vecs.bottom_left_x = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.bottom_left_y = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.bottom_left_z = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.extents_x = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.extents_y = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.extents_z = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.r = (float*)malloc(sizeof(float) * scene.num_boxes);

        scene.box_vecs.inv_rot[0][0] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[0][1] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[0][2] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[1][0] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[1][1] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[1][2] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[2][0] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[2][1] = (float*)malloc(sizeof(float) * scene.num_boxes);
        scene.box_vecs.inv_rot[2][2] = (float*)malloc(sizeof(float) * scene.num_boxes);

        // tori
        scene.torus_vecs.center_x = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.center_y = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.center_z = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.r1 = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.r2 = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.r = (float*)malloc(sizeof(float) * scene.num_tori);

        scene.torus_vecs.inv_rot[0][0] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[0][1] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[0][2] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[1][0] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[1][1] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[1][2] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[2][0] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[2][1] = (float*)malloc(sizeof(float) * scene.num_tori);
        scene.torus_vecs.inv_rot[2][2] = (float*)malloc(sizeof(float) * scene.num_tori);

        // cones
        scene.cone_vecs.center_x = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.center_y = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.center_z = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.r1 = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.r2 = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.height = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.r = (float*)malloc(sizeof(float) * scene.num_cones);

        scene.cone_vecs.inv_rot[0][0] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[0][1] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[0][2] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[1][0] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[1][1] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[1][2] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[2][0] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[2][1] = (float*)malloc(sizeof(float) * scene.num_cones);
        scene.cone_vecs.inv_rot[2][2] = (float*)malloc(sizeof(float) * scene.num_cones);

        // octahedra
        scene.octa_vecs.center_x = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.center_y = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.center_z = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.s = (float*)malloc(sizeof(float) * scene.num_octahedra);

        scene.octa_vecs.inv_rot[0][0] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[0][1] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[0][2] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[1][0] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[1][1] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[1][2] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[2][0] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[2][1] = (float*)malloc(sizeof(float) * scene.num_octahedra);
        scene.octa_vecs.inv_rot[2][2] = (float*)malloc(sizeof(float) * scene.num_octahedra);

        // second pass to actually load the shapes
        int sphere_idx = 0;
        int plane_idx = 0;
        int box_idx = 0;
        int torus_idx = 0;
        int cone_idx = 0;
        int octa_idx = 0;

        for (int i = 0; i < num_shapes; i++) {
            json current_shape = j["objects"][i];
            std::string current = current_shape["kind"].get<std::string>();
            if (current == "sphere") {
                scene.spheres[sphere_idx] = load_sphere(current_shape);

                // vectorized data layout
                scene.sphere_vecs.center_x[sphere_idx] = scene.spheres[sphere_idx].center.x;
                scene.sphere_vecs.center_y[sphere_idx] = scene.spheres[sphere_idx].center.y;
                scene.sphere_vecs.center_z[sphere_idx] = scene.spheres[sphere_idx].center.z;
                scene.sphere_vecs.radius[sphere_idx] = scene.spheres[sphere_idx].radius;

                sphere_idx++;

            } else if (current == "plane") {
                scene.planes[plane_idx++] = load_plane(current_shape);
            } else if (current == "box") {
                scene.boxes[box_idx] = load_box(current_shape);

                // vectorized data layout
                scene.box_vecs.bottom_left_x[box_idx] = scene.boxes[box_idx].bottom_left.x;
                scene.box_vecs.bottom_left_y[box_idx] = scene.boxes[box_idx].bottom_left.y;
                scene.box_vecs.bottom_left_z[box_idx] = scene.boxes[box_idx].bottom_left.z;

                scene.box_vecs.extents_x[box_idx] = scene.boxes[box_idx].extents.x;
                scene.box_vecs.extents_y[box_idx] = scene.boxes[box_idx].extents.y;
                scene.box_vecs.extents_z[box_idx] = scene.boxes[box_idx].extents.z;
                scene.box_vecs.r[box_idx] = scene.boxes[box_idx].r;

                scene.box_vecs.inv_rot[0][0][box_idx] = scene.boxes[box_idx].inv_rot.val[0][0];
                scene.box_vecs.inv_rot[0][1][box_idx] = scene.boxes[box_idx].inv_rot.val[0][1];
                scene.box_vecs.inv_rot[0][2][box_idx] = scene.boxes[box_idx].inv_rot.val[0][2];
                scene.box_vecs.inv_rot[1][0][box_idx] = scene.boxes[box_idx].inv_rot.val[1][0];
                scene.box_vecs.inv_rot[1][1][box_idx] = scene.boxes[box_idx].inv_rot.val[1][1];
                scene.box_vecs.inv_rot[1][2][box_idx] = scene.boxes[box_idx].inv_rot.val[1][2];
                scene.box_vecs.inv_rot[2][0][box_idx] = scene.boxes[box_idx].inv_rot.val[2][0];
                scene.box_vecs.inv_rot[2][1][box_idx] = scene.boxes[box_idx].inv_rot.val[2][1];
                scene.box_vecs.inv_rot[2][2][box_idx] = scene.boxes[box_idx].inv_rot.val[2][2];

                box_idx++;

            } else if (current == "torus") {
                scene.tori[torus_idx] = load_torus(current_shape);

                // vectorized data layout
                scene.torus_vecs.center_x[torus_idx] = scene.tori[torus_idx].center.x;
                scene.torus_vecs.center_y[torus_idx] = scene.tori[torus_idx].center.y;
                scene.torus_vecs.center_z[torus_idx] = scene.tori[torus_idx].center.z;
                scene.torus_vecs.r1[torus_idx] = scene.tori[torus_idx].r1;
                scene.torus_vecs.r2[torus_idx] = scene.tori[torus_idx].r2;
                scene.torus_vecs.r[torus_idx] = scene.tori[torus_idx].r;

                scene.torus_vecs.inv_rot[0][0][torus_idx] = scene.tori[torus_idx].inv_rot.val[0][0];
                scene.torus_vecs.inv_rot[0][1][torus_idx] = scene.tori[torus_idx].inv_rot.val[0][1];
                scene.torus_vecs.inv_rot[0][2][torus_idx] = scene.tori[torus_idx].inv_rot.val[0][2];
                scene.torus_vecs.inv_rot[1][0][torus_idx] = scene.tori[torus_idx].inv_rot.val[1][0];
                scene.torus_vecs.inv_rot[1][1][torus_idx] = scene.tori[torus_idx].inv_rot.val[1][1];
                scene.torus_vecs.inv_rot[1][2][torus_idx] = scene.tori[torus_idx].inv_rot.val[1][2];
                scene.torus_vecs.inv_rot[2][0][torus_idx] = scene.tori[torus_idx].inv_rot.val[2][0];
                scene.torus_vecs.inv_rot[2][1][torus_idx] = scene.tori[torus_idx].inv_rot.val[2][1];
                scene.torus_vecs.inv_rot[2][2][torus_idx] = scene.tori[torus_idx].inv_rot.val[2][2];

                torus_idx++;

            } else if (current == "cone") {
                scene.cones[cone_idx] = load_cone(current_shape);

                // vectorized data layout
                scene.cone_vecs.center_x[cone_idx] = scene.cones[cone_idx].center.x;
                scene.cone_vecs.center_y[cone_idx] = scene.cones[cone_idx].center.y;
                scene.cone_vecs.center_z[cone_idx] = scene.cones[cone_idx].center.z;
                scene.cone_vecs.r1[cone_idx] = scene.cones[cone_idx].r1;
                scene.cone_vecs.r2[cone_idx] = scene.cones[cone_idx].r2;
                scene.cone_vecs.height[cone_idx] = scene.cones[cone_idx].height;
                scene.cone_vecs.r[cone_idx] = scene.cones[cone_idx].r;

                scene.cone_vecs.inv_rot[0][0][cone_idx] = scene.cones[cone_idx].inv_rot.val[0][0];
                scene.cone_vecs.inv_rot[0][1][cone_idx] = scene.cones[cone_idx].inv_rot.val[0][1];
                scene.cone_vecs.inv_rot[0][2][cone_idx] = scene.cones[cone_idx].inv_rot.val[0][2];
                scene.cone_vecs.inv_rot[1][0][cone_idx] = scene.cones[cone_idx].inv_rot.val[1][0];
                scene.cone_vecs.inv_rot[1][1][cone_idx] = scene.cones[cone_idx].inv_rot.val[1][1];
                scene.cone_vecs.inv_rot[1][2][cone_idx] = scene.cones[cone_idx].inv_rot.val[1][2];
                scene.cone_vecs.inv_rot[2][0][cone_idx] = scene.cones[cone_idx].inv_rot.val[2][0];
                scene.cone_vecs.inv_rot[2][1][cone_idx] = scene.cones[cone_idx].inv_rot.val[2][1];
                scene.cone_vecs.inv_rot[2][2][cone_idx] = scene.cones[cone_idx].inv_rot.val[2][2];

                cone_idx++;

            } else if (current == "octahedron") {
                scene.octahedra[octa_idx] = load_octa(current_shape);

                // vectorized data layout
                scene.octa_vecs.center_x[octa_idx] = scene.octahedra[octa_idx].center.x;
                scene.octa_vecs.center_y[octa_idx] = scene.octahedra[octa_idx].center.y;
                scene.octa_vecs.center_z[octa_idx] = scene.octahedra[octa_idx].center.z;
                scene.octa_vecs.s[octa_idx] = scene.octahedra[octa_idx].s;

                scene.octa_vecs.inv_rot[0][0][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[0][0];
                scene.octa_vecs.inv_rot[0][1][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[0][1];
                scene.octa_vecs.inv_rot[0][2][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[0][2];
                scene.octa_vecs.inv_rot[1][0][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[1][0];
                scene.octa_vecs.inv_rot[1][1][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[1][1];
                scene.octa_vecs.inv_rot[1][2][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[1][2];
                scene.octa_vecs.inv_rot[2][0][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[2][0];
                scene.octa_vecs.inv_rot[2][1][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[2][1];
                scene.octa_vecs.inv_rot[2][2][octa_idx] = scene.octahedra[octa_idx].inv_rot.val[2][2];

                octa_idx++;
            } else {
                throw std::runtime_error("Unknown shape " + current);
            }
        }
    }

    // }}}
    void load_scene(std::string& input) {
        std::ifstream i;
        i.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try {
            i.open(input);
        } catch (std::system_error& e) {
            std::cerr << "Failed to open input file '" << input << "': " << strerror(errno) << std::endl;
            throw;
        }

        json j;
        i >> j;

        // load camera paremters
        scene.cam = load_camera(j);
        // load lights
        load_light(j);
        // don't have to save shapes since it uses static variable
        load_shapes(j);
    }

} // namespace impl::vec3
