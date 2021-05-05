/**
 * Contains all functions related to the scene in the reference implementation.
 *
 * All other implementations will use this to construct their scenes.
 */
#include "impl_ref/scene.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace impl::ref {

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

    shape make_shape(vec color, float reflection, float shininess, const m44 matrix, distance_fun f, normal_fun n,
        void* data, size_t data_size) {
        shape shap;
        shap.distance = f;
        shap.normal = n;
        shap.matrix = matrix;
        shap.inv_matrix = m44_inv(matrix);
        shap.color = color;
        shap.reflection = reflection;
        shap.shininess = shininess;
        memcpy(&shap.data, data, data_size);
        return shap;
    }

    vec shape_normal(const shape s, const vec pos) {
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        // Some shapes can calculate this directly
        vec normal = vec_normalize({
            s.distance(s, vec_add(pos, delta1)) - s.distance(s, vec_sub(pos, delta1)),
            s.distance(s, vec_add(pos, delta2)) - s.distance(s, vec_sub(pos, delta2)),
            s.distance(s, vec_add(pos, delta3)) - s.distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    // Sphere {{{
    float sphere_distance(const shape s, const vec from) {
        INS_INC(sphere);
        sphere sp = *((sphere*)s.data);
        INS_ADD;
        return vec_length(vec_sub(sp.center, from)) - sp.radius;
    }

    vec sphere_normal(sphere s, vec pos) {
        return vec_normalize(vec_sub(pos, s.center));
    }

    shape make_sphere(float x, float y, float z, float r, vec color, float reflection, float shininess) {
        sphere s = {{x, y, z}, r};
        return make_shape(color, reflection, shininess, get_transf_matrix({x, y, z}, {0, 0, 0}), sphere_distance,
            shape_normal, &s, sizeof(s));
    }

    shape load_sphere(json& j) {
        float r;
        vec pos = load_pos(j);
        r = j["params"]["radius"];
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];
        return make_sphere(pos.x, pos.y, pos.z, r, color, reflection, shininess);
    }

    // }}}

    // Box {{{
    float box_distance(const shape s, const vec from) {
        INS_INC(box);
        box b = *((box*)s.data);
        vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
        vec q = vec_sub(vec_abs(pos), b.extents);
        return FADD(vec_length(vec_max(q, 0)), min(0.0f, max(max(q.x, q.y), q.z)));
    }

    shape make_box(vec bottom_left, vec extents, vec rot, vec color, float reflection, float shininess) {
        box s = {bottom_left, extents};
        return make_shape(color, reflection, shininess, get_transf_matrix(bottom_left, rot), box_distance, shape_normal,
            &s, sizeof(s));
    }

    shape load_box(json& j) {
        vec pos = load_pos(j);
        vec extents = load_vec(j["params"]["extents"]);
        vec rot = load_rot(j);
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];
        return make_box(pos, extents, rot, color, reflection, shininess);
    }

    // }}}

    // Plane {{{

    float plane_distance(const shape s, const vec from) {
        INS_INC(plane);
        plane p = *((plane*)s.data);
        return vec_dot(p.normal, vec_sub(from, p.point));
    }

    shape make_plane(vec normal, vec point, vec color, float reflection, float shininess) {
        plane p = {vec_normalize(normal), point};
        return make_shape(color, reflection, shininess, identity, plane_distance, shape_normal, &p, sizeof(p));
    }

    shape load_plane(json& j) {
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
    float torus_distance(const shape s, const vec from) {
        INS_INC(torus);
        torus t = *((torus*)s.data);
        vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
        vec2 posxz = {pos.x, pos.z};
        INS_ADD;
        vec2 q = {vec2_length(posxz) - t.r1, pos.y};

        INS_ADD;
        return vec2_length(q) - t.r2;
    }

    shape make_torus(vec center, float r1, float r2, vec rot, vec color, float reflection, float shininess) {
        torus t = {center, r1, r2};
        return make_shape(
            color, reflection, shininess, get_transf_matrix(center, rot), torus_distance, shape_normal, &t, sizeof(t));
    }

    shape load_torus(json& j) {
        float r1, r2;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        r1 = j["params"]["r1"];
        r2 = j["params"]["r2"];

        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_torus(pos, r1, r2, rot, color, reflection, shininess);
    }

    // }}}

    // Cone {{{
    float cone_distance(const shape shap, const vec from) {
        INS_INC(cone);
        cone c = *((cone*)shap.data);
        vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));

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

    shape make_cone(
        vec center, float r1, float r2, float height, vec rot, vec color, float reflection, float shininess) {
        cone c = {center, r1, r2, height};
        return make_shape(
            color, reflection, shininess, get_transf_matrix(center, rot), cone_distance, shape_normal, &c, sizeof(c));
    }

    shape load_cone(json& j) {
        // TODO: adjust cone definition according to mail!
        float r1, r2, height;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        r1 = j["params"][0];
        r2 = j["params"][1];
        height = j["params"][2];

        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_cone(pos, r1, r2, height, rot, color, reflection, shininess);
    }

    // }}}

    // Octahedron {{{
    float octahedron_distance(const shape shap, const vec from) {
        INS_INC(octa);
        octa o = *((octa*)shap.data);
        vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));
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

    shape make_octahedron(vec center, float s, vec rot, vec color, float reflection, float shininess) {
        octa o = {center, s};
        return make_shape(color, reflection, shininess, get_transf_matrix(center, rot), octahedron_distance,
            shape_normal, &o, sizeof(o));
    }

    shape load_octa(json& j) {
        float s;
        vec pos = load_pos(j);
        vec rot = load_rot(j);

        s = j["params"]["s"];
        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_octahedron(pos, s, rot, color, reflection, shininess);
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
        scene.num_shapes = j["objects"].size();
        scene.shapes = (shape*)malloc(sizeof(shape) * scene.num_shapes);

        for (int i = 0; i < scene.num_shapes; i++) {
            shape new_shape;
            json current_shape = j["objects"][i];
            std::string current = current_shape["kind"].get<std::string>();
            if (current == "sphere") {
                new_shape = load_sphere(current_shape);
            } else if (current == "plane") {
                new_shape = load_plane(current_shape);
            } else if (current == "box") {
                new_shape = load_box(current_shape);
            } else if (current == "torus") {
                new_shape = load_torus(current_shape);
            } else if (current == "cone") {
                new_shape = load_cone(current_shape);
            } else if (current == "octahedron") {
                new_shape = load_octa(current_shape);
            } else {
                throw std::runtime_error("Unknown shape " + current);
            }
            scene.shapes[i] = new_shape;
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

} // namespace impl::ref
