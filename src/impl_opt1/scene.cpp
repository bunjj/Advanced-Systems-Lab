/**
 * Contains all functions related to the scene in the reference implementation.
 *
 * All other implementations will use this to construct their scenes.
 */
#include "impl_opt1/scene.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace impl::opt1 {

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
    float sphere_distance(const sphere sp, const vec from) {
        INS_INC(sphere);
        INS_ADD;
        return vec_length(vec_sub(sp.center, from)) - sp.radius;
    }

    vec sphere_normal(sphere s, vec pos) {
        // return vec_normalize(vec_sub(pos, s.center));

        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            sphere_distance(s, vec_add(pos, delta1)) - sphere_distance(s, vec_sub(pos, delta1)),
            sphere_distance(s, vec_add(pos, delta2)) - sphere_distance(s, vec_sub(pos, delta2)),
            sphere_distance(s, vec_add(pos, delta3)) - sphere_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    sphere make_sphere(float x, float y, float z, float r, vec color, float reflection, float shininess) {
        sphere s;
        s.center = {x, y, z};
        s.radius = r;
        m44 matrix = get_transf_matrix({x, y, z}, {0, 0, 0});
        s.inv_matrix = m44_inv(matrix);
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    sphere load_sphere(json& j) {
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
    float box_distance(const box b, const vec from) {
        INS_INC(box);
        vec pos = vec4_to_vec(m44_mul_vec(b.inv_matrix, vec4_from_point(from)));
        vec q = vec_sub(vec_abs(pos), b.extents);
        return FADD(vec_length(vec_max(q, 0)), min(0.0f, max(max(q.x, q.y), q.z)));
    }

    vec box_normal(box s, vec pos) {
        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            box_distance(s, vec_add(pos, delta1)) - box_distance(s, vec_sub(pos, delta1)),
            box_distance(s, vec_add(pos, delta2)) - box_distance(s, vec_sub(pos, delta2)),
            box_distance(s, vec_add(pos, delta3)) - box_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    box make_box(vec bottom_left, vec extents, vec rot, vec color, float reflection, float shininess) {
        box s;
        s.bottom_left = bottom_left;
        s.extents = extents;
        m44 matrix = get_transf_matrix(bottom_left, rot);
        s.inv_matrix = m44_inv(matrix);
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    box load_box(json& j) {
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

    float plane_distance(const plane p, const vec from) {
        INS_INC(plane);
        return vec_dot(p.normal, vec_sub(from, p.point));
    }

    vec plane_normal(plane s, vec pos) {
        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            plane_distance(s, vec_add(pos, delta1)) - plane_distance(s, vec_sub(pos, delta1)),
            plane_distance(s, vec_add(pos, delta2)) - plane_distance(s, vec_sub(pos, delta2)),
            plane_distance(s, vec_add(pos, delta3)) - plane_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    plane make_plane(vec normal, vec point, vec color, float reflection, float shininess) {
        plane s;
        s.normal = vec_normalize(normal);
        s.point = point;
        s.inv_matrix = identity;
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
    float torus_distance(const torus t, const vec from) {
        INS_INC(torus);
        vec pos = vec4_to_vec(m44_mul_vec(t.inv_matrix, vec4_from_point(from)));
        vec2 posxz = {pos.x, pos.z};
        INS_ADD;
        vec2 q = {vec2_length(posxz) - t.r1, pos.y};

        INS_ADD;
        return vec2_length(q) - t.r2;
    }

    vec torus_normal(torus s, vec pos) {
        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            torus_distance(s, vec_add(pos, delta1)) - torus_distance(s, vec_sub(pos, delta1)),
            torus_distance(s, vec_add(pos, delta2)) - torus_distance(s, vec_sub(pos, delta2)),
            torus_distance(s, vec_add(pos, delta3)) - torus_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    torus make_torus(vec center, float r1, float r2, vec rot, vec color, float reflection, float shininess) {
        torus s;
        s.center = center;
        s.r1 = r1;
        s.r2 = r2;
        m44 matrix = get_transf_matrix(center, rot);
        s.inv_matrix = m44_inv(matrix);
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    torus load_torus(json& j) {
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
    float cone_distance(const cone c, const vec from) {
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

    vec cone_normal(cone s, vec pos) {
        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            cone_distance(s, vec_add(pos, delta1)) - cone_distance(s, vec_sub(pos, delta1)),
            cone_distance(s, vec_add(pos, delta2)) - cone_distance(s, vec_sub(pos, delta2)),
            cone_distance(s, vec_add(pos, delta3)) - cone_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    cone make_cone(
        vec center, float r1, float r2, float height, vec rot, vec color, float reflection, float shininess) {
        cone s;
        s.center = center;
        s.r1 = r1;
        s.r2 = r2;
        s.height = height;
        m44 matrix = get_transf_matrix(center, rot);
        s.inv_matrix = m44_inv(matrix);
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
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

        vec color = load_vec(j["color"]);
        float reflection = j["reflection"];
        float shininess = j["shininess"];

        return make_cone(pos, r1, r2, height, rot, color, reflection, shininess);
    }

    // }}}

    // Octahedron {{{
    float octahedron_distance(const octa o, const vec from) {
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

    vec octahedron_normal(octa s, vec pos) {
        // Numerical approximation of the normal (for now)
        static const float delta = 10e-5;
        vec delta1 = {delta, 0, 0};
        vec delta2 = {0, delta, 0};
        vec delta3 = {0, 0, delta};

        INS_INC1(add, 3);
        vec normal = vec_normalize({
            octahedron_distance(s, vec_add(pos, delta1)) - octahedron_distance(s, vec_sub(pos, delta1)),
            octahedron_distance(s, vec_add(pos, delta2)) - octahedron_distance(s, vec_sub(pos, delta2)),
            octahedron_distance(s, vec_add(pos, delta3)) - octahedron_distance(s, vec_sub(pos, delta3)),
        });

        return normal;
    }

    octa make_octahedron(vec center, float s_param, vec rot, vec color, float reflection, float shininess) {
        octa s;
        s.center = center;
        s.s = s_param;
        m44 matrix = get_transf_matrix(center, rot);
        s.inv_matrix = m44_inv(matrix);
        s.color = color;
        s.reflection = reflection;
        s.shininess = shininess;
        return s;
    }

    octa load_octa(json& j) {
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
                scene.spheres[sphere_idx++] = load_sphere(current_shape);
            } else if (current == "plane") {
                scene.planes[plane_idx++] = load_plane(current_shape);
            } else if (current == "box") {
                scene.boxes[box_idx++] = load_box(current_shape);
            } else if (current == "torus") {
                scene.tori[torus_idx++] = load_torus(current_shape);
            } else if (current == "cone") {
                scene.cones[cone_idx++] = load_cone(current_shape);
            } else if (current == "octahedron") {
                scene.octahedra[octa_idx++] = load_octa(current_shape);
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

} // namespace impl::opt1
