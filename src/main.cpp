#include <assert.h>
#include <dbg.h>
#include <inttypes.h>
#include <stdio.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>

#include "geometry.h"
#include "instrument.h"
#include "loader.h"

flops_t flops_counter;

using json = nlohmann::json;

// Scene {{{
struct camera {
    float fov;
    vec pos;
    vec rotation;
};

struct light {
    vec pos;
    vec emission;
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
};

shape make_shape(vec color, float reflection, float shininess, const m44 matrix, distance_fun f, normal_fun n, void* data, size_t data_size) {
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


vec shape_normal(const shape s, const vec pos){
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
    return make_shape(color, reflection, shininess, get_transf_matrix({x, y, z}, {0, 0, 0}), sphere_distance, shape_normal, &s, sizeof(s));
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
    return make_shape(color, reflection, shininess, get_transf_matrix(bottom_left, rot), box_distance, shape_normal, &s, sizeof(s));
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
    return make_shape(color, reflection, shininess, get_transf_matrix(center, rot), torus_distance, shape_normal, &t, sizeof(t));
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

shape make_cone(vec center, float r1, float r2, float height, vec rot, vec color, float reflection, float shininess) {
    cone c = {center, r1, r2, height};
    return make_shape(color, reflection, shininess, get_transf_matrix(center, rot), cone_distance, shape_normal, &c, sizeof(c));
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
    return make_shape(color, reflection, shininess, get_transf_matrix(center, rot), octahedron_distance, shape_normal, &o, sizeof(o));
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

static shape* shapes;
static int num_shapes;
static light* lights;
static int num_lights;
// }}}

// Sphere Tracing {{{

// max distance
static float D = 100;
static float EPS = 0.001;

struct hit {
    bool is_hit;
    float distance;
    int steps;
    vec color;
};

static float sphere_trace_shadow(vec point, vec light_dir, float max_distance) {
    // TODO if we start with t = 0 this function causes some pixels to be black
    // because it erroneously detects a collision with the original object
    float t = EPS;

    while (t < max_distance) {
        vec pos = vec_add(point, vec_scale(light_dir, t));

        float min_distance = INFINITY;

        for (int k = 0; k < num_shapes; k++) {
            float distance = shapes[k].distance(shapes[k], pos);

            INS_CMP;
            if (distance < min_distance) {
                min_distance = distance;

                INS_CMP;
                INS_MUL;
                if (min_distance <= EPS * t) {
                    return 0.f;
                }
            }
        }

        t = FADD(t, min_distance);
    }

    return 1.f;
}

// https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
static float sphere_trace_softshadow(vec point, vec light_dir, float max_distance) {
    float t = EPS;

    float k = 3.f;
    float res = 1.f;
    while (t < max_distance) {
        vec pos = vec_add(point, vec_scale(light_dir, t));

        float min_distance = INFINITY;

        for (int k = 0; k < num_shapes; k++) {
            float distance = shapes[k].distance(shapes[k], pos);

            INS_CMP;
            if (distance < min_distance) {
                min_distance = distance;

                INS_CMP;
                INS_MUL;
                if (min_distance <= EPS * t) {
                    return 0.0f;
                }
            }
        }

        INS_MUL;
        INS_DIV;
        res = min( res, k * min_distance / t);
        t = FADD(t, min_distance);
    }
    return res;
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

            INS_CMP;
            if (distance < min_distance) {
                min_distance = distance;
                shape_idx = k;

                INS_CMP;
                if (min_distance <= EPS) {
                    break;
                }
            }
        }

        INS_CMP;
        if (min_distance <= EPS) {
            shape s = shapes[shape_idx];

            /* Prepare shading parameters */
            INS_MUL;
            float alpha = s.shininess;      // shininess parameter
            float ks = s.reflection * 0.4;  // specular parameter 
            float kd = 1.f;                 // diffuse parameter
            float ka = 0.0075f;             // ambient parameter

            vec wi;                         // incident direction
            vec wr;                         // reflected direction
            vec wo = vec_scale(dir, -1.f);  // outgoing direction
            vec wn = s.normal(s, pos);      // normal direction

            vec Li;                         // incident Light
            vec Lo = {0, 0, 0};             // outgoing Light
            vec La = {0, 0, 0};             // ambient Light

            for (int i = 0; i < num_lights; i++) {

                wi = vec_sub(lights[i].pos, pos);   // unnormalized 
                
                INS_DIV;
                float dist2 = vec_dot2(wi);     // squared distance of light
                float dist = FSQRT(dist2);      // distance of light
                wi = vec_scale(wi, 1/dist);     // normalize incident direction

                // incoming light 
                INS_INC1(mul, 2);
                INS_DIV;
                Li = vec_scale(lights[i].emission, 1 / (4 * M_PI_F * dist2)); // incident light
                La = vec_add(La, Li);   // incident light contributes to ambient light

                INS_CMP;
                if (vec_dot(wn, wi) > 0) {    

                    float shadow = sphere_trace_softshadow(pos, wi, dist);
                    if (shadow > EPS) {

                        Li = vec_scale(Li, shadow);
                        
                        // diffuse
                        INS_MUL;
                        vec f_diffuse = vec_scale(s.color, kd * vec_dot(wn, wi)); // fraction of reflected light
                        Lo = vec_add(Lo, vec_mul(Li, f_diffuse)); // diffuse contribution to outgoing light

                        // specular
                        INS_INC1(mul, 2);
                        INS_POW;
                        wr = vec_sub( vec_scale(wn, 2 * vec_dot(wn, wi) ), wi); // reflected direction
                        float f_specular = ks * pow(max(0.f, vec_dot(wr, wo)), alpha);  // fraction of reflected light TODO: normalization?
                        Lo = vec_add(Lo, vec_scale(Li, f_specular)); // specular contribution to outgoing light
                    }
                }
            }
            
            vec f_ambient = vec_scale(s.color, ka); // fraction of reflected ambient light
            Lo = vec_add(Lo, vec_mul(La, f_ambient)); // ambient contribution to outgoing light

            // fog
            INS_INC1(mul, 3);
            INS_POW; // TODO: count exponentials separately?
            Lo = vec_scale(Lo, std::exp(-4e-6*t*t*t ));

            return {true, t, steps, Lo};
        }

        t = FADD(t, min_distance);
        steps++;
    }

    return hit{false, t, steps, {0, 0, 0}};
}

// }}}

static void dump_image_ldr(std::ostream& out, int width, int height, const float* pixels, float exposure) {
    out << "P6\n" << width << " " << height << "\n255\n";

    // TODO: gamma correction?
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                unsigned char channel = std::clamp(exposure * pixels[3 * (width * j + i) + k], 0.f, 1.f) * 255.f;
                out << channel;
            }
        }
    }
}

static void dump_image_hdr(std::ostream& out, int width, int height, const float* pixels, float exposure) {
    out << "PF\n" << width << " " << height << "\n-1.0\n";

    for (int j = height-1; j > -1; j--) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < 3; k++) {
                float channel = exposure * pixels[3 * (width * j + i) + k];
                out.write((char *) &channel, sizeof(float));
            }
        }
    }   
}

static void dump_image(std::ostream& out, int width, int height, const float* pixels, bool hdr=false){
    
    /* TODO: Possible Strategies to determine exposure:
        - user defined by exposure parameter in camera
        - choose exposure empirically 
        - map max_val to 1
        - map 95th percentile to 1
    */

    float exposure = 1250.f; // 

   /*
    float max_val = 0.f;
    for (int i = 0; i < width * height * 3; i++){
        max_val = std::max(max_val, pixels[i]);
    }
    float exposure = 1.f / max_val * 1.76; 
    */

    if (hdr) {
        // store in binary .pfm format
        dump_image_hdr(out, width, height, pixels, exposure);
    } else {
        // store in binary .ppm format
        dump_image_ldr(out, width, height, pixels, exposure);
    }
}

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

static std::vector<light> load_light(json& j) {
    std::vector<light> lights;
    json light = j["pointlight"];

    if (light.is_array()) {
        for (auto& l : light) {
            lights.push_back(load_single_light(l));
        }
    } else {
        lights.push_back(load_single_light(light));
    }

    return lights;
}

static void load_shapes(json& j) {
    num_shapes = j["objects"].size();
    shapes = (shape*)malloc(sizeof(shape) * num_shapes);

    for (int i = 0; i < num_shapes; i++) {
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
            std::cerr << "Unknown shape " << current << std::endl;
            break;
        }
        shapes[i] = new_shape;
    }
}

// }}}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input> <output>\n", argv[0]);
        return EXIT_FAILURE;
    }

    std::string input = argv[1];
    std::string output = argv[2];

    std::ifstream i;
    i.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        i.open(input);
    } catch (std::system_error& e) {
        std::cerr << "Failed to open input file '" << input << "': " << strerror(errno) << std::endl;
        return EXIT_FAILURE;
    }

    json j;
    i >> j;

    ins_rst();

    // Height of the resulting image in pixels
    int height = 1080;
    int width = 1920;

    // load camera paremters
    camera cam = load_camera(j);
    // load light (as defined in scene, not with intensity parameter)
    std::vector<light> lights_vector = load_light(j);

    lights = lights_vector.data();
    num_lights = lights_vector.size();

    // don't have to save shapes since it uses static variable
    load_shapes(j);

    m44 camera_matrix = get_transf_matrix(cam.pos, cam.rotation);

    float aspect_ratio = static_cast<float>(width) / height;

    float fov_factor = tanf(TO_RAD(cam.fov / 2));

    auto pixels = std::make_unique<float[]>(height * width * 3);

    vec origin = vec4_to_vec(m44_mul_vec(camera_matrix, vec4_from_point({0, 0, 0})));

    ins_dump("Setup");
    ins_rst();

    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            /*
             * Position of the pixel in camera space.
             *
             * We assume that the camera is looking towards positive z and the
             * image plane is one unit away from the camera (z = -1 in this
             * case).
             */
            INS_INC1(add, 4);
            INS_INC1(mul, 5);
            INS_INC1(div, 2);
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

    std::ofstream o;
    o.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    
    bool hdr;
    std::string fileformat = output.substr(output.find_last_of("."));
    if (fileformat.compare(".pfm") == 0) {
        hdr = true;
    } else if (fileformat.compare(".ppm") == 0) {
        hdr = false;
    } else {
        std::cerr << "Unsupported file format '" << fileformat << "', defaulting to '.pfm'" << std::endl;
        hdr = true;
        output.append(".pfm");
    }

    try {
        o.open(output);
    } catch (std::system_error& e) {
        std::cerr << "Failed to open output file '" << output << "': " << strerror(errno) << std::endl;
        return EXIT_FAILURE;
    }

    dump_image(o, width, height, pixels.get(), hdr);

    ins_dump(NULL); // TODO: dump before dumping image, in case something goes wrong?

    return 0;
}
