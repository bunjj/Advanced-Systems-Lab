#include <dbg.h>
#include <inttypes.h>
#include <stdio.h>

#include <cmath>
#include <iostream>
#include <memory>

static float M_PI_F = M_PI;

#define TO_RAD(angle) ((angle) / 180.0f * M_PI_F)

struct vec {
    float x;
    float y;
    float z;
};

std::ostream& operator<<(std::ostream& out, const vec& v) {
    out << "{" << v.x << ", " << v.y << ", " << v.z << "}";
    return out;
}

struct vec4 {
    float x;
    float y;
    float z;
    float t;
};

vec4 vec4_init(float *values) {
    return {values[0], values[1], values[2], values[3]};
}

float vec4_dot(vec4 v1, vec4 v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.t * v2.t;
}

vec4 vec4_from_point(vec p) {
    return {p.x, p.y, p.z, 1};
}

vec4 vec4_from_dir(vec p) {
    return {p.x, p.y, p.z, 0};
}

vec vec4_to_vec(vec4 p) {
    return {p.x, p.y, p.z};
}


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
};

static const m44 identity = m44{.val = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}};

m44 get_rot_matrix(float x, float y, float z) {
    float c = TO_RAD(x);
    float b = TO_RAD(y);
    float a = TO_RAD(z);

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
 * Calculates m * v
 */
vec4 m44_mul_vec(m44 m, vec4 v) {
    float x = vec4_dot(vec4_init(m.val[0]), v);
    float y = vec4_dot(vec4_init(m.val[1]), v);
    float z = vec4_dot(vec4_init(m.val[2]), v);
    float t = vec4_dot(vec4_init(m.val[3]), v);

    return {x, y, z, t};
}

struct light {
    vec pos;
    vec color;
    float intensity;
};

struct hit {
    bool is_hit;
    float distance;
    int steps;
    vec color;
};

struct sphere {
    vec center;
    float radius;
};

float vec_norm(vec v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

float vec_length(vec v) { return sqrtf(vec_norm(v)); }

#define VEC_OP(v1, v2, OP) \
    vec { (v1).x OP(v2).x, (v1).y OP(v2).y, (v1).z OP(v2).z }

vec vec_add(vec v1, vec v2) { return VEC_OP(v1, v2, +); }

vec vec_sub(vec v1, vec v2) { return VEC_OP(v1, v2, -); }

vec vec_mul(vec v1, float factor) {
    return VEC_OP(v1, (vec{factor, factor, factor}), *);
}

vec vec_normalize(vec v) {
    float len = vec_length(v);
    return vec{v.x / len, v.y / len, v.z / len};
}

float vec_dot(vec v1, vec v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

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

static sphere spheres[] = {
    {{1, -5, -20}, 7},
    {{0, 1, -15}, 2},
    {{0, 120, -200}, 100},
};

int num_spheres = 3;

/*
 * Static light sources
 */
static light lights[] = {
    {{0, 10, 0}, {1.0, 0.9, 0.7}, 3000},
    {{20, 10, -15}, {0, 0, 1.0}, 6000},
    {{-20, 10, -15}, {0.847, 0.2588, 0.2588}, 6000},
    {{-100, 40, -80}, {0, 0.9, 0.7}, 90000},
    {{100, 40, -80}, {1.0, 0, 0.7}, 90000},
};

int num_lights = 5;

// max distance
static float D = 2048;
static float EPS = 0.001;

static bool sphere_trace_shadow(vec point, vec light_dir, float max_distance) {
    // TODO if we start with t = 0 this function causes some pixels to be black
    // because it erroneously detects a collision with the original object
    float t = EPS;

    while (t < max_distance) {
        vec pos = vec_add(point, vec_mul(light_dir, t));

        float min_distance = INFINITY;
        int sphere = -1;

        for (int k = 0; k < num_spheres; k++) {
            float distance =
                vec_length(vec_sub(pos, spheres[k].center)) - spheres[k].radius;

            if (distance < min_distance) {
                min_distance = distance;
                sphere = k;

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
        vec pos = vec_add(origin, vec_mul(dir, t));

        float min_distance = INFINITY;
        int sphere = -1;

        for (int k = 0; k < num_spheres; k++) {
            float distance =
                vec_length(vec_sub(pos, spheres[k].center)) - spheres[k].radius;

            if (distance < min_distance) {
                min_distance = distance;
                sphere = k;
            }
        }

        if (min_distance < EPS) {
            vec normal = vec_normalize(vec_sub(pos, spheres[sphere].center));
            vec color{0, 0, 0};

            for (int i = 0; i < num_lights; i++) {
                vec light_point = vec_sub(lights[i].pos, pos);

                if (vec_dot(light_point, normal) > 0) {
                    vec light_dir = vec_normalize(light_point);
                    // Squared distance from light to point
                    float dist_sq = vec_norm(light_point);

                    if (!sphere_trace_shadow(pos, light_dir, sqrt(dist_sq))) {
                        color =
                            vec_add(color, vec_mul(lights[i].color,
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

m44 get_camera_matrix(vec pos, float rot_x, float rot_y, float rot_z) {
    m44 camera_matrix = get_rot_matrix(rot_x, rot_y, rot_z);

    camera_matrix.val[0][3] = pos.x;
    camera_matrix.val[1][3] = pos.y;
    camera_matrix.val[2][3] = pos.z;

    return camera_matrix;
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

    vec camera_pos = {0, 0, 0};
    m44 camera_matrix = get_camera_matrix(camera_pos, 0, 0, 0);

    float aspect_ratio = static_cast<float>(width) / height;

    // Field of view in degrees
    float fov = 80;
    float fov_factor = tanf(TO_RAD(fov / 2));

    auto pixels = std::make_unique<float[]>(height * width * 3);

    vec origin = vec4_to_vec(m44_mul_vec(camera_matrix, vec4_from_point({0, 0, 0})));

    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            /*
             * Position of the pixel in camera space.
             *
             * We assume that the camera is  looking towards negative z and the
             * image plane is one unit away from the camera (z = -1 in this
             * case).
             */
            float x = (2 * (px + 0.5) / width - 1) * aspect_ratio * fov_factor;
            float y = (1 - 2 * (py + 0.5) / height) * fov_factor;
            float z = -1;

            // Direction in camera space.
            vec dir = vec_normalize({x, y, z});

            vec4 world_dir = m44_mul_vec(camera_matrix, vec4_from_dir(dir));

            auto h = sphere_trace(origin, vec4_to_vec(world_dir));

            vec color = h.is_hit ? h.color : vec{0.1, 0.1, 0.1};

            pixels[3 * (width * py + px)] = color.x;
            pixels[3 * (width * py + px) + 1] = color.y;
            pixels[3 * (width * py + px) + 2] = color.z;
        }
    }

    dump_image(width, height, pixels.get());

    return 0;
}
