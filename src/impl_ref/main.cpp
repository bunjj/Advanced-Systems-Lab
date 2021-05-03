#include "impl.hpp"
#include "impl_ref/scene.hpp"

namespace impl::ref {
    // max distance
    static const float D = 2048;
    static const float EPS = 0.001;

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

            for (int k = 0; k < scene.num_shapes; k++) {
                float distance = scene.shapes[k].distance(scene.shapes[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    if (min_distance <= EPS * t) {
                        return true;
                    }
                }
            }

            t = FADD(t, min_distance);
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

            for (int k = 0; k < scene.num_shapes; k++) {
                float distance = scene.shapes[k].distance(scene.shapes[k], pos);

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
                shape s = scene.shapes[shape_idx];
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
                vec color{0, 0, 0};

                vec diffuse = {0.f, 0.f, 0.f};
                vec specular = {0.f, 0.f, 0.f};
                for (int i = 0; i < scene.num_lights; i++) {
                    vec light_point = vec_sub(scene.lights[i].pos, pos);

                    INS_CMP;
                    if (vec_dot(light_point, normal) > 0) {
                        vec light_dir = vec_normalize(light_point);
                        // Squared distance from light to point
                        float dist_sq = vec_dot2(light_point);

                        if (!sphere_trace_shadow(pos, light_dir, FSQRT(dist_sq))) {
                            INS_INC1(mul, 2);
                            INS_DIV;
                            vec light_intensity =
                                vec_scale(scene.lights[i].color, scene.lights[i].intensity / (4 * M_PI_F * dist_sq));

                            // diffuse
                            diffuse =
                                vec_add(diffuse, vec_scale(light_intensity, max(0.f, vec_dot(normal, light_dir))));

                            // specular
                            // TODO: how to choose n?
                            float n = 4.f;
                            INS_MUL;
                            vec r =
                                vec_add(vec_scale(normal, 2 * vec_dot(normal, vec_scale(light_dir, -1.f))), light_dir);
                            INS_POW;
                            specular = vec_add(specular, vec_scale(light_intensity, pow(max(0.f, vec_dot(r, dir)), n)));
                        }
                    }
                }

                // TODO: not clear what the unit of the shininess parameter in the scenes is
                // TODO: is ks + kd == 1 really a requirement?
                INS_MUL;
                float ks = s.shininess * 0.01f; // specular parameter
                INS_ADD;
                float kd = 1 - ks; // diffuse parameter

                vec object_color = s.color;

                // scale object color s.t. intensity of most prominent color is 1
                float max_color = max(max(object_color.x, object_color.y), object_color.z);
                INS_DIV;
                float scale_factor = 1.f / max_color;
                object_color = vec_scale(object_color, scale_factor);

                INS_INC1(mul, 3);
                diffuse = {diffuse.x * object_color.x, diffuse.y * object_color.y, diffuse.z * object_color.z};
                color = vec_add(color, vec_add(vec_scale(diffuse, kd), vec_scale(specular, ks)));

                return {true, t, steps, color};
            }

            t = FADD(t, min_distance);
            steps++;
        }

        return hit{false, t, steps, {0, 0, 0}};
    }

    void render_init() {}

    void render(int width, int height, float* pixels) {
        m44 camera_matrix = get_transf_matrix(scene.cam.pos, scene.cam.rotation);
        INS_DIV;
        INS_DIV;
        INS_MUL;
        INS_TAN;
        float fov_factor = tanf(TO_RAD(scene.cam.fov / 2));

        INS_DIV;
        float aspect_ratio = static_cast<float>(width) / height;
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

                auto h = sphere_trace(scene.cam.pos, vec4_to_vec(world_dir));

                vec color = h.is_hit ? h.color : vec{0, 0, 0};

                pixels[3 * (width * py + px)] = color.x;
                pixels[3 * (width * py + px) + 1] = color.y;
                pixels[3 * (width * py + px) + 2] = color.z;
            }
        }
    }

} // namespace impl::ref
