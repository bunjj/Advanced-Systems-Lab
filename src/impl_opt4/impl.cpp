#include "impl_opt4/impl.hpp"

#include "impl_opt4/scene.hpp"

namespace impl::opt4 {
    // max distance
    static const float D = 100;
    static const float EPS = 0.001;

    struct hit {
        bool is_hit;
        float distance;
        int steps;
        vec color;
    };

    // https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    static float sphere_trace_softshadow(Ray r_world, float max_distance) {
        float t = EPS;

        float sharpness = 3.f; // sharpness of shadows
        float res = 1.f;
        while (t < max_distance) {
            vec p_world = trace_ray(r_world, t);

            float min_distance = INFINITY;

            // spheres
            for (int k = 0; k < scene.num_spheres; k++) {
                float distance = sphere_distance_short(scene.spheres[k], p_world, min_distance);

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

            // planes
            for (int k = 0; k < scene.num_planes; k++) {
                float distance = plane_distance(scene.planes[k], p_world);

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

            // boxes
            for (int k = 0; k < scene.num_boxes; k++) {
                box s = scene.boxes[k];
                vec p_obj = invtransform_point(s.inv_rot, s.bottom_left, p_world);
                float distance = box_distance_short(scene.boxes[k], p_obj, min_distance);

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

            // tori
            for (int k = 0; k < scene.num_tori; k++) {
                torus s = scene.tori[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);  
                float distance = torus_distance_short(scene.tori[k], p_obj, min_distance);

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

            // cones
            for (int k = 0; k < scene.num_cones; k++) {
                cone s = scene.cones[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);
                float distance = cone_distance_short(scene.cones[k], p_obj, min_distance);

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

            // octahedra
            for (int k = 0; k < scene.num_octahedra; k++) {
                octa s = scene.octahedra[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);
                float distance = octahedron_distance_short(scene.octahedra[k], p_obj, min_distance);

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
            res = min(res, sharpness * min_distance / t);
            t = FADD(t, min_distance);
        }
        return res;
    }

    static vec shade(vec normal, float shininess, float reflection, vec color, vec pos, vec dir, float t) {
        /* Prepare shading parameters */
        INS_MUL;
        float alpha = shininess;     // shininess parameter
        float ks = reflection * 0.4; // specular parameter
        float kd = 1.f;              // diffuse parameter
        float ka = 0.0075f;          // ambient parameter
        float sigma_a = 4e-6f;       // atmospheric absorbtion coeff

        vec wi;                        // incident direction
        vec wr;                        // reflected direction
        vec wo = vec_scale(dir, -1.f); // outgoing direction
        vec wn = normal;               // normal direction

        vec Li;             // incident Light
        vec Lo = {0, 0, 0}; // outgoing Light
        vec La = {0, 0, 0}; // ambient Light

        for (int i = 0; i < scene.num_lights; i++) {
            wi = vec_sub(scene.lights[i].pos, pos); // unnormalized

            INS_DIV;
            float dist2 = vec_dot2(wi);   // squared distance of light
            float dist = FSQRT(dist2);    // distance of light
            wi = vec_scale(wi, 1 / dist); // normalize incident direction

            // incoming light
            INS_INC1(mul, 2);
            INS_DIV;
            Li = vec_scale(scene.lights[i].emission, 1 / (4 * M_PI_F * dist2)); // incident light
            La = vec_add(La, Li); // incident light contributes to ambient light

            INS_CMP;
            if (vec_dot(wn, wi) > 0) {
                Ray shadowray = {pos, wi};
                float shadow = sphere_trace_softshadow(shadowray, dist);
                INS_CMP;
                if (shadow > EPS) {
                    Li = vec_scale(Li, shadow);

                    // diffuse
                    INS_MUL;
                    vec f_diffuse = vec_scale(color, kd * vec_dot(wn, wi)); // fraction of reflected light
                    Lo = vec_add(Lo, vec_mul(Li, f_diffuse));               // diffuse contribution to outgoing light

                    // specular
                    INS_INC1(mul, 2);
                    INS_POW;
                    wr = vec_sub(vec_scale(wn, 2 * vec_dot(wn, wi)), wi); // reflected direction
                    float f_specular =
                        ks * pow(max(0.f, vec_dot(wr, wo)), alpha); // fraction of reflected light TODO: normalization?
                    Lo = vec_add(Lo, vec_scale(Li, f_specular));    // specular contribution to outgoing light
                }
            }
        }

        vec f_ambient = vec_scale(color, ka);     // fraction of reflected ambient light
        Lo = vec_add(Lo, vec_mul(La, f_ambient)); // ambient contribution to outgoing light

        // atmospheric effect using exponential decay
        INS_INC1(mul, 3);
        INS_POW;
        Lo = vec_scale(Lo, powf(M_E, -sigma_a * t * t * t));

        return Lo;
    }

    static hit sphere_trace(Ray r_world) {
        float t = 0;
        int steps = 0;
        vec color;

        while (t < D) {
            vec p_world = trace_ray(r_world, t);

            float min_distance = INFINITY;

            // spheres
            for (int k = 0; k < scene.num_spheres; k++) {
                float distance = sphere_distance_short(scene.spheres[k], p_world, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        sphere s = scene.spheres[k];
                        vec normal = sphere_normal(s, p_world);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // planes
            for (int k = 0; k < scene.num_planes; k++) {
                float distance = plane_distance(scene.planes[k], p_world);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        plane s = scene.planes[k];
                        vec normal = plane_normal(s, p_world);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // boxes
            for (int k = 0; k < scene.num_boxes; k++) {
                box s = scene.boxes[k];
                vec p_obj = invtransform_point(s.inv_rot, s.bottom_left, p_world);
                float distance = box_distance_short(scene.boxes[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        box s = scene.boxes[k];
                        vec normal = box_normal(s, p_obj);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // tori
            for (int k = 0; k < scene.num_tori; k++) {
                torus s = scene.tori[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);  
                float distance = torus_distance_short(scene.tori[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        vec normal = torus_normal(s, p_obj);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // cones
            for (int k = 0; k < scene.num_cones; k++) {
                cone s = scene.cones[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);
                float distance = cone_distance_short(scene.cones[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        cone s = scene.cones[k];
                        vec normal = cone_normal(s, p_obj);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // octahedra
            for (int k = 0; k < scene.num_octahedra; k++) {
                octa s = scene.octahedra[k];
                vec p_obj = invtransform_point(s.inv_rot, s.center, p_world);
                float distance = octahedron_distance_short(scene.octahedra[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        octa s = scene.octahedra[k];
                        vec normal = octahedron_normal(s, p_obj);
                        color = shade(normal, s.shininess, s.reflection, s.color, p_world, r_world.d, t);
                        return {true, t, steps, color};
                    }
                }
            }

            t = FADD(t, min_distance);
            steps++;
        }

        return hit{false, t, steps, {0, 0, 0}};
    }

    void render_init(std::string input) {
        // convert scene from reference format to custom format
        // std::cout << "Reference scene already loaded (" << input << "), converting to required format" << std::endl;
        // from_ref_scene();

        // alternatively, we can just read the scene again from the json file
        std::cout << "Loading scene again" << std::endl;
        load_scene(input);
    }

    void render(int width, int height, float* pixels) {
        m33 camera_rotation = get_rot_matrix_33(scene.cam.rotation);
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
                Ray r_cam = {{0.f,0.f,0.f}, vec_normalize({x, y, z})};
                Ray r_world = transform_ray(camera_rotation, scene.cam.pos, r_cam);

                auto h = sphere_trace(r_world);

                vec color = h.is_hit ? h.color : vec{0, 0, 0};

                pixels[3 * (width * py + px)] = color.x;
                pixels[3 * (width * py + px) + 1] = color.y;
                pixels[3 * (width * py + px) + 2] = color.z;
            }
        }
    }

} // namespace impl::opt4
