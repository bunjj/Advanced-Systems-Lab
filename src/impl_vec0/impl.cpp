#include "impl_vec0/impl.hpp"

#include "impl_vec0/scene.hpp"

namespace impl::vec0 {
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
    static float sphere_trace_softshadow(vec point, vec light_dir, float max_distance) {
        float t = EPS;

        float sharpness = 3.f; // sharpness of shadows
        float res = 1.f;
        while (t < max_distance) {
            vec pos = vec_add(point, vec_scale(light_dir, t));

            float min_distance = INFINITY;

            int k;

            // spheres
            for (k = 0; k < scene.num_spheres - 7; k += 8) {

                float dists[8];

                // store square distances between sphere center and point in dists (for early termination)
                sphere_distance_vectorized(k, dists, scene.sphere_vecs.center_x, scene.sphere_vecs.center_y, scene.sphere_vecs.center_z, scene.sphere_vecs.radius, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            return 0.0f;
                        }
                    }
                }

            }

            // remaining sphere iterations
            for (; k < scene.num_spheres; k++) {
                float distance = sphere_distance(scene.spheres[k], pos);

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
                float distance = plane_distance(scene.planes[k], pos);

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
            for (k = 0; k < scene.num_boxes - 7; k += 8) {

                float dists[8];

                // compute first part of distance function and store intermediate results (for early termination)
                box_distance_vectorized(k, dists, scene.box_vecs.bottom_left_x, scene.box_vecs.bottom_left_y, scene.box_vecs.bottom_left_z, scene.box_vecs.extents_x, scene.box_vecs.extents_y, scene.box_vecs.extents_z, scene.box_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            return 0.0f;
                        }
                    }
                }
            }

            // remaining box iterations
            for (; k < scene.num_boxes; k++) {
                float distance = box_distance(scene.boxes[k], pos);

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
            for (k = 0; k < scene.num_tori - 7; k += 8) {

                float dists[8];

                torus_distance_vectorized(k, dists, scene.torus_vecs.center_x, scene.torus_vecs.center_y, scene.torus_vecs.center_z, scene.torus_vecs.r1, scene.torus_vecs.r2, scene.torus_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            return 0.0f;
                        }
                    }
                }

            }

            // remaining torus iterations
            for (; k < scene.num_tori; k++) {
                float distance = torus_distance(scene.tori[k], pos);

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
            for (k = 0; k < scene.num_cones - 7; k += 8) {

                float dists[8];

                cone_distance_vectorized(k, dists, scene.cone_vecs.center_x, scene.cone_vecs.center_y, scene.cone_vecs.center_z, scene.cone_vecs.r1, scene.cone_vecs.r2, scene.cone_vecs.height, scene.cone_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            return 0.0f;
                        }
                    }
                }
            }

            // remaining cone iterations
            for (; k < scene.num_cones; k++) {
                float distance = cone_distance(scene.cones[k], pos);

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
            for (k = 0; k < scene.num_octahedra - 7; k += 8) {

                float dists[8];

                octahedron_distance_vectorized(k, dists, scene.octa_vecs.center_x, scene.octa_vecs.center_y, scene.octa_vecs.center_z, scene.octa_vecs.s, scene.octa_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            return 0.0f;
                        }
                    }
                }
            }

            // remaining octahedra iterations
            for (; k < scene.num_octahedra; k++) {
                float distance = octahedron_distance(scene.octahedra[k], pos);

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
                float shadow = sphere_trace_softshadow(pos, wi, dist);
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

    static hit sphere_trace(vec origin, vec dir) {
        float t = 0;
        int steps = 0;
        vec color;

        while (t < D) {
            vec pos = vec_add(origin, vec_scale(dir, t));

            float min_distance = INFINITY;

            int k;

            // spheres
            for (k = 0; k < scene.num_spheres - 7; k += 8) {

                float dists[8];

                // store square distances between sphere center and point in dists (for early termination)
                sphere_distance_vectorized(k, dists, scene.sphere_vecs.center_x, scene.sphere_vecs.center_y, scene.sphere_vecs.center_z, scene.sphere_vecs.radius, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            sphere s = scene.spheres[k+i];
                            vec normal = sphere_normal(s, pos);
                            color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                            return {true, t, steps, color};
                        }
                    }
                }

            }

            // remaining sphere iterations
            for (; k < scene.num_spheres; k++) {
                float distance = sphere_distance(scene.spheres[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        sphere s = scene.spheres[k];
                        vec normal = sphere_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // planes
            for (int k = 0; k < scene.num_planes; k++) {
                float distance = plane_distance(scene.planes[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        plane s = scene.planes[k];
                        vec normal = plane_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // boxes
            for (k = 0; k < scene.num_boxes - 7; k += 8) {

                float dists[8];

                // compute first part of distance function and store intermediate results (for early termination)
                box_distance_vectorized(k, dists, scene.box_vecs.bottom_left_x, scene.box_vecs.bottom_left_y, scene.box_vecs.bottom_left_z, scene.box_vecs.extents_x, scene.box_vecs.extents_y, scene.box_vecs.extents_z, scene.box_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            box s = scene.boxes[k+i];
                            vec normal = box_normal(s, pos);
                            color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                            return {true, t, steps, color};
                        }
                    }
                }
            }

            // remaining box iterations
            for (; k < scene.num_boxes; k++) {
                float distance = box_distance(scene.boxes[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        box s = scene.boxes[k];
                        vec normal = box_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // tori
            for (k = 0; k < scene.num_tori - 7; k += 8) {

                float dists[8];

                torus_distance_vectorized(k, dists, scene.torus_vecs.center_x, scene.torus_vecs.center_y, scene.torus_vecs.center_z, scene.torus_vecs.r1, scene.torus_vecs.r2, scene.torus_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            torus s = scene.tori[k+i];
                            vec normal = torus_normal(s, pos);
                            color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                            return {true, t, steps, color};
                        }
                    }
                }

            }

            // remaining torus iterations
            for (; k < scene.num_tori; k++) {
                float distance = torus_distance(scene.tori[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        torus s = scene.tori[k];
                        vec normal = torus_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // cones
            for (k = 0; k < scene.num_cones - 7; k += 8) {

                float dists[8];

                // cone_distance_vectorized(k, dists, scene.cone_vecs.center_x, scene.cone_vecs.center_y, scene.cone_vecs.center_z, scene.cone_vecs.r1, scene.cone_vecs.r2, scene.cone_vecs.height, scene.cone_vecs.inv_rot, pos);

                cone_distance_vectorized(k, dists, scene.cone_vecs.center_x, scene.cone_vecs.center_y, scene.cone_vecs.center_z, scene.cone_vecs.r1, scene.cone_vecs.r2, scene.cone_vecs.height, scene.cone_vecs.inv_rot, pos);

                for (int i = 0; i < 8; i++) {
                    INS_CMP;
                    if (dists[i] < min_distance) {
                        min_distance = dists[i];

                        INS_CMP;
                        INS_MUL;
                        if (min_distance <= EPS * t) {
                            cone s = scene.cones[k+i];
                            vec normal = cone_normal(s, pos);
                            color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                            return {true, t, steps, color};
                        }
                    }
                }
            }

            // remaining cone iterations
            for (; k < scene.num_cones; k++) {
                float distance = cone_distance(scene.cones[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        cone s = scene.cones[k];
                        vec normal = cone_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                        return {true, t, steps, color};
                    }
                }
            }

            // octahedra
            for (k = 0; k < scene.num_octahedra - 7; k += 8) {

                float dists[8];

                octahedron_distance_vectorized(k, dists, scene.octa_vecs.center_x, scene.octa_vecs.center_y, scene.octa_vecs.center_z, scene.octa_vecs.s, scene.octa_vecs.inv_rot, pos);

                    for (int i = 0; i < 8; i++) {
                        INS_CMP;
                        if (dists[i] < min_distance) {
                            min_distance = dists[i];

                            INS_CMP;
                            INS_MUL;
                            if (min_distance <= EPS * t) {
                                octa s = scene.octahedra[k+i];
                                vec normal = octahedron_normal(s, pos);
                                color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
                                return {true, t, steps, color};
                            }
                        }
                    }
            }

            // remaining octahedron iterations
            for (; k < scene.num_octahedra; k++) {
                float distance = octahedron_distance(scene.octahedra[k], pos);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        octa s = scene.octahedra[k];
                        vec normal = octahedron_normal(s, pos);
                        color = shade(normal, s.shininess, s.reflection, s.color, pos, dir, t);
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

} // namespace impl::vec0
