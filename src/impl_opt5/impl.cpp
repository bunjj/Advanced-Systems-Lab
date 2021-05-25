#include "impl_opt5/impl.hpp"

#include "impl_opt5/scene.hpp"

namespace impl::opt5 {
    // max distance
    static const float D = 100;
    static const float EPS = 0.001;

    static Ray* r_boxes;
    static Ray* r_tori;
    static Ray* r_cones;
    static Ray* r_octahedra;

    static Ray* r_boxes_shade;
    static Ray* r_tori_shade;
    static Ray* r_cones_shade;
    static Ray* r_octahedra_shade;

    // https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    static float sphere_trace_softshadow(Ray r_world, float max_distance) {
        
        /* Precompute ray in all object coordinates */
        for (int k = 0; k < scene.num_boxes; k++) {
            r_boxes_shade[k] = invtransform_ray(scene.boxes[k].inv_rot, scene.boxes[k].bottom_left, r_world);
        }
        for (int k = 0; k < scene.num_tori; k++) {
            r_tori_shade[k] = invtransform_ray(scene.tori[k].inv_rot, scene.tori[k].center, r_world);
        }         
        for (int k = 0; k < scene.num_cones; k++) {
            r_cones_shade[k] = invtransform_ray(scene.cones[k].inv_rot, scene.cones[k].center, r_world);
        }
        for (int k = 0; k < scene.num_octahedra; k++) {
            r_octahedra_shade[k] = invtransform_ray(scene.octahedra[k].inv_rot, scene.octahedra[k].center, r_world);
        }

        /* Actual Shadow Tracing */
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
                vec p_obj = trace_ray(r_boxes_shade[k], t);
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
                vec p_obj = trace_ray(r_tori_shade[k], t);
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
                vec p_obj = trace_ray(r_cones_shade[k], t);
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
                vec p_obj = trace_ray(r_octahedra_shade[k], t);
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

    static vec sphere_trace(vec d) {
        /* Precompute ray direction in all object coordinates */  
        for (int k = 0; k < scene.num_boxes; k++) {
            r_boxes[k].d = m33_mul_vec(scene.boxes[k].inv_rot, d);
        }
        for (int k = 0; k < scene.num_tori; k++) {
            r_tori[k].d = m33_mul_vec(scene.tori[k].inv_rot, d);
        }         
        for (int k = 0; k < scene.num_cones; k++) {
            r_cones[k].d = m33_mul_vec(scene.cones[k].inv_rot, d);
        }
        for (int k = 0; k < scene.num_octahedra; k++) {
            r_octahedra[k].d = m33_mul_vec(scene.octahedra[k].inv_rot, d);
        }

        /* Actual Sphere Tracing */
        float t = 0;

        while (t < D) {
            vec p_world = vec_add(scene.cam.pos, vec_scale(d, t));

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
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
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
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                    }
                }
            }

            // boxes
            for (int k = 0; k < scene.num_boxes; k++) {
                vec p_obj = trace_ray(r_boxes[k], t);
                float distance = box_distance_short(scene.boxes[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        box s = scene.boxes[k];
                        vec normal = box_normal(s, p_obj);
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                    }
                }
            }

            // tori
            for (int k = 0; k < scene.num_tori; k++) {
                vec p_obj = trace_ray(r_tori[k], t);
                float distance = torus_distance_short(scene.tori[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        torus s = scene.tori[k];
                        vec normal = torus_normal(s, p_obj);
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                    }
                }
            }

            // cones
            for (int k = 0; k < scene.num_cones; k++) {
                vec p_obj = trace_ray(r_cones[k], t);
                float distance = cone_distance_short(scene.cones[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        cone s = scene.cones[k];
                        vec normal = cone_normal(s, p_obj);
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                    }
                }
            }

            // octahedra
            for (int k = 0; k < scene.num_octahedra; k++) {
                vec p_obj = trace_ray(r_octahedra[k], t);
                float distance = octahedron_distance_short(scene.octahedra[k], p_obj, min_distance);

                INS_CMP;
                if (distance < min_distance) {
                    min_distance = distance;

                    INS_CMP;
                    INS_MUL;
                    if (min_distance <= EPS * t) {
                        octa s = scene.octahedra[k];
                        vec normal = octahedron_normal(s, p_obj);
                        return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                    }
                }
            }

            t = FADD(t, min_distance);
        }

        return {0, 0, 0};
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

        r_boxes = (struct Ray*) malloc(sizeof(Ray) * scene.num_boxes);
        r_tori = (struct Ray*) malloc(sizeof(Ray) * scene.num_tori);
        r_cones = (struct Ray*) malloc(sizeof(Ray) * scene.num_cones);
        r_octahedra = (struct Ray*) malloc(sizeof(Ray) * scene.num_octahedra);

        r_boxes_shade = (struct Ray*) malloc(sizeof(Ray) * scene.num_boxes);
        r_tori_shade = (struct Ray*) malloc(sizeof(Ray) * scene.num_tori);
        r_cones_shade = (struct Ray*) malloc(sizeof(Ray) * scene.num_cones);
        r_octahedra_shade = (struct Ray*) malloc(sizeof(Ray) * scene.num_octahedra);

        /* Precompute camera ray origin in object coordinates */
        for (int k = 0; k < scene.num_boxes; k++) {
            r_boxes[k].o = invtransform_point(scene.boxes[k].inv_rot, scene.boxes[k].bottom_left, scene.cam.pos);
        }
        for (int k = 0; k < scene.num_tori; k++) {
            r_tori[k].o = invtransform_point(scene.tori[k].inv_rot, scene.tori[k].center, scene.cam.pos);
        }         
        for (int k = 0; k < scene.num_cones; k++) {
            r_cones[k].o = invtransform_point(scene.cones[k].inv_rot, scene.cones[k].center, scene.cam.pos);
        }
        for (int k = 0; k < scene.num_octahedra; k++) {
            r_octahedra[k].o = invtransform_point(scene.octahedra[k].inv_rot, scene.octahedra[k].center, scene.cam.pos);
        }

        INS_DIV;
        float aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
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
                vec d = m33_mul_vec(camera_rotation, vec_normalize({x, y, z}));

                vec color = sphere_trace(d);

                memcpy(pixels + (3 * (width * py + px)), &color, 12);
            }
        }
    }
} // namespace impl::opt5
