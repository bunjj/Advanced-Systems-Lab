#include "impl_vec4/impl.hpp"

#include "impl_vec4/scene.hpp"

namespace impl::vec4 {
    // max distance
    static const float D = 100.f;
    static const float EPS = 0.001f;
    static const float PI4 = 4 * M_PI_F;
    static const float M_E_F = M_E;

    static const float sharpness = 3.f; // sharpness of shadows

    static float* boxes_o_x;
    static float* boxes_o_y;
    static float* boxes_o_z;

    static float* boxes_d_x;
    static float* boxes_d_y;
    static float* boxes_d_z;

    static float* boxes_o_shade_x;
    static float* boxes_o_shade_y;
    static float* boxes_o_shade_z;

    static float* tori_o_x;
    static float* tori_o_y;
    static float* tori_o_z;

    static float* tori_d_x;
    static float* tori_d_y;
    static float* tori_d_z;

    static float* tori_o_shade_x;
    static float* tori_o_shade_y;
    static float* tori_o_shade_z;

    static float* cones_o_x;
    static float* cones_o_y;
    static float* cones_o_z;

    static float* cones_d_x;
    static float* cones_d_y;
    static float* cones_d_z;

    static float* cones_o_shade_x;
    static float* cones_o_shade_y;
    static float* cones_o_shade_z;

    static float* octahedra_o_x;
    static float* octahedra_o_y;
    static float* octahedra_o_z;

    static float* octahedra_d_x;
    static float* octahedra_d_y;
    static float* octahedra_d_z;

    static float* octahedra_o_shade_x;
    static float* octahedra_o_shade_y;
    static float* octahedra_o_shade_z;

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
        int k = 0;
        for (; k < scene.num_boxes - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.box_vecs.inv_rot, scene.box_vecs.bottom_left_x,
                scene.box_vecs.bottom_left_y, scene.box_vecs.bottom_left_z, r_world.o);
            _mm256_store_ps(boxes_o_shade_x + k, pos.x);
            _mm256_store_ps(boxes_o_shade_y + k, pos.y);
            _mm256_store_ps(boxes_o_shade_z + k, pos.z);

            vec256 dir = m33_mul_vec_vectorized(k, scene.box_vecs.inv_rot, r_world.d);
            _mm256_store_ps(boxes_d_x + k, dir.x);
            _mm256_store_ps(boxes_d_y + k, dir.y);
            _mm256_store_ps(boxes_d_z + k, dir.z);
        }

        for (; k < scene.num_boxes; k++) {
            r_boxes_shade[k] = invtransform_ray(scene.boxes[k].inv_rot, scene.boxes[k].bottom_left, r_world);
            boxes_o_shade_x[k] = r_boxes_shade[k].o.x;
            boxes_o_shade_y[k] = r_boxes_shade[k].o.y;
            boxes_o_shade_z[k] = r_boxes_shade[k].o.z;
            boxes_d_x[k] = r_boxes_shade[k].d.x;
            boxes_d_y[k] = r_boxes_shade[k].d.y;
            boxes_d_z[k] = r_boxes_shade[k].d.z;
        }

        k = 0;
        for (; k < scene.num_tori - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.torus_vecs.inv_rot, scene.torus_vecs.center_x,
                scene.torus_vecs.center_y, scene.torus_vecs.center_z, r_world.o);
            _mm256_store_ps(tori_o_shade_x + k, pos.x);
            _mm256_store_ps(tori_o_shade_y + k, pos.y);
            _mm256_store_ps(tori_o_shade_z + k, pos.z);
            vec256 dir = m33_mul_vec_vectorized(k, scene.torus_vecs.inv_rot, r_world.d);
            _mm256_store_ps(tori_d_x + k, dir.x);
            _mm256_store_ps(tori_d_y + k, dir.y);
            _mm256_store_ps(tori_d_z + k, dir.z);
        }

        for (; k < scene.num_tori; k++) {
            r_tori_shade[k] = invtransform_ray(scene.tori[k].inv_rot, scene.tori[k].center, r_world);
            tori_o_shade_x[k] = r_tori_shade[k].o.x;
            tori_o_shade_y[k] = r_tori_shade[k].o.y;
            tori_o_shade_z[k] = r_tori_shade[k].o.z;
            tori_d_x[k] = r_tori_shade[k].d.x;
            tori_d_y[k] = r_tori_shade[k].d.y;
            tori_d_z[k] = r_tori_shade[k].d.z;
        }

        k = 0;
        for (; k < scene.num_cones - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.cone_vecs.inv_rot, scene.cone_vecs.center_x,
                scene.cone_vecs.center_y, scene.cone_vecs.center_z, r_world.o);
            _mm256_store_ps(cones_o_shade_x + k, pos.x);
            _mm256_store_ps(cones_o_shade_y + k, pos.y);
            _mm256_store_ps(cones_o_shade_z + k, pos.z);
            vec256 dir = m33_mul_vec_vectorized(k, scene.cone_vecs.inv_rot, r_world.d);
            _mm256_store_ps(cones_d_x + k, dir.x);
            _mm256_store_ps(cones_d_y + k, dir.y);
            _mm256_store_ps(cones_d_z + k, dir.z);
        }

        for (; k < scene.num_cones; k++) {
            r_cones_shade[k] = invtransform_ray(scene.cones[k].inv_rot, scene.cones[k].center, r_world);
            cones_o_shade_x[k] = r_cones_shade[k].o.x;
            cones_o_shade_y[k] = r_cones_shade[k].o.y;
            cones_o_shade_z[k] = r_cones_shade[k].o.z;
            cones_d_x[k] = r_cones_shade[k].d.x;
            cones_d_y[k] = r_cones_shade[k].d.y;
            cones_d_z[k] = r_cones_shade[k].d.z;
        }

        k = 0;
        for (; k < scene.num_octahedra - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.octa_vecs.inv_rot, scene.octa_vecs.center_x,
                scene.octa_vecs.center_y, scene.octa_vecs.center_z, r_world.o);
            _mm256_store_ps(octahedra_o_shade_x + k, pos.x);
            _mm256_store_ps(octahedra_o_shade_y + k, pos.y);
            _mm256_store_ps(octahedra_o_shade_z + k, pos.z);
            vec256 dir = m33_mul_vec_vectorized(k, scene.octa_vecs.inv_rot, r_world.d);
            _mm256_store_ps(octahedra_d_x + k, dir.x);
            _mm256_store_ps(octahedra_d_y + k, dir.y);
            _mm256_store_ps(octahedra_d_z + k, dir.z);
        }

        for (; k < scene.num_octahedra; k++) {
            r_octahedra_shade[k] = invtransform_ray(scene.octahedra[k].inv_rot, scene.octahedra[k].center, r_world);
            octahedra_o_shade_x[k] = r_octahedra_shade[k].o.x;
            octahedra_o_shade_y[k] = r_octahedra_shade[k].o.y;
            octahedra_o_shade_z[k] = r_octahedra_shade[k].o.z;
            octahedra_d_x[k] = r_octahedra_shade[k].d.x;
            octahedra_d_y[k] = r_octahedra_shade[k].d.y;
            octahedra_d_z[k] = r_octahedra_shade[k].d.z;
        }

        /* Actual Shadow Tracing */
        float t = EPS;

        float res = 1.f;
        while (t < max_distance) {
            vec p_world = trace_ray(r_world, t);

            float min_distance = INFINITY;

            int k;

            // spheres
            for (k = 0; k < scene.num_spheres - 7; k += 8) {
                __m256 dists;

                int not_terminate_early =
                    sphere_distance_short_vectorized(k, &dists, scene.sphere_vecs.center_x, scene.sphere_vecs.center_y,
                        scene.sphere_vecs.center_z, scene.sphere_vecs.radius, p_world, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    if (dist < min_distance) {
                        min_distance = dist;
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
            for (k = 0; k < scene.num_boxes - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(
                    k, boxes_o_shade_x, boxes_o_shade_y, boxes_o_shade_z, boxes_d_x, boxes_d_y, boxes_d_z, t);
                __m256 dists;

                int not_terminate_early = box_distance_short_vectorized(k, &dists, scene.box_vecs.extents_x,
                    scene.box_vecs.extents_y, scene.box_vecs.extents_z, scene.box_vecs.r, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    if (dist < min_distance) {
                        min_distance = dist;
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
            for (k = 0; k < scene.num_tori - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(
                    k, tori_o_shade_x, tori_o_shade_y, tori_o_shade_z, tori_d_x, tori_d_y, tori_d_z, t);
                __m256 dists;

                int not_terminate_early = torus_distance_short_vectorized(
                    k, &dists, scene.torus_vecs.r1, scene.torus_vecs.r2, scene.torus_vecs.r, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    if (dist < min_distance) {
                        min_distance = dist;
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
            for (k = 0; k < scene.num_cones - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(
                    k, cones_o_shade_x, cones_o_shade_y, cones_o_shade_z, cones_d_x, cones_d_y, cones_d_z, t);
                __m256 dists;

                int not_terminate_early =
                    cone_distance_short_vectorized(k, &dists, scene.cone_vecs.r1, scene.cone_vecs.r2,
                        scene.cone_vecs.height, scene.cone_vecs.r, scene.cone_vecs.k2d2inv, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    if (dist < min_distance) {
                        min_distance = dist;
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
            for (k = 0; k < scene.num_octahedra - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(k, octahedra_o_shade_x, octahedra_o_shade_y, octahedra_o_shade_z,
                    octahedra_d_x, octahedra_d_y, octahedra_d_z, t);
                __m256 dists;

                int not_terminate_early =
                    octahedron_distance_short_vectorized(k, &dists, scene.octa_vecs.s, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    if (dist < min_distance) {
                        min_distance = dist;
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

    static const float kd = 1.f;             // diffuse parameter
    static const float ka = 0.0075f;         // ambient parameter
    static const float neg_sigma_a = -4e-6f; // atmospheric absorbtion coeff

    static vec shade(vec normal, float shininess, float reflection, vec color, vec pos, vec dir, float t) {
        /* Prepare shading parameters */
        INS_MUL;
        float ks = reflection * 0.4; // specular parameter

        vec wo = vec_scale(dir, -1.f); // outgoing direction

        vec Lo = {0, 0, 0}; // outgoing Light
        vec La = {0, 0, 0}; // ambient Light

        for (int i = 0; i < scene.num_lights; i++) {
            vec wi = vec_sub(scene.lights[i].pos, pos); // incident direction unnormalized

            INS_DIV;
            float dist2 = vec_dot2(wi);   // squared distance of light
            float dist = FSQRT(dist2);    // distance of light
            wi = vec_scale(wi, 1 / dist); // normalize incident direction

            // incoming light
            INS_MUL;
            INS_DIV;
            vec Li = vec_scale(scene.lights[i].emission, 1 / (PI4 * dist2)); // incident light
            La = vec_add(La, Li); // incident light contributes to ambient light

            INS_CMP;
            if (vec_dot(normal, wi) > 0) {
                float shadow = sphere_trace_softshadow({pos, wi}, dist);
                INS_CMP;
                if (shadow > EPS) {
                    Li = vec_scale(Li, shadow);

                    // diffuse
                    INS_MUL;
                    vec f_diffuse = vec_scale(color, kd * vec_dot(normal, wi)); // fraction of reflected light
                    Lo = vec_add(Lo, vec_mul(Li, f_diffuse)); // diffuse contribution to outgoing light

                    // specular
                    INS_INC1(mul, 2);
                    INS_POW;

                    vec wr = vec_sub(vec_scale(normal, 2 * vec_dot(normal, wi)), wi); // reflected direction
                    float f_specular = ks * pow(max(0.f, vec_dot(wr, wo)),
                                                shininess);      // fraction of reflected light TODO: normalization?
                    Lo = vec_add(Lo, vec_scale(Li, f_specular)); // specular contribution to outgoing light
                }
            }
        }

        vec f_ambient = vec_scale(color, ka);     // fraction of reflected ambient light
        Lo = vec_add(Lo, vec_mul(La, f_ambient)); // ambient contribution to outgoing light

        // atmospheric effect using exponential decay
        INS_INC1(mul, 3);
        INS_POW;
        Lo = vec_scale(Lo, powf(M_E_F, neg_sigma_a * t * t * t));

        return Lo;
    }

    static vec sphere_trace(vec d) {
        /* Precompute ray direction in all object coordinates */
        int k = 0;
        for (; k < scene.num_boxes - 7; k += 8) {
            vec256 dir = m33_mul_vec_vectorized(k, scene.box_vecs.inv_rot, d);
            _mm256_store_ps(boxes_d_x + k, dir.x);
            _mm256_store_ps(boxes_d_y + k, dir.y);
            _mm256_store_ps(boxes_d_z + k, dir.z);
        }

        for (; k < scene.num_boxes; k++) {
            r_boxes[k].d = m33_mul_vec(scene.boxes[k].inv_rot, d);
            boxes_d_x[k] = r_boxes[k].d.x;
            boxes_d_y[k] = r_boxes[k].d.y;
            boxes_d_z[k] = r_boxes[k].d.z;
        }

        k = 0;
        for (; k < scene.num_tori - 7; k += 8) {
            vec256 dir = m33_mul_vec_vectorized(k, scene.torus_vecs.inv_rot, d);
            _mm256_store_ps(tori_d_x + k, dir.x);
            _mm256_store_ps(tori_d_y + k, dir.y);
            _mm256_store_ps(tori_d_z + k, dir.z);
        }

        for (; k < scene.num_tori; k++) {
            r_tori[k].d = m33_mul_vec(scene.tori[k].inv_rot, d);
            tori_d_x[k] = r_tori[k].d.x;
            tori_d_y[k] = r_tori[k].d.y;
            tori_d_z[k] = r_tori[k].d.z;
        }

        k = 0;
        for (; k < scene.num_cones - 7; k += 8) {
            vec256 dir = m33_mul_vec_vectorized(k, scene.cone_vecs.inv_rot, d);
            _mm256_store_ps(cones_d_x + k, dir.x);
            _mm256_store_ps(cones_d_y + k, dir.y);
            _mm256_store_ps(cones_d_z + k, dir.z);
        }

        for (; k < scene.num_cones; k++) {
            r_cones[k].d = m33_mul_vec(scene.cones[k].inv_rot, d);
            cones_d_x[k] = r_cones[k].d.x;
            cones_d_y[k] = r_cones[k].d.y;
            cones_d_z[k] = r_cones[k].d.z;
        }

        k = 0;
        for (; k < scene.num_octahedra - 7; k += 8) {
            vec256 dir = m33_mul_vec_vectorized(k, scene.octa_vecs.inv_rot, d);
            _mm256_store_ps(octahedra_d_x + k, dir.x);
            _mm256_store_ps(octahedra_d_y + k, dir.y);
            _mm256_store_ps(octahedra_d_z + k, dir.z);
        }

        for (; k < scene.num_octahedra; k++) {
            r_octahedra[k].d = m33_mul_vec(scene.octahedra[k].inv_rot, d);
            octahedra_d_x[k] = r_octahedra[k].d.x;
            octahedra_d_y[k] = r_octahedra[k].d.y;
            octahedra_d_z[k] = r_octahedra[k].d.z;
        }

        float t = 0;
        while (t < D) {
            vec p_world = vec_add(scene.cam.pos, vec_scale(d, t));

            float min_distance = INFINITY;

            int k;

            /*
             * Temporarily stores 8 distances returned from distance function.
             * Only needed when we need to figure out which shape was hit.
             */
            static float dists_tmp[8] __attribute__((aligned(32)));

            // spheres
            for (k = 0; k < scene.num_spheres - 7; k += 8) {
                __m256 dists;

                int not_terminate_early =
                    sphere_distance_short_vectorized(k, &dists, scene.sphere_vecs.center_x, scene.sphere_vecs.center_y,
                        scene.sphere_vecs.center_z, scene.sphere_vecs.radius, p_world, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    INS_MUL;
                    if (dist > EPS * t) {
                        min_distance = dist;
                    } else {
                        _mm256_store_ps(dists_tmp, dists);
                        for (int i = 0; i < 8; i++) {
                            INS_CMP;
                            if (dists_tmp[i] == dist) {
                                sphere s = scene.spheres[k + i];
                                vec normal = sphere_normal(s, p_world);
                                return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                            }
                        }
                    }
                }
            }

            // remaining sphere iterations
            for (; k < scene.num_spheres; k++) {
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
            for (k = 0; k < scene.num_boxes - 7; k += 8) {
                vec256 p_obj =
                    trace_ray_vectorized(k, boxes_o_x, boxes_o_y, boxes_o_z, boxes_d_x, boxes_d_y, boxes_d_z, t);
                __m256 dists;

                int not_terminate_early = box_distance_short_vectorized(k, &dists, scene.box_vecs.extents_x,
                    scene.box_vecs.extents_y, scene.box_vecs.extents_z, scene.box_vecs.r, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    INS_MUL;
                    if (dist > EPS * t) {
                        min_distance = dist;
                    } else {
                        _mm256_store_ps(dists_tmp, dists);
                        for (int i = 0; i < 8; i++) {
                            INS_CMP;
                            if (dists_tmp[i] == dist) {
                                box s = scene.boxes[k + i];
                                vec normal = box_normal(s, vec256_get(p_obj, i));
                                return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                            }
                        }
                    }
                }
            }

            // remaining box iterations
            for (; k < scene.num_boxes; k++) {
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
            for (k = 0; k < scene.num_tori - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(k, tori_o_x, tori_o_y, tori_o_z, tori_d_x, tori_d_y, tori_d_z, t);
                __m256 dists;

                int not_terminate_early = torus_distance_short_vectorized(
                    k, &dists, scene.torus_vecs.r1, scene.torus_vecs.r2, scene.torus_vecs.r, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    INS_MUL;
                    if (dist > EPS * t) {
                        min_distance = dist;
                    } else {
                        _mm256_store_ps(dists_tmp, dists);
                        for (int i = 0; i < 8; i++) {
                            INS_CMP;
                            if (dists_tmp[i] == dist) {
                                torus s = scene.tori[k + i];
                                vec normal = torus_normal(s, vec256_get(p_obj, i));
                                return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                            }
                        }
                    }
                }
            }

            // remaining torus iterations
            for (; k < scene.num_tori; k++) {
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
            for (k = 0; k < scene.num_cones - 7; k += 8) {
                vec256 p_obj =
                    trace_ray_vectorized(k, cones_o_x, cones_o_y, cones_o_z, cones_d_x, cones_d_y, cones_d_z, t);
                __m256 dists;

                int not_terminate_early =
                    cone_distance_short_vectorized(k, &dists, scene.cone_vecs.r1, scene.cone_vecs.r2,
                        scene.cone_vecs.height, scene.cone_vecs.r, scene.cone_vecs.k2d2inv, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    INS_MUL;
                    if (dist > EPS * t) {
                        min_distance = dist;
                    } else {
                        _mm256_store_ps(dists_tmp, dists);
                        for (int i = 0; i < 8; i++) {
                            INS_CMP;
                            if (dists_tmp[i] == dist) {
                                cone s = scene.cones[k + i];
                                vec normal = cone_normal(s, vec256_get(p_obj, i));
                                return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                            }
                        }
                    }
                }
            }

            // remaining cone iterations
            for (; k < scene.num_cones; k++) {
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
            for (k = 0; k < scene.num_octahedra - 7; k += 8) {
                vec256 p_obj = trace_ray_vectorized(
                    k, octahedra_o_x, octahedra_o_y, octahedra_o_z, octahedra_d_x, octahedra_d_y, octahedra_d_z, t);
                __m256 dists;

                int not_terminate_early =
                    octahedron_distance_short_vectorized(k, &dists, scene.octa_vecs.s, p_obj, min_distance);

                if (not_terminate_early) {
                    float dist = min_vectorized(dists);

                    INS_CMP;
                    INS_MUL;
                    if (dist > EPS * t) {
                        min_distance = dist;
                    } else {
                        _mm256_store_ps(dists_tmp, dists);
                        for (int i = 0; i < 8; i++) {
                            INS_CMP;
                            if (dists_tmp[i] == dist) {
                                octa s = scene.octahedra[k + i];
                                vec normal = octahedron_normal(s, vec256_get(p_obj, i));
                                return shade(normal, s.shininess, s.reflection, s.color, p_world, d, t);
                            }
                        }
                    }
                }
            }

            // remaining octahedron iterations
            for (; k < scene.num_octahedra; k++) {
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
        load_scene(input);
    }

    void render(int width, int height, float* pixels) {
        m33 camera_rotation = get_rot_matrix_33(scene.cam.rotation);
        INS_DIV;
        INS_DIV;
        INS_MUL;
        INS_TAN;
        float fov_factor = tanf(TO_RAD(scene.cam.fov / 2));

        r_boxes = (struct Ray*)malloc(sizeof(Ray) * scene.num_boxes);
        r_tori = (struct Ray*)malloc(sizeof(Ray) * scene.num_tori);
        r_cones = (struct Ray*)malloc(sizeof(Ray) * scene.num_cones);
        r_octahedra = (struct Ray*)malloc(sizeof(Ray) * scene.num_octahedra);

        r_boxes_shade = (struct Ray*)malloc(sizeof(Ray) * scene.num_boxes);
        r_tori_shade = (struct Ray*)malloc(sizeof(Ray) * scene.num_tori);
        r_cones_shade = (struct Ray*)malloc(sizeof(Ray) * scene.num_cones);
        r_octahedra_shade = (struct Ray*)malloc(sizeof(Ray) * scene.num_octahedra);

        boxes_o_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_o_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_o_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);

        boxes_d_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_d_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_d_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);

        boxes_o_shade_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_o_shade_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);
        boxes_o_shade_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_boxes);

        tori_o_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_o_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_o_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);

        tori_d_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_d_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_d_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);

        tori_o_shade_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_o_shade_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);
        tori_o_shade_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_tori);

        cones_o_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_o_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_o_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);

        cones_d_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_d_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_d_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);

        cones_o_shade_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_o_shade_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);
        cones_o_shade_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_cones);

        octahedra_o_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_o_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_o_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);

        octahedra_d_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_d_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_d_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);

        octahedra_o_shade_x = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_o_shade_y = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);
        octahedra_o_shade_z = (float*)aligned_alloc_wrapper(32, 4 * scene.num_octahedra);

        /* Precompute camera ray origin in object coordinates */
        int k = 0;
        for (; k < scene.num_boxes - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.box_vecs.inv_rot, scene.box_vecs.bottom_left_x,
                scene.box_vecs.bottom_left_y, scene.box_vecs.bottom_left_z, scene.cam.pos);
            _mm256_store_ps(boxes_o_x + k, pos.x);
            _mm256_store_ps(boxes_o_y + k, pos.y);
            _mm256_store_ps(boxes_o_z + k, pos.z);
        }

        for (; k < scene.num_boxes; k++) {
            r_boxes[k].o = invtransform_point(scene.boxes[k].inv_rot, scene.boxes[k].bottom_left, scene.cam.pos);
        }

        k = 0;
        for (; k < scene.num_tori - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.torus_vecs.inv_rot, scene.torus_vecs.center_x,
                scene.torus_vecs.center_y, scene.torus_vecs.center_z, scene.cam.pos);
            _mm256_store_ps(tori_o_x + k, pos.x);
            _mm256_store_ps(tori_o_y + k, pos.y);
            _mm256_store_ps(tori_o_z + k, pos.z);
        }

        for (; k < scene.num_tori; k++) {
            r_tori[k].o = invtransform_point(scene.tori[k].inv_rot, scene.tori[k].center, scene.cam.pos);
        }

        k = 0;
        for (; k < scene.num_cones - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.cone_vecs.inv_rot, scene.cone_vecs.center_x,
                scene.cone_vecs.center_y, scene.cone_vecs.center_z, scene.cam.pos);
            _mm256_store_ps(cones_o_x + k, pos.x);
            _mm256_store_ps(cones_o_y + k, pos.y);
            _mm256_store_ps(cones_o_z + k, pos.z);
        }
        for (; k < scene.num_cones; k++) {
            r_cones[k].o = invtransform_point(scene.cones[k].inv_rot, scene.cones[k].center, scene.cam.pos);
        }

        k = 0;
        for (; k < scene.num_octahedra - 7; k += 8) {
            vec256 pos = invtransform_point_vectorized(k, scene.octa_vecs.inv_rot, scene.octa_vecs.center_x,
                scene.octa_vecs.center_y, scene.octa_vecs.center_z, scene.cam.pos);
            _mm256_store_ps(octahedra_o_x + k, pos.x);
            _mm256_store_ps(octahedra_o_y + k, pos.y);
            _mm256_store_ps(octahedra_o_z + k, pos.z);
        }

        for (; k < scene.num_octahedra; k++) {
            r_octahedra[k].o = invtransform_point(scene.octahedra[k].inv_rot, scene.octahedra[k].center, scene.cam.pos);
        }

        float fwidth = width;
        float fheight = height;

        INS_DIV;
        float aspect_ratio = fwidth / fheight;

        INS_MUL;
        float fov_aspect = aspect_ratio * fov_factor;

        for (int py = 0; py < height; py++) {
            for (int px = 0; px < width; px++) {
                /*
                 * Position of the pixel in camera space.
                 *
                 * We assume that the camera is looking towards positive z and the
                 * image plane is one unit away from the camera (z = -1 in this
                 * case).
                 */
                INS_INC1(add, 2);
                INS_INC1(mul, 2);
                INS_INC1(div, 2);

                /*
                 * 2 * (px + 0.5) == 2 * px + 1
                 *
                 * Now we only use integer operations
                 */
                float px_1 = static_cast<float>(2 * px + 1);
                float py_1 = static_cast<float>(2 * py + 1);

                float x = (px_1 / fwidth - 1) * fov_aspect;
                float y = (1 - py_1 / fheight) * fov_factor;
                float z = 1;

                // Direction in camera space.
                vec d = m33_mul_vec(camera_rotation, vec_normalize({x, y, z}));

                vec color = sphere_trace(d);

                memcpy(pixels + (3 * (width * py + px)), &color, 12);
            }
        }
    }
} // namespace impl::vec4
