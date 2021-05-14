#include "impl_ref/geometry.h"
#include "impl_ref/scene.hpp"
#include "impl_ref/impl.hpp"

#include "impl_opt0/impl.hpp"


namespace impl::opt0 { 


    // forward declaration
    enum shape_type {
        SHAPE_SPHERE = 0,
        SHAPE_PLANE,
        SHAPE_BOX,
        SHAPE_TORUS,
        SHAPE_CONE,
        SHAPE_OCTA,
    };

    // import vector geometry functions
    using impl::ref::min;
    using impl::ref::clamp;

    using impl::ref::vec;
    using impl::ref::vec_normalize;
    using impl::ref::vec_sub;
    using impl::ref::vec_add;
    using impl::ref::vec_div;
    using impl::ref::vec_abs;

    using impl::ref::vec2;
    using impl::ref::vec2_scale;
    using impl::ref::vec2_length;

    using impl::ref::m44;
    using impl::ref::vec4_from_point;
    using impl::ref::m44_mul_vec;
    using impl::ref::vec4_to_vec;

    // import shapes
    using impl::ref::scene;
    using impl::ref::shape;
    using impl::ref::sphere;
    using impl::ref::box;
    using impl::ref::plane;
    using impl::ref::torus;
    using impl::ref::cone;
    using impl::ref::octa;

    /**
     * Calculates m[0:2][0:2] * v
     */
    static inline vec m44_rotate_only(m44 m, vec v) {
        INS_INC1(mul, 9);
        INS_INC1(add, 6);
        vec res = {0,0,0};
        res.x = m.val[0][0] * v.x + m.val[0][1] * v.y + m.val[0][2] * v.z;
        res.y = m.val[1][0] * v.x + m.val[1][1] * v.y + m.val[1][2] * v.z;
        res.z = m.val[2][0] * v.x + m.val[2][1] * v.y + m.val[2][2] * v.z;

        return res;
    }  

    vec sphere_normal(const shape s, const vec pos) {
        INS_INC(sphere_n);
        sphere sp = *((sphere*)s.data);
        return vec_normalize(vec_sub(pos, sp.center));
    }

    vec box_normal(const shape s, const vec from) {
        INS_INC(box_n);
        box b = *((box*)s.data);
        vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec sign = vec_div(pos, abs);
        vec q = vec_sub(abs, b.extents);

        // argmax(q.x, q.y, q.z)
        vec n_obj = {0,0,0};
        INS_INC1(cmp, 3);
        if (q.x > q.y && q.x > q.z && q.x > 0) {
            n_obj = {1,0,0};
            INS_INC1(cmp, 2);
        } else if (q.y > q.z && q.y > 0) {
            n_obj = {0,1,0};
            INS_INC1(cmp, 1);
        } else if (q.z > 0) {
            n_obj = {0,0,1};
        }

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        vec n_world = m44_rotate_only(s.matrix, n_obj);
        return n_world;
    }


    vec plane_normal(const shape s, const vec) {
        INS_INC(plane_n);

        plane p = *((plane*)s.data);
        return p.normal;
    }


    vec torus_normal(const shape s, const vec from) {
        INS_INC(torus_n);
        torus t = *((torus*)s.data);
        vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
        
        vec2 posxz = {pos.x, pos.z};
        
        INS_ADD;
        INS_DIV;
        posxz = vec2_scale(posxz, 1 - t.r1/vec2_length(posxz));
        vec q = {posxz.x, pos.y, posxz.y};
        vec n_obj = vec_normalize(q);

        vec n_world = m44_rotate_only(s.matrix, n_obj);
        return n_world;
    }


    vec cone_normal(const shape shap, const vec from) {
        INS_INC(cone_n);
        cone c = *((cone*)shap.data);
        vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));

        float r1 = c.r1;
        float r2 = c.r2;
        float h = c.height;

        // transform into rotation invariant subspace around y-axis
        vec2 posxz = {pos.x, pos.z};
        vec2 q = {vec2_length(posxz), pos.y};
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

        // invert transform from rotation invariant subspace
        vec n_obj = {0,0,0};
        INS_CMP;
        if (vec2_dot2(ca) < vec2_dot2(cb)){
            INS_DIV;
            posxz = vec2_scale(posxz, ca.x/vec2_length(posxz));
            n_obj = {posxz.x, ca.y, posxz.y};
        } else {
            INS_DIV;
            posxz = vec2_scale(posxz, cb.x/vec2_length(posxz));
            n_obj = {posxz.x, cb.y, posxz.y};
        }
        n_obj = vec_normalize(n_obj);

        vec n_world = m44_rotate_only(shap.matrix, n_obj);
        return n_world;
    }

    vec octahedron_normal(const shape shap, const vec from) {
        INS_INC(octa_n);
        vec pos = vec4_to_vec(m44_mul_vec(shap.inv_matrix, vec4_from_point(from)));

        // transform into upper right quadrant
        vec abs = vec_abs(pos);
        vec sign = vec_div(pos, abs);

        vec n_obj = {1,1,1};
        n_obj = vec_normalize(n_obj);

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        vec n_world = m44_rotate_only(shap.matrix, n_obj);
        return n_world;
    }


    // replace shape.normal() with new function
    void render_init(std::string) {

        for (int k = 0; k < scene.num_shapes; k++) {

            switch (scene.shapes[k].type) {

                case SHAPE_SPHERE: 
                    scene.shapes[k].normal = sphere_normal; break;

                case SHAPE_BOX:
                    scene.shapes[k].normal = box_normal; break; 

                case SHAPE_PLANE:
                    scene.shapes[k].normal = plane_normal; break;
                
                case SHAPE_TORUS:
                    scene.shapes[k].normal = torus_normal; break;
                
                case SHAPE_CONE:
                    scene.shapes[k].normal = cone_normal; break;
                
                case SHAPE_OCTA:
                    scene.shapes[k].normal = octahedron_normal; break;

                default: break;
            }
        }
    }
    
    // forward call to reference implementation
    void render(int width, int height, float* pixels) {
        return impl::ref::render(width, height, pixels);
    }

} // namespace impl::opt0
