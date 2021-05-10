#include "impl_ref/geometry.h"
#include "impl_ref/scene.hpp"
#include "impl_ref/impl.hpp"

#include "impl_normals/impl.hpp"


namespace impl::normals {   

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
    using impl::ref::vec;
    using impl::ref::vec_normalize;
    using impl::ref::vec_sub;
    using impl::ref::vec_add;
    using impl::ref::vec_div;
    using impl::ref::vec_abs;

    using impl::ref::m44;
    using impl::ref::vec4_from_point;
    using impl::ref::m44_mul_vec;
    using impl::ref::vec4_to_vec;
    using impl::ref::m44_rotate_only;

    // import shapes
    using impl::ref::scene;
    using impl::ref::shape;
    using impl::ref::sphere;
    using impl::ref::box;
    using impl::ref::plane;


    vec sphere_normal(const shape s, const vec pos) {
        sphere sp = *((sphere*)s.data);
        return vec_normalize(vec_sub(pos, sp.center));
    }

    vec box_normal(const shape s, const vec from) {
        box b = *((box*)s.data);
        vec pos = vec4_to_vec(m44_mul_vec(s.inv_matrix, vec4_from_point(from)));
        vec abs = vec_abs(pos);

        // transform into upper right quadrant
        vec sign = vec_div(pos, abs);
        vec q = vec_sub(abs, b.extents);

        // argmax(q.x, q.y, q.z)
        vec n_obj = {0,0,0};
        INS_INC1(cmp, 3);
        if (q.x > q.y && q.x > q.z && q.x > 0) {
            n_obj = {1,0,0};
        } else if (q.y > q.z && q.y > 0) {
            n_obj = {0,1,0};
        } else if (q.z > 0) {
            n_obj = {0,0,1};
        }

        // invert transform from upper right quadrant, before abs()
        n_obj = vec_mul(sign, n_obj);

        vec n_world = m44_rotate_only(s.matrix, n_obj);
        return n_world;
    }


    vec plane_normal(const shape s, const vec pos) {
        plane p = *((plane*)s.data);
        return p.normal;
    }


    // replace shape.normal() with new function
    void render_init() {

        for (int k = 0; k < scene.num_shapes; k++) {

            switch (scene.shapes[k].type) {

                case SHAPE_SPHERE: 
                    scene.shapes[k].normal = sphere_normal; break;

                case SHAPE_BOX:
                    scene.shapes[k].normal = box_normal; break; 

                case SHAPE_PLANE:
                    scene.shapes[k].normal = plane_normal; break;

                default: break;
            }
        }
    }
    
    // forward call to reference implementation
    void render(int width, int height, float* pixels) {
        return impl::ref::render(width, height, pixels);
    }

} // namespace impl::normals
