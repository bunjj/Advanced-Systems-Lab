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

    // import shapes
    using impl::ref::scene;
    using impl::ref::shape;
    using impl::ref::sphere;
    using impl::ref::plane;


    vec sphere_normal(const shape s, const vec pos) {
        sphere sp = *((sphere*)s.data);
        return vec_normalize(vec_sub(pos, sp.center));
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
                    std::cout<< "worked" << std::endl;
                    scene.shapes[k].normal = sphere_normal; break;

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
