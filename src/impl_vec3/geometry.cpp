#include "impl_vec3/geometry.h"

namespace impl::vec3 {
    float M_PI_F = M_PI;

    std::ostream& operator<<(std::ostream& out, const vec& v) {
        out << "{" << v.x << ", " << v.y << ", " << v.z << "}";
        return out;
    }

    m33 get_rot_matrix_33(vec rot) {
        INS_INC1(div, 3);
        INS_INC1(mul, 3);
        float c = TO_RAD(rot.x);
        float b = TO_RAD(rot.y);
        float a = TO_RAD(rot.z);

        float ca = cosf(a);
        float cb = cosf(b);
        float cc = cosf(c);
        float sa = sinf(a);
        float sb = sinf(b);
        float sc = sinf(c);

        m33 m = identity_33;

        INS_INC1(add, 5);
        INS_INC1(mul, 16);

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
     * Computes only inverse for rotations, i.e computes the transpose. have to cleanup the code
     */
    m33 m33_inv(m33 m) {
        m33 m_t = identity_33;
        m_t.val[0][0] = m.val[0][0];
        m_t.val[0][1] = m.val[1][0];
        m_t.val[0][2] = m.val[2][0];
        m_t.val[1][0] = m.val[0][1];
        m_t.val[1][1] = m.val[1][1];
        m_t.val[1][2] = m.val[2][1];
        m_t.val[2][0] = m.val[0][2];
        m_t.val[2][1] = m.val[1][2];
        m_t.val[2][2] = m.val[2][2];

        return m_t;
    }
} // namespace impl::vec3
