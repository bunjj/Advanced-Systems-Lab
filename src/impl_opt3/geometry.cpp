#include "impl_opt3/geometry.h"

namespace impl::opt3 {
    float M_PI_F = M_PI;

    std::ostream& operator<<(std::ostream& out, const vec& v) {
        out << "{" << v.x << ", " << v.y << ", " << v.z << "}";
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const m44& m) {
        int width = 8;
        out << std::fixed << std::setprecision(2) << std::right << std::endl;
        // clang-format off
    out << "┌ " << std::setw(width) << m.val[0][0] << " " << std::setw(width) << m.val[0][1] << " " << std::setw(width) << m.val[0][2] << " " << std::setw(width) << m.val[0][3] << " ┐" << std::endl;
    out << "│ " << std::setw(width) << m.val[1][0] << " " << std::setw(width) << m.val[1][1] << " " << std::setw(width) << m.val[1][2] << " " << std::setw(width) << m.val[1][3] << " │" << std::endl;
    out << "│ " << std::setw(width) << m.val[2][0] << " " << std::setw(width) << m.val[2][1] << " " << std::setw(width) << m.val[2][2] << " " << std::setw(width) << m.val[2][3] << " │" << std::endl;
    out << "└ " << std::setw(width) << m.val[3][0] << " " << std::setw(width) << m.val[3][1] << " " << std::setw(width) << m.val[3][2] << " " << std::setw(width) << m.val[3][3] << " ┘" << std::endl;
        // clang-format on
        return out;
    }

    m44 m44_inv(m44 m) {
        float a11 = m.val[0][0];
        float a12 = m.val[0][1];
        float a13 = m.val[0][2];
        float a14 = m.val[0][3];
        float a21 = m.val[1][0];
        float a22 = m.val[1][1];
        float a23 = m.val[1][2];
        float a24 = m.val[1][3];
        float a31 = m.val[2][0];
        float a32 = m.val[2][1];
        float a33 = m.val[2][2];
        float a34 = m.val[2][3];
        float a41 = m.val[3][0];
        float a42 = m.val[3][1];
        float a43 = m.val[3][2];
        float a44 = m.val[3][3];

        INS_INC1(add, 18);
        INS_INC1(mul, 36);

        float A2323 = a33 * a44 - a34 * a43;
        float A1323 = a32 * a44 - a34 * a42;
        float A1223 = a32 * a43 - a33 * a42;
        float A0323 = a31 * a44 - a34 * a41;
        float A0223 = a31 * a43 - a33 * a41;
        float A0123 = a31 * a42 - a32 * a41;
        float A2313 = a23 * a44 - a24 * a43;
        float A1313 = a22 * a44 - a24 * a42;
        float A1213 = a22 * a43 - a23 * a42;
        float A2312 = a23 * a34 - a24 * a33;
        float A1312 = a22 * a34 - a24 * a32;
        float A1212 = a22 * a33 - a23 * a32;
        float A0313 = a21 * a44 - a24 * a41;
        float A0213 = a21 * a43 - a23 * a41;
        float A0312 = a21 * a34 - a24 * a31;
        float A0212 = a21 * a33 - a23 * a31;
        float A0113 = a21 * a42 - a22 * a41;
        float A0112 = a21 * a32 - a22 * a31;

        INS_INC1(add, 7);
        INS_INC1(mul, 16);

        float det = a11 * (a22 * A2323 - a23 * A1323 + a24 * A1223) - a12 * (a21 * A2323 - a23 * A0323 + a24 * A0223) +
                    a13 * (a21 * A1323 - a22 * A0323 + a24 * A0123) - a14 * (a21 * A1223 - a22 * A0223 + a23 * A0123);
        INS_INC1(div, 1);
        det = 1 / det;

        m44 inv;

        INS_INC1(add, 40);
        INS_INC1(mul, 64);

        inv.val[0][0] = det * (a22 * A2323 - a23 * A1323 + a24 * A1223);
        inv.val[0][1] = det * -(a12 * A2323 - a13 * A1323 + a14 * A1223);
        inv.val[0][2] = det * (a12 * A2313 - a13 * A1313 + a14 * A1213);
        inv.val[0][3] = det * -(a12 * A2312 - a13 * A1312 + a14 * A1212);
        inv.val[1][0] = det * -(a21 * A2323 - a23 * A0323 + a24 * A0223);
        inv.val[1][1] = det * (a11 * A2323 - a13 * A0323 + a14 * A0223);
        inv.val[1][2] = det * -(a11 * A2313 - a13 * A0313 + a14 * A0213);
        inv.val[1][3] = det * (a11 * A2312 - a13 * A0312 + a14 * A0212);
        inv.val[2][0] = det * (a21 * A1323 - a22 * A0323 + a24 * A0123);
        inv.val[2][1] = det * -(a11 * A1323 - a12 * A0323 + a14 * A0123);
        inv.val[2][2] = det * (a11 * A1313 - a12 * A0313 + a14 * A0113);
        inv.val[2][3] = det * -(a11 * A1312 - a12 * A0312 + a14 * A0112);
        inv.val[3][0] = det * -(a21 * A1223 - a22 * A0223 + a23 * A0123);
        inv.val[3][1] = det * (a11 * A1223 - a12 * A0223 + a13 * A0123);
        inv.val[3][2] = det * -(a11 * A1213 - a12 * A0213 + a13 * A0113);
        inv.val[3][3] = det * (a11 * A1212 - a12 * A0212 + a13 * A0112);
        return inv;
    }

    m44 get_rot_matrix(vec rot) {
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

        m44 m = identity;

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

    m33 get_rot_matrix_33(vec rot){
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

    m33 get_rot_from_trans(m44 trans){
        m33 m = identity_33;
        m.val[0][0] = trans.val[0][0];
        m.val[0][1] = trans.val[0][1];
        m.val[0][2] = trans.val[0][2];
        m.val[1][0] = trans.val[1][0];
        m.val[1][1] = trans.val[1][1];
        m.val[1][2] = trans.val[1][2];
        m.val[2][0] = trans.val[2][0];
        m.val[2][1] = trans.val[2][1];
        m.val[2][2] = trans.val[2][2];

        return m; 
    }

    /**
     * Computes only inverse for rotations, i.e computes the transpose. have to cleanup the code
     */
    m33 m33_inv(m33 m){
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

    /**
     * Creates a object-space to world-space transformation matrix
     *
     * pos is the origin of the object space and the direction of the axis is given
     * as rotation along the three axis in degrees.
     */
    m44 get_transf_matrix(vec pos, vec rot) {
        m44 camera_matrix = get_rot_matrix(rot);

        camera_matrix.val[0][3] = pos.x;
        camera_matrix.val[1][3] = pos.y;
        camera_matrix.val[2][3] = pos.z;

        return camera_matrix;
    }

} // namespace impl::opt3
