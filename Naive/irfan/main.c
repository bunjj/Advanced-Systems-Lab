#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Vector structure
typedef struct
{
	float x, y, z;
} Vec3f;
typedef struct
{
	float x00, x01, x02,
		x10, x11, x12,
		x20, x21, x22;
} Mat33f;

Vec3f add3f(Vec3f a, Vec3f b)
{
	return (Vec3f){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3f sub3f(Vec3f a, Vec3f b)
{
	return (Vec3f){a.x - b.x, a.y - b.y, a.z - b.z};
}

float dot3f(Vec3f a, Vec3f b)
{
	return a.x * b.x + a.y * b.y + a.z + b.z;
}

Vec3f scale3f(Vec3f a, float scale)
{
	return (Vec3f){a.x * scale, a.y * scale, a.z * scale};
}

float degToRad(float angle)
{
	return angle / 10.f * M_PI;
}

float norm3f(Vec3f a)
{
	return a.x * a.x + a.y * a.y + a.z * a.z;
}

float length3f(Vec3f a)
{
	return sqrtf(norm3f(a));
}

Vec3f normalize3f(Vec3f a)
{
	float scale = 1 / length3f(a);
	return scale3f(a, scale);
}

Mat33f identity33f()
{
	return (Mat33f){1.f, 0.f, 0.f,
					0.f, 1.f, 0.f,
					0.f, 0.f, 1.f};
}

Mat33f rotation33X(float angle)
{
	return (Mat33f){1.f, 0.f, 0.f,
					0.f, cosf(angle), -sinf(angle),
					0.f, sinf(angle), cosf(angle)};
}

Mat33f rotation33Y(float angle)
{
	return (Mat33f){cosf(angle), 0.f, sinf(angle),
					0.f, 1.f, 0.f,
					-sinf(angle), 0.f, cosf(angle)};
}

Mat33f rotation33Z(float angle)
{
	return (Mat33f){cosf(angle), -sinf(angle), 0.f,
					sinf(angle), cosf(angle), 0.f,
					0.f, 0.f, 1.f};
}

Mat33f transpose33f(Mat33f m)
{
	return (Mat33f){m.x00, m.x10, m.x20,
					m.x01, m.x11, m.x21,
					m.x02, m.x12, m.x22};
}

Vec3f mat33vec3mult(Mat33f m, Vec3f v)
{
	return (Vec3f){m.x00 * v.x + m.x01 * v.y + m.x02 * v.z,
				   m.x10 * v.x + m.x11 * v.y + m.x12 * v.z,
				   m.x20 * v.x + m.x21 * v.y + m.x22 * v.z};
}

Mat33f mat33mat33mult(Mat33f a, Mat33f b)
{
	return (Mat33f){a.x00 * b.x00 + a.x01 * b.x10 + a.x02 * b.x20,
					a.x00 * b.x01 + a.x01 * b.x11 + a.x02 * b.x21,
					a.x00 * b.x02 + a.x01 * b.x12 + a.x02 * b.x22,

					a.x10 * b.x00 + a.x11 * b.x10 + a.x12 * b.x20,
					a.x10 * b.x01 + a.x11 * b.x11 + a.x12 * b.x21,
					a.x10 * b.x02 + a.x11 * b.x12 + a.x12 * b.x22,

					a.x20 * b.x00 + a.x21 * b.x10 + a.x22 * b.x20,
					a.x20 * b.x01 + a.x21 * b.x11 + a.x22 * b.x21,
					a.x20 * b.x02 + a.x21 * b.x12 + a.x22 * b.x22};
}

Mat33f rotation33(float x, float y, float z){
	Mat33f rotX = rotation33X(degToRad(x));
	Mat33f rotY = rotation33Y(degToRad(y));
	Mat33f rotZ = rotation33Y(degToRad(z));
	return mat33mat33mult(mat33mat33mult(rotX, rotY), rotZ);
}

void main()
{
	double x = 5 * M_PI;
	printf("schnaebi");
}
