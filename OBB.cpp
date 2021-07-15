#include "OBB.h"
#include "CommonFunctions.hpp"

using namespace JC::Math;

OBB::Mtx44 OBB::Mtx44::CreateWithRotation(float degrees, float axisX, float axisY, float axisZ)
{
	Mtx44 newmat4;
	double mag = sqrt(axisX * axisX + axisY * axisY + axisZ * axisZ);
	if (Math::Fabs((float)mag) < EPSILON)
		throw "divide by zero";
	double x = axisX / mag, y = axisY / mag, z = axisZ / mag;
	double c = cos(degrees * PI / 180), s = sin(degrees * PI / 180);
	newmat4.a[0] = (float)(x * x * (1.f - c) + c);
	newmat4.a[1] = (float)(y * x * (1.f - c) + z * s);
	newmat4.a[2] = (float)(x * z * (1.f - c) - y * s);
	newmat4.a[3] = 0;
	newmat4.a[4] = (float)(x * y * (1.f - c) - z * s);
	newmat4.a[5] = (float)(y * y * (1.f - c) + c);
	newmat4.a[6] = (float)(y * z * (1.f - c) + x * s);
	newmat4.a[7] = 0;
	newmat4.a[8] = (float)(x * z * (1.f - c) + y * s);
	newmat4.a[9] = (float)(y * z * (1.f - c) - x * s);
	newmat4.a[10] = (float)(z * z * (1.f - c) + c);
	newmat4.a[11] = 0;
	newmat4.a[12] = 0;
	newmat4.a[13] = 0;
	newmat4.a[14] = 0;
	newmat4.a[15] = 1;
	return newmat4;
}

OBB::Mtx44 OBB::Mtx44::operator*(const Mtx44 & rhs) const
{
	Mtx44 ret;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			ret.a[i * 4 + j] = a[0 * 4 + j] * rhs.a[i * 4 + 0] + a[1 * 4 + j] * rhs.a[i * 4 + 1] + a[2 * 4 + j] * rhs.a[i * 4 + 2] + a[3 * 4 + j] * rhs.a[i * 4 + 3];
	return ret;
}

Vector3 OBB::Mtx44::operator*(const Vector3 & rhs) const
{
	float b[4];
	for (int i = 0; i < 4; i++)
		b[i] = a[0 * 4 + i] * rhs.x + a[1 * 4 + i] * rhs.y + a[2 * 4 + i] * rhs.z;
	return Vector3(b[0], b[1], b[2]);
}

OBB::Mtx44 OBB::Mtx44::getRotation(Vector3& degreeAxis, RotationOrder rotOrder) {
	Mtx44 result;
	Mtx44 x, y, z;
	x = Mtx44::CreateWithRotation(degreeAxis.x, 1, 0, 0);
	y = Mtx44::CreateWithRotation(degreeAxis.y, 0, 1, 0);
	z = Mtx44::CreateWithRotation(degreeAxis.z, 0, 0, 1);

	switch (rotOrder) {
	case XYZ:
		result = z * y * x;
	case ZYX:
		result = x * y * z;
	}

	return result;
}

OBB::OBB()
{
	isinit = false;
}

OBB::OBB(Vector3& obj1_pos, Vector3& obj1_rot, RotationOrder obj1_rotOrder,
	Vector3& obj1_scale, Vector3& obj1_maxpoints, Vector3& obj1_minpoints,
	Vector3& obj2_pos, Vector3& obj2_rot, RotationOrder obj2_rotOrder,
	Vector3& obj2_scale, Vector3& obj2_maxpoints, Vector3& obj2_minpoints,
	bool RadCheckEnabled, float obj1_rad, float obj2_rad):
	obj1_pos(obj1_pos),
	obj1_rot(obj1_rot),
	obj1_rotOrder(obj1_rotOrder),
	obj1_scale(obj1_scale),
	obj1_maxpoints(obj1_maxpoints),
	obj1_minpoints(obj1_minpoints),
	obj2_pos(obj2_pos),
	obj2_rot(obj2_rot),
	obj2_rotOrder(obj2_rotOrder),
	obj2_scale(obj2_scale),
	obj2_maxpoints(obj2_maxpoints),
	obj2_minpoints(obj2_minpoints),
	RadCheckEnabled(RadCheckEnabled),
	obj1_rad(obj1_rad),
	obj2_rad(obj2_rad),
	isinit(true)
{
}

OBB::~OBB()
{
}

void OBB::init(Vector3& obj1_pos, Vector3& obj1_rot, RotationOrder obj1_rotOrder,
	Vector3& obj1_scale, Vector3& obj1_maxpoints, Vector3& obj1_minpoints,
	Vector3& obj2_pos, Vector3& obj2_rot, RotationOrder obj2_rotOrder,
	Vector3& obj2_scale, Vector3& obj2_maxpoints, Vector3& obj2_minpoints,
	bool RadCheckEnabled, float obj1_rad, float obj2_rad) {
	this->obj1_pos = obj1_pos;
	this->obj1_rot = obj1_rot;
	this->obj1_rotOrder = obj1_rotOrder;
	this->obj1_scale = obj1_scale;
	this->obj1_maxpoints = obj1_maxpoints;
	this->obj1_minpoints = obj1_minpoints;
	this->obj2_pos = obj2_pos;
	this->obj2_rot = obj2_rot;
	this->obj2_rotOrder = obj2_rotOrder;
	this->obj2_scale = obj2_scale;
	this->obj2_maxpoints = obj2_maxpoints;
	this->obj2_minpoints = obj2_minpoints;
	this->RadCheckEnabled = RadCheckEnabled;
	this->obj1_rad = obj1_rad;
	this->obj2_rad = obj2_rad;
	isinit = true;
}

bool OBB::getResult()
{
	if (!isinit) {
		return false;
	}
	if (RadCheckEnabled) {
		if (!checkRad())
			return false;
	}
	Vector3 obj1_coords[8], obj2_coords[8];
	std::vector<Vector3> normals;
	normals.reserve(4);
	float obj1_prj[8], obj2_prj[8];
	float obj1_min, obj1_max, obj2_min, obj2_max;
	setcoords(obj1_coords, obj1_rot, obj1_rotOrder, obj1_minpoints, obj1_maxpoints, obj1_pos, obj1_scale);
	setcoords(obj2_coords, obj2_rot, obj2_rotOrder, obj2_minpoints, obj2_maxpoints, obj2_pos, obj2_scale);
	setnormals(normals, obj1_rot, obj1_rotOrder, obj1_maxpoints, obj1_pos, obj1_scale);
	setnormals(normals, obj2_rot, obj2_rotOrder, obj2_maxpoints, obj2_pos, obj2_scale);
	int num_normals = normals.size();
	for (int i = 0; i < num_normals; ++i) {
		getprojection(obj1_coords, normals[i], obj1_prj);
		getprojection(obj2_coords, normals[i], obj2_prj);
		getminmax(obj1_min, obj1_max, obj2_min, obj2_max, obj1_prj, obj2_prj);
		if (obj1_max < obj2_min || obj2_max < obj1_min) {
			return false;
		}
	}
	return true;
}

void OBB::setcoords(Vector3* coords, Vector3& rot, RotationOrder rotOrder, Vector3& minPoints, Vector3& maxPoints, Vector3& pos, Vector3& scale)
{
	Mtx44 rotation = Mtx44::getRotation(rot, rotOrder);
	coords[0].Set(maxPoints.x, minPoints.y, maxPoints.z);
	coords[1].Set(minPoints.x, minPoints.y, maxPoints.z);
	coords[2].Set(maxPoints.x, minPoints.y, minPoints.z);
	coords[3].Set(minPoints.x, minPoints.y, minPoints.z);
	coords[4].Set(maxPoints.x, maxPoints.y, maxPoints.z);
	coords[5].Set(minPoints.x, maxPoints.y, maxPoints.z);
	coords[6].Set(maxPoints.x, maxPoints.y, minPoints.z);
	coords[7].Set(minPoints.x, maxPoints.y, minPoints.z);
	for (int i = 0; i < 8; ++i) {
		coords[i].x *= scale.x;
		coords[i].y *= scale.y;
		coords[i].z *= scale.z;
		coords[i] = rotation * coords[i];
		coords[i] += pos;
	}
}

void OBB::setnormals(std::vector<Vector3>& normals, Vector3& rot, RotationOrder rotOrder, Vector3& maxPoints, Vector3& pos, Vector3& scale)
{
	Vector3 temp1, temp2, temp3;
	Mtx44 rotation = Mtx44::getRotation(rot, rotOrder);
	temp1.Set(maxPoints.x, 0, 0);
	temp2.Set(0, 0, maxPoints.z);
	temp3.Set(0, maxPoints.y , 0);
	temp1 = rotation * temp1;
	temp2 = rotation * temp2;
	temp3 = rotation * temp3;
	normals.push_back(temp1);
	normals.push_back(temp2);
	normals.push_back(temp3);
}

void OBB::getprojection(Vector3* coords, Vector3& normal, float* proj)
{
	for (int i = 0; i < 8; ++i) {
		proj[i] = coords[i].Dot(normal);
	}
}

void OBB::getminmax(float& minA, float& maxA, float& minB, float& maxB, float* prj1, float* prj2)
{
	minA = prj1[0]; maxA = prj1[0]; minB = prj2[0]; maxB = prj2[0];
	for (int i = 0; i < 8; ++i) {
		minA = Math::Min(minA, prj1[i]);
		maxA = Math::Max(maxA, prj1[i]);
		minB = Math::Min(minB, prj2[i]);
		maxB = Math::Max(maxB, prj2[i]);
	}
}

bool OBB::checkRad()
{
	Vector3 p1 = obj1_pos, p2 = obj2_pos;
	p1.y = p2.y = 0;
	float displacement = (p1 - p2).LengthSquared();
	float radiusSquared = (obj1_rad + obj1_rad) * (obj2_rad + obj2_rad);
	return displacement < radiusSquared ? true : false;
}
