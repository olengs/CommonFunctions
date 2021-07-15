#pragma once
#include "TemplatedFunc.hpp"
#include <iostream>
#include <vector>
namespace jc_m = JC::Math;

class OBB {

public:
	enum RotationOrder {
		XYZ,
		ZYX,
	};
	class Mtx44 {
	public:
		static Mtx44 CreateWithRotation(float degrees, float axisX, float axisY, float axisZ);
		static Mtx44 getRotation(jc_m::Vector3& degreeAxis, RotationOrder rotOrder);
		Mtx44 operator*(const Mtx44& rhs) const;
		jc_m::Vector3 operator*(const jc_m::Vector3& rhs) const;
		float a[16];
	};

	OBB();
	OBB(jc_m::Vector3& obj1_pos, jc_m::Vector3& obj1_rot, RotationOrder obj1_rotOrder,
		jc_m::Vector3& obj1_scale, jc_m::Vector3& obj1_maxpoints, jc_m::Vector3& obj1_minpoints,
		jc_m::Vector3& obj2_pos, jc_m::Vector3& obj2_rot, RotationOrder obj2_rotOrder,
		jc_m::Vector3& obj2_scale, jc_m::Vector3& obj2_maxpoints, jc_m::Vector3& obj2_minpoints,
		bool RadCheckEnabled = false, float obj1_rad = 0.f, float obj2_rad = 0.f);
	~OBB();
	void init(jc_m::Vector3& obj1_pos, jc_m::Vector3& obj1_rot, RotationOrder obj1_rotOrder,
		jc_m::Vector3& obj1_scale, jc_m::Vector3& obj1_maxpoints, jc_m::Vector3& obj1_minpoints,
		jc_m::Vector3& obj2_pos, jc_m::Vector3& obj2_rot, RotationOrder obj2_rotOrder,
		jc_m::Vector3& obj2_scale, jc_m::Vector3& obj2_maxpoints, jc_m::Vector3& obj2_minpoints,
		bool RadCheckEnabled = false, float obj1_rad = 0.f, float obj2_rad = 0.f);
	bool getResult();

private:


	static void setcoords(jc_m::Vector3* coords, jc_m::Vector3& rot, RotationOrder rotOrder, jc_m::Vector3& minPoints, jc_m::Vector3& maxPoints, jc_m::Vector3& pos, jc_m::Vector3& scale);
	static void setnormals(std::vector<jc_m::Vector3>& normals, jc_m::Vector3& rot, RotationOrder rotOrder, jc_m::Vector3& maxPoints, jc_m::Vector3& pos, jc_m::Vector3& scale);
	static void getprojection(jc_m::Vector3* coords, jc_m::Vector3& normal, float* proj);
	static void getminmax(float& minA, float& maxA, float& minB, float& maxB, float* prj1, float* prj2);
	bool checkRad();

	jc_m::Vector3 obj1_pos;
	jc_m::Vector3 obj1_rot;
	RotationOrder obj1_rotOrder;
	jc_m::Vector3 obj1_scale;
	jc_m::Vector3 obj1_maxpoints;
	jc_m::Vector3 obj1_minpoints;

	jc_m::Vector3 obj2_pos;
	jc_m::Vector3 obj2_rot;
	RotationOrder obj2_rotOrder;
	jc_m::Vector3 obj2_scale;
	jc_m::Vector3 obj2_maxpoints;
	jc_m::Vector3 obj2_minpoints;

	bool RadCheckEnabled;
	float obj1_rad;
	float obj2_rad;

	bool isinit;

};