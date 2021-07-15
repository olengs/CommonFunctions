#pragma once
#include <iostream>
#include <fstream>
#include "CommonFunctions.hpp"

struct OBJDimensions {
	float minX, minY, minZ, maxX, maxY, maxZ;
	float radius;

	OBJDimensions(float minX = 0.f, float minY = 0.f, float minZ = 0.f, float maxX = 0.f, float maxY = 0.f, float maxZ = 0.f, float radius = 0.f)
	{
		radius = 0.f;
		Set(minX, minY, minZ, maxX, maxY, maxZ, radius);
	}
	void Set(float minX, float minY, float minZ, float maxX, float maxY, float maxZ, float radius)
	{
		this->minX = minX; this->minY = minY; this->minZ = minZ;
		this->maxX = maxX; this->maxY = maxY; this->maxZ = maxZ;
		this->radius = radius;
	}

	void Loadcoord(const char* file_path)
	{
		std::ifstream fileStream(file_path, std::ios::binary);
		if (!fileStream.is_open())
		{
			std::cout << "Unable to open " << file_path << ". Are you in the right directory ?\n";
			return;
		}
		while (!fileStream.eof())
		{
			char buf[256];
			fileStream.getline(buf, 256);
			if (strncmp("v ", buf, 2) == 0)
			{
				JC::Math::Vector3 tempcoord;
				sscanf_s((buf + 2), "%f%f%f", &tempcoord.x, &tempcoord.y, &tempcoord.z);
				Set(JC::Math::Math::Min(tempcoord.x, minX), JC::Math::Math::Min(tempcoord.y, minY), JC::Math::Math::Min(tempcoord.z, minZ),
					JC::Math::Math::Max(tempcoord.x, maxX), JC::Math::Math::Max(tempcoord.y, maxY), JC::Math::Math::Max(tempcoord.z, maxZ), JC::Math::Math::Max(radius, tempcoord.Length()));
			}
		}
	}

	friend std::ostream& operator<<(std::ostream& os, OBJDimensions& d) {
		std::cout << "Dimensions: min: [" << d.minX << "," << d.minY << "," << d.minZ << "], max: [" << d.maxX << "," << d.maxY << "," << d.maxZ << "] radius : " << d.radius <<  " \n";
		return os;
	}
	
};