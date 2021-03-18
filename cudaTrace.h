#pragma once
#include "util.h"

struct trianglec
{
	float V0_[3];
	float V1_[3];
	float V2_[3];
	float N_[3];
};

//generates the camera rays
void generateRays(float P0t[3], float P1t[3], float P2t[3], float Et[3]);
//extends and shades all the rays in a loop (20 iterations)
void traceRays(int N);
//read the framebuffer from the gpu
float* getFrameBuffer();

//functions to get all the data on the gpu
void initializeObjs();
void addMat(float* albedo, float reflectivity, bool isLight, bool isDialectric, float index, int matIndex);
void allocateMats();
void intializeSkydome(unsigned char* data, int nx, int ny, float* center, float r2);
void initializeObjs();
void initializeRayBuffer(const int nx, const int ny);
void initializeBVH(BVHutil* BVHpool, int N);
void initializeTriangles(trianglec* triangles, int N, int matIndex);


