#include "include/cudaTrace.h"
#include "include/vec3.h"
#include <stdlib.h>
#include <iostream>

constexpr float eps = 0.000005;
constexpr float aspectRatio = 512.0 / 384.0;

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif

// -----------------------------------------------------------
// Data structures
// SoA used for most entities, except for the triangles because
// they are retrieved scattered through memory and need all 3
// vertices everytime anyways.
// -----------------------------------------------------------

struct rays
{
	vec3* O_;
	vec3* D_;
	float* t_;
	int* primIdx_;
	int* pixelId_;
	int* objId_;
};

struct materials {
	vec3* albedo;
	float* reflectivity;//0 = diffuse, 1= mirror
	bool* isLight;//if its light the albedo will be the emittance
	bool* isDialectric;//if its a dialectric the index will be used for the fresnel equation
	float* index;
};

struct triangleCount {
	int index_;
	int count_;
};

struct triangle {
	vec3 V0_;
	vec3 V1_;
	vec3 V2_;
};

struct object {
	BVHutil* BVHpool;
	triangle* Triangles;
	vec3* normals;
	int* matIndex;
};

std::ostream& operator<<(std::ostream& os, const triangle& dt) {
	for (int i = 0; i < 3; ++i)
		os << dt.V0_[i] << " ";
	os << "\n";
	for (int i = 0; i < 3; ++i)
		os << dt.V1_[i] << " ";
	os << "\n";
	for (int i = 0; i < 3; ++i)
		os << dt.V2_[i] << " ";
	os << "\n";
	return os;
}

std::ostream& operator<<(std::ostream& os, const BVHutil& dt) {
	os << "--------------BVHUTIL--------------\n";
	for (int i = 0; i < 3; ++i)
		os << dt.minBounds_[i] << " ";
	os << "\n";
	for (int i = 0; i < 3; ++i)
		os << dt.maxBounds_[i] << " ";
	os << "\n";
	os << dt.isLeaf_ << " ";
	os << dt.left_ << " " << dt.right_ << "\n";
	os << dt.first_ << " " << dt.count_ << "\n";
	os << "----------------------------------\n";
	return os;
}

int nx;
int ny;
int nOfRays;

// -----------------------------------------------------------
// Buffers
// -----------------------------------------------------------
constexpr int numberOfObjects = 2;// need to change this if more objects are added
object* h_objects = new object [numberOfObjects];
object* d_objects;

materials h_mats;
materials d_mats;

rays d_rays;
rays d_newrays;

materials d_materials;

vec3* h_fb;
vec3* d_fb;

//for rng
unsigned int* h_seeds;
unsigned int* d_seeds;

vec3* h_T;
vec3* d_T;
vec3* h_E;
vec3* d_E;

unsigned char* d_skydomeData;
__device__ int d_skydomeNx; 
__device__ int d_skydomeNy;
__device__ vec3 d_skydomeCenter;//0,0,0
__device__ float d_skydomeR2;//10000000

// -----------------------------------------------------------
// Bunch of device functions required for path tracing
// -----------------------------------------------------------

__device__ vec3 intersectionPoint(vec3& O, vec3& D, vec3& N, float t) {
	return ((O + t * D) + N * eps);
}

// RNG - Marsaglia's xor32
__host__ __device__ unsigned int& aRandomUInt(unsigned int& seed)
{
	seed ^= seed << 13;
	seed ^= seed >> 17;
	seed ^= seed << 5;
	return seed;
}
__host__ __device__ float aRandomFloat(unsigned int& seed) { return aRandomUInt(seed) * 2.3283064365387e-10f; }
__host__ __device__ float aRand(unsigned int& seed, float range) { return aRandomFloat(seed) * range; }

__host__ __device__ unsigned int wang_hash(unsigned int seed) {
	seed = (seed ^ 61) ^ (seed >> 16);
	seed *= 9;
	seed = seed ^ (seed >> 4);
	seed *= 0x27d4eb2d;
	seed = seed ^ (seed >> 15);
	return seed;
}

__device__ vec3 getRandomDirection(vec3& N, float r1, float r2) {
	//setup a tangent space with normal vector N, tangent vector T, and bitangent vector B
	vec3 W;
	if (fabs(N[0]) > 0.99)
		W = vec3(0, 1, 0);
	else
		W = vec3(1, 0, 0);
	vec3 T = cross(N, W);
	T.make_unit_vector();
	vec3 B = cross(T, N);
	B.make_unit_vector();

	//x,y, and z are random points on the unit hemisphere
	float x = cos(2 * M_PI * r1) * sqrtf(1.0f - r2 * r2);
	float y = sin(2 * M_PI * r1) * sqrtf(1.0f - r2 * r2);
	float z = r2;
	vec3 P(x, y, z);
	P.make_unit_vector();

	//transform the vector to tangent space
	float Px = P[0] * T[0] + P[1] * B[0] + P[2] * N[0];
	float Py = P[0] * T[1] + P[1] * B[1] + P[2] * N[1];
	float Pz = P[0] * T[2] + P[1] * B[2] + P[2] * N[2];

	vec3 Protated(Px, Py, Pz);
	Protated.make_unit_vector();

	return Protated;
}

__device__ bool rayBoxIntersect(vec3& O, vec3& D, float mins[3], float maxs[3]) {
	float tmin = -INFINITY, tmax = INFINITY;

	float tx1 = (mins[0] - O.x()) / D.x();
	float tx2 = (maxs[0] - O.x()) / D.x();

	tmin = max(tmin, min(tx1, tx2));
	tmax = min(tmax, max(tx1, tx2));

	float ty1 = (mins[1] - O.y()) / D.y();
	float ty2 = (maxs[1] - O.y()) / D.y();

	tmin = max(tmin, min(ty1, ty2));
	tmax = min(tmax, max(ty1, ty2));

	float tz1 = (mins[2] - O.z()) / D.z();
	float tz2 = (maxs[2] - O.z()) / D.z();

	tmin = max(tmin, min(tz1, tz2));
	tmax = min(tmax, max(tz1, tz2));

	return (tmax >= tmin && tmax > 0);
}

__device__ float rayTriangleIntersect(vec3& O, vec3& D, vec3& V0, vec3& V1, vec3& V2) {
	vec3 V0V1 = V1 - V0;
	vec3 V0V2 = V2 - V0;
	vec3 Pvec = cross(D, V0V2);
	float det = dot(V0V1, Pvec);

	if (fabs(det) < eps) return INFINITY;
	float invDet = 1 / det;

	vec3 Tvec = O - V0;
	float u = dot(Tvec, Pvec) * invDet;
	if (u < 0 || u > 1) return INFINITY;

	vec3 Qvec = cross(Tvec, V0V1);
	float v = dot(D, Qvec) * invDet;
	if (v < 0 || u + v > 1) return INFINITY;

	float t = dot(V0V2, Qvec) * invDet;

	return t;
}

__device__ float raySphereIntersect(vec3& O, vec3& D, vec3& center, float radius2) {
	vec3 C = center - O;
	float t = dot(C, D);
	vec3 Q = C - t * D;
	float p2 = dot(Q, Q);
	if (p2 > radius2)
		return INFINITY; // r2 = r * r
	t -= sqrt(radius2 - p2);
	return t;
}

__device__ void get_sphere_uv(const vec3& tp, float& u, float& v) {
	float phi = atan2(tp[2], tp[0]);
	float theta = asin(tp[1]);
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / M_PI;
}

__device__ vec3 getSkydomeColor(float u, float v, unsigned char* skydome) {
	if (u < 0)		u = 0;
	if(v < 0)		v = 0;
	if (u > 1) u = 1;
	if (v > 1) v = 1;

	int i = ((u)*(float)d_skydomeNx);
	int j = (1.0f - v)* (float)d_skydomeNy - 0.001;

	if (i < 0) i = 0;
	if (j < 0) j = 0;

	if (i > d_skydomeNx - 1) i = d_skydomeNx - 1;
	if (j > d_skydomeNy - 1) j = d_skydomeNy - 1;

	float r = int(skydome[3 * i + 3 * d_skydomeNx * j]) / 255.0;
	float g = int(skydome[3 * i + 3 * d_skydomeNx * j + 1]) / 255.0;
	float b = int(skydome[3 * i + 3 * d_skydomeNx * j + 2]) / 255.0;

	return { r,g,b };
}

__device__ vec3 reflect(vec3& i, vec3& n) { return i - 2.0f * n * dot(n, i); }

// -----------------------------------------------------------
// 1) Camera ray generations 
// 2) Ray extension (find intersections)
// 3) Shading and new ray generations
// no seperate shade and connect kernel because not doing NEE, and
// also because not using a complicated material class
// -----------------------------------------------------------
__global__ void generateCameraRays(rays r, int max_x, int max_y, vec3 P0, vec3 P1, vec3 P2, vec3 E) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	if ((i >= max_x) || (j >= max_y)) return;
	int pixel_index = j * max_x + i;

	float u = float(i) / float(max_y);
	float v = float(j) / float(max_x) * aspectRatio;

	r.O_[pixel_index] = E;

	vec3 Puv = P0 + u * (P1 - P0) + v * (P2 - P0);
	Puv -= E;
	Puv /= Puv.length();

	r.D_[pixel_index] = Puv;
	r.t_[pixel_index] = INFINITY;
	r.pixelId_[pixel_index] = pixel_index;
}

__global__ void extendRays(rays r, object* objects, int nObjects) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= nObjects)
		return;

	vec3 O = r.O_[i];
	vec3 D = r.D_[i];
	float dist = INFINITY;
	int primIndex = -1;
	int bvhObjectIndex = -1;

	for (int objIndex = 0;objIndex < numberOfObjects; ++objIndex) {
		BVHutil currentNode = objects[objIndex].BVHpool[0];
		float* mins = currentNode.minBounds_;
		float* maxs = currentNode.maxBounds_;

		if (rayBoxIntersect(O, D, mins, maxs)) {
			int stack[64];
			triangleCount t_stack[64];
			int t_stackIndex = -1;

			int stackIndex = 1;
			stack[0] = currentNode.left_;
			stack[1] = currentNode.right_;

			while (stackIndex >= 0 || t_stackIndex >= 0) {
				if (stackIndex >= 0) {
					currentNode = objects[objIndex].BVHpool[stack[stackIndex--]];//pop last node, decrease stackIndex
					BVHutil leftNode = objects[objIndex].BVHpool[currentNode.left_];
					mins = leftNode.minBounds_;
					maxs = leftNode.maxBounds_;
					if (rayBoxIntersect(O, D, mins, maxs)) {
						if (leftNode.isLeaf_) {
							triangleCount ttemp;
							ttemp.index_ = leftNode.first_;
							ttemp.count_ = leftNode.count_;
							t_stack[++t_stackIndex] = ttemp;
						}
						else
							stack[++stackIndex] = currentNode.left_;
					}

					BVHutil rightNode = objects[objIndex].BVHpool[currentNode.right_];
					mins = rightNode.minBounds_;
					maxs = rightNode.maxBounds_;
					if (rayBoxIntersect(O, D, mins, maxs)) {
						if (rightNode.isLeaf_) {
							triangleCount ttemp;
							ttemp.index_ = rightNode.first_;
							ttemp.count_ = rightNode.count_;
							t_stack[++t_stackIndex] = ttemp;
						}
						else
							stack[++stackIndex] = currentNode.right_;
					}
				}

				if (t_stackIndex >= 0) {
				
					triangleCount t_pair = t_stack[t_stackIndex];
					triangle triangle = objects[objIndex].Triangles[t_pair.index_ + t_pair.count_ - 1];

					float t = rayTriangleIntersect(O, D, triangle.V0_, triangle.V1_, triangle.V2_);
					if (t < dist && t >= 0) {
						dist = t;
						primIndex = t_pair.index_ + t_pair.count_ - 1;
						bvhObjectIndex = objIndex;
					}

					if (t_pair.count_ > 1) {
						--t_stack[t_stackIndex].count_;
					}
					else --t_stackIndex;
				}
			}
		}
	}

	r.objId_[i] = bvhObjectIndex;
	float lightDist = raySphereIntersect(O, D, { 0,-100,-100 }, 10000);
	if (lightDist < dist && lightDist >= 0) {//it hit the light
		r.t_[i] = lightDist;
		r.primIdx_[i] = -1;
	}
	else {
		r.t_[i] = dist;
		r.primIdx_[i] = primIndex;
	}
}

__device__ unsigned int extensionRayIdx;//extension ray buffer index

__global__ void shadeRays(rays rayBuffer, rays newRayBuffer, object* objects, vec3* T, unsigned int* seeds, int nObjects, vec3* fb, materials mats, unsigned char* skydome) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= nObjects)
		return;

	float t = rayBuffer.t_[i];
	int pixelIndex = rayBuffer.pixelId_[i];

	if (t == INFINITY) {//hit blackness
		float u, v;
		vec3 D = rayBuffer.D_[i];
		get_sphere_uv(D, u, v);
		vec3 color = getSkydomeColor(u, v, skydome);
		vec3 E = T[pixelIndex] * color;// vec3{ 1,1,1 };// color;
		fb[pixelIndex] = E;
		return;
	}

	//primIndex of -1 means the ray hit a light
	int primIndex = rayBuffer.primIdx_[i];
	
	if (primIndex < 0) {//hit a lightsource
		vec3 emittance = { 2,2,2 };
		vec3 E = T[pixelIndex] * emittance;

		fb[pixelIndex] = E;
		return;
	}

	vec3 O = rayBuffer.O_[i];
	vec3 D = rayBuffer.D_[i];

	//objID is the id of the object in which the closest intersection was found
	int objID = rayBuffer.objId_[i];
	vec3 N = objects[objID].normals[primIndex];

	int matIndex = objects[objID].matIndex[primIndex];
	float reflectivity = mats.reflectivity[matIndex];

	unsigned int seed = seeds[i];
	float r = aRandomFloat(seed);//use this to determine if it will reflect or not
	vec3 newDirection;

	if (r < reflectivity) {
		//reflective bounce
		newDirection = reflect(D, N);
	}
	else {
		//diffuse bounce
		float r1 = aRandomFloat(seed);
		float r2 = aRandomFloat(seed);
		newDirection = getRandomDirection(N, r1, r2);
	}
	
	vec3 I = intersectionPoint(O, D, N, t);
	I += N * eps;

	newDirection.make_unit_vector();

	int ei = atomicAdd(&extensionRayIdx, 1);
	newRayBuffer.D_[ei] = newDirection;
	newRayBuffer.O_[ei] = I;
	newRayBuffer.pixelId_[ei] = pixelIndex;

	//in case of a diffuse bounce update T, the brdf is a fixed albedo (color) currently.
	vec3 BRDF = mats.albedo[matIndex];
	//vec3 BRDF(0.8, 0.5, 0.7);
	BRDF /= M_PI;
	float NdotR = dot(N, newDirection);
	if (NdotR < 0)
		NdotR = 0;
	vec3 newt = M_PI * 2.0f * BRDF * NdotR;
	T[pixelIndex] *= newt;
	seeds[i] = seed;
}

// -----------------------------------------------------------
// 1) Camera ray generations 
// 2) tracing of the rays
// 3) retrieve the frame buffer
// -----------------------------------------------------------

void generateRays(float P0t[3], float P1t[3], float P2t[3], float Et[3]) {
	//reset the T to {1,1,1} and the framebuffer to {0,0,0} 
	cudaMemcpy(d_T, h_T, sizeof(vec3) * nOfRays, cudaMemcpyHostToDevice);
	cudaMemcpy(d_fb, h_E, sizeof(vec3) * nOfRays, cudaMemcpyHostToDevice);

	vec3 P0(P0t);
	vec3 P1(P1t);
	vec3 P2(P2t);
	vec3 E(Et);

	const int tx = 8;
	const int ty = 8;

	dim3 blocks(nx / tx + 1, ny / ty + 1);
	dim3 threads(tx, ty);

	generateCameraRays << <blocks, threads >> > (d_rays, nx, ny, P0, P1, P2, E);
}

void traceRays(int N) {
	int i = 0;
	int currentNOfRays = nOfRays;

	//currently fixed to 20 iterations
	while (i < 15) {
		int threadsPerBlock = 256;
		int blocksPerGrid = (currentNOfRays + threadsPerBlock - 1) / threadsPerBlock;
		//extending the camera ray
		extendRays << <blocksPerGrid, threadsPerBlock >> > (d_rays,d_objects, currentNOfRays);

		unsigned int zero = 0;
		cudaMemcpyToSymbol(extensionRayIdx, &zero, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);

		shadeRays << <blocksPerGrid, threadsPerBlock >> > (d_rays, d_newrays, d_objects, d_T, d_seeds, currentNOfRays, d_fb, d_mats,d_skydomeData);
		unsigned int nOfNewRays = 0;
		cudaMemcpyFromSymbol(&nOfNewRays, extensionRayIdx, sizeof(unsigned int), 0, cudaMemcpyDeviceToHost);

		//switcheroo
		rays tempRays = d_newrays;
		d_rays = d_newrays;
		d_newrays = tempRays;

		currentNOfRays = nOfNewRays;
		++i;
	}
}

float* getFrameBuffer() {
	size_t size = sizeof(vec3) * nx * ny;
	cudaMemcpy(h_fb, d_fb, size, cudaMemcpyDeviceToHost);

	float* fb = (float*)malloc(sizeof(float) * 3 * nx * ny);
	for (int i = 0, j = 0; i < nx * ny; ++i, j += 3) {
		fb[j] = h_fb[i][0];
		fb[j + 1] = h_fb[i][1];
		fb[j + 2] = h_fb[i][2];
	}

	return fb;
}

// -----------------------------------------------------------
// Intialization of all data required on the gpu
// -----------------------------------------------------------

void initializeRayBuffer(const int x, const int y) {
	nx = x; ny = y;
	const int nPixels = x * y;
	nOfRays = nPixels;

	//initialize the seeds for each pixel (which will be used by all the rays).
	//a pixel or primary ray does not need to use the same seed every time can use whatever is convenient
	h_seeds = (unsigned int*)malloc(sizeof(unsigned int) * nOfRays);
	h_T = (vec3*)malloc(sizeof(vec3) * nOfRays);
	h_E = (vec3*)malloc(sizeof(vec3) * nOfRays);
	unsigned int seed = 0x12345678;

	for (int i = 0; i < nOfRays; ++i) {
		h_seeds[i] = wang_hash(i);//change this with wang hash
		h_T[i] = { 1,1,1 };
		h_E[i] = { 0,0,0 };
	}

	cudaMalloc(&d_seeds, sizeof(unsigned int) * nOfRays);
	cudaMemcpy(d_seeds, h_seeds, sizeof(unsigned int) * nOfRays, cudaMemcpyHostToDevice);

	cudaMalloc(&d_T, sizeof(vec3) * nOfRays);

	cudaMalloc(&d_rays.D_, sizeof(vec3) * nPixels);
	cudaMalloc(&d_rays.O_, sizeof(vec3) * nPixels);
	cudaMalloc(&d_rays.t_, sizeof(float) * nPixels);
	cudaMalloc(&d_rays.primIdx_, sizeof(int) * nPixels);
	cudaMalloc(&d_rays.pixelId_, sizeof(int) * nPixels);
	cudaMalloc(&d_rays.objId_, sizeof(int) * nPixels);

	cudaMalloc(&d_newrays.D_, sizeof(vec3) * nPixels);
	cudaMalloc(&d_newrays.O_, sizeof(vec3) * nPixels);
	cudaMalloc(&d_newrays.t_, sizeof(float) * nPixels);
	cudaMalloc(&d_newrays.primIdx_, sizeof(int) * nPixels);
	cudaMalloc(&d_newrays.pixelId_, sizeof(int) * nPixels);
	cudaMalloc(&d_newrays.objId_, sizeof(int) * nPixels);

	h_fb = (vec3*)malloc(sizeof(vec3) * nPixels);
	cudaMalloc(&d_fb, sizeof(vec3) * nPixels);
}

void addMat(float* albedo, float reflectivity, bool isLight, bool isDialectric, float index, int matIndex) {
	h_mats.albedo[matIndex][0] = albedo[0];
	h_mats.albedo[matIndex][1] = albedo[1];
	h_mats.albedo[matIndex][2] = albedo[2];
	h_mats.reflectivity[matIndex] = reflectivity;
	h_mats.isLight[matIndex] = isLight;
	h_mats.isDialectric[matIndex] = isDialectric;
	h_mats.index[matIndex] = index;
}

void initializeObjs() {
	cudaMalloc((void**)&d_objects, sizeof(*d_objects) * numberOfObjects);

	h_mats.albedo = new vec3[numberOfObjects];
	h_mats.reflectivity = new float[numberOfObjects];
	h_mats.isDialectric = new bool[numberOfObjects];
	h_mats.isLight = new bool[numberOfObjects];
	h_mats.index = new float[numberOfObjects];
}

void allocateMats() {
	cudaMalloc(&d_mats.albedo, sizeof(vec3) * numberOfObjects);
	cudaMalloc(&d_mats.reflectivity, sizeof(float) * numberOfObjects);
	cudaMalloc(&d_mats.isDialectric, sizeof(bool) * numberOfObjects);
	cudaMalloc(&d_mats.isLight, sizeof(bool) * numberOfObjects);
	cudaMalloc(&d_mats.index, sizeof(float) * numberOfObjects);

	cudaMemcpy(d_mats.albedo, h_mats.albedo, sizeof(vec3) * numberOfObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mats.reflectivity, h_mats.reflectivity, sizeof(float) * numberOfObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mats.isDialectric, h_mats.isDialectric, sizeof(bool) * numberOfObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mats.isLight, h_mats.isLight, sizeof(bool) * numberOfObjects, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mats.index, h_mats.index, sizeof(float) * numberOfObjects, cudaMemcpyHostToDevice);
}

int bvhIndex = 0;

void initializeBVH(BVHutil* BVHpool, int N) {
	BVHutil* d_BVHpool;
	cudaMalloc((void**)&d_BVHpool, N * sizeof(BVHutil));

	h_objects[bvhIndex].BVHpool = BVHpool;

	cudaMemcpy(d_BVHpool, BVHpool, N * sizeof(BVHutil), cudaMemcpyHostToDevice);
	cudaMemcpy(&(d_objects[bvhIndex].BVHpool), &d_BVHpool, sizeof(d_objects[bvhIndex].BVHpool), cudaMemcpyHostToDevice);
}

void intializeSkydome(unsigned char* data, int nx, int ny, float* center, float r2) {
	cudaMemcpyToSymbol(d_skydomeNx, &nx, sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_skydomeNy, &ny, sizeof(int), 0, cudaMemcpyHostToDevice);
	vec3 cent{ center[0],center[1],center[2] };

	cudaMemcpyToSymbol(d_skydomeCenter, &cent, sizeof(vec3), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_skydomeR2, &r2, sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaMalloc(&d_skydomeData, sizeof(unsigned char*) * nx * ny * 3);
	cudaMemcpy(d_skydomeData, data, sizeof(unsigned char) * nx * ny * 3, cudaMemcpyHostToDevice);
}

void initializeTriangles(trianglec* triangles, int N, int matIndex) {
	size_t size = N * sizeof(triangle);

	h_objects[bvhIndex].Triangles = (triangle*)malloc(size);
	h_objects[bvhIndex].normals = (vec3*)malloc(sizeof(vec3) * N);
	h_objects[bvhIndex].matIndex = (int*)malloc(sizeof(int) * N);
	
	for (int i = 0; i < N; ++i) {
		triangle temp;
		temp.V0_ = triangles[i].V0_;
		temp.V1_ = triangles[i].V1_;
		temp.V2_ = triangles[i].V2_;
		h_objects[bvhIndex].normals[i] = triangles[i].N_;
		h_objects[bvhIndex].Triangles[i] = temp;
		h_objects[bvhIndex].matIndex[i] = matIndex;
	}

	triangle* d_triangles;
	vec3* d_normals;
	int* d_matIndex;

	cudaMalloc((void**)&d_normals, sizeof(vec3) * N);
	cudaMalloc((void**)&d_matIndex, sizeof(int) * N);
	cudaMalloc((void**)&d_triangles, size);

	cudaMemcpy(d_normals, h_objects[bvhIndex].normals, sizeof(vec3) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(d_matIndex, h_objects[bvhIndex].matIndex, sizeof(int) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(d_triangles, h_objects[bvhIndex].Triangles, size, cudaMemcpyHostToDevice);

	cudaMemcpy(&(d_objects[bvhIndex].normals), &d_normals, sizeof(d_objects[bvhIndex].normals), cudaMemcpyHostToDevice);
	cudaMemcpy(&(d_objects[bvhIndex].matIndex), &d_matIndex, sizeof(d_objects[bvhIndex].matIndex), cudaMemcpyHostToDevice);
	cudaMemcpy(&(d_objects[bvhIndex].Triangles), &d_triangles, sizeof(d_objects[bvhIndex].Triangles), cudaMemcpyHostToDevice);

	bvhIndex++;
}