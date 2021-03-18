#pragma once

namespace Tmpl8 {

class Game
{
public:
	void SetTarget( Surface* surface ) { screen = surface; }
	void Init();
	void PresentFloatBuffer(float nFrames);
	void Shutdown();
	void Tick( float deltaTime );
	void MouseUp( int button ) { /* implement if you want to detect mouse button presses */ }
	void MouseDown( int button ) { /* implement if you want to detect mouse button presses */ }
	void MouseMove( int x, int y ) { /* implement if you want to detect mouse movement */ }
	void KeyUp( int key ) { /* implement if you want to handle keys */ }
	void KeyDown(int key);// { /* implement if you want to handle keys */ }
	Surface* screen;
	float3* floatBuffer;
};

}; // namespace Tmpl8

constexpr float eps = 0.0005f;//a small number

void get_sphere_uv(const float3& p, const float3 c, float& u, float& v);
bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1);

//Textures
class Texture {
public:
	virtual float3 getColor(const float u, const float v, const float3& p) = 0;
};

class ConstantTexture : public Texture {
public:
	ConstantTexture(const float3& color) : color_(color) {};

	virtual float3 getColor(const float u, const float v, const float3& p) {
		return color_;
	}
private:
	float3 color_;
};

class CheckerTexture : public Texture {
public:
	CheckerTexture(const float3& color1, const float3& color2) : color1_(color1), color2_(color2) {}

	virtual float3 getColor(float u, float v, const float3& p) {
		float sines = sin(4 * p.x) * sin(4 * p.y) * sin(4 * p.z);
		if (sines < 0)
			return color1_;
		else
			return color2_;
	}

private:
	float3 color1_;
	float3 color2_;
};

class ImageTexture : public Texture {
public:
	ImageTexture(const std::string& filename); 

	float3 getColor(const float u, const float v, const float3& p) {
		int i = (u)* nx;
		int j = (1 - v) * ny - 0.001;

		if (i < 0) i = 0;
		if (j < 0) j = 0;

		if (i > nx - 1) i = nx - 1;
		if (j > ny - 1) j = ny - 1;

		float r = int(data[3 * i + 3 * nx * j]) / 255.0;
		float g = int(data[3 * i + 3 * nx * j + 1]) / 255.0;
		float b = int(data[3 * i + 3 * nx * j + 2]) / 255.0;

		return { r,g,b };
	}

	~ImageTexture() { delete[] data; }

	unsigned char* data;
	int nx, ny;
private:
	
};

//Materials
enum class MaterialType { diffuse, reflective, dialectric, light };

class Material {
public:
	Texture* texture_;
	float reflectivity_;//in case of reflective
	float index_;//in case of dialectric
	float3 emittance_;//in case of light source
	MaterialType matType_;

	Material(const float3& color, const MaterialType type, const float reflectivity, const float index, const float3& emittance = { 1,1,1 }) :
		texture_(new ConstantTexture(color)), matType_(type), reflectivity_(reflectivity), index_(index), emittance_(emittance) {
	}

	Material(Texture* texture, const MaterialType type, const float reflectivity, const float index) :
		texture_(texture), matType_(type), reflectivity_(reflectivity), index_(index) {
	}

	float3 getColor(const float u, const float v, const float3& p) {
		return texture_->getColor(u, v, p);
	}
};

class Ray {
public:
	float3 O_;
	float3 D_;

	float t_;//closest intersection distance
	float3 N_;//normal of intersected surface
	bool inside;

	float u;
	float v;

	Material* mat_;

	Ray(float3& O, float3& D) : O_(O), D_(D) {
		t_ = std::numeric_limits<float>::max();
	}
};

//Primitives
class Primitive {
public:
	Material* mat_;

	Primitive(Material* mat) : mat_(mat) {}

	virtual bool intersect(Ray& ray) const = 0;
};

class Torus : public Primitive {
public:
	float R_;//inner radius
	float r_;//tube radius
	float3 C_;//center location

	Torus(const float R, const float r, const float3 C, Material* mat) : R_(R), r_(r), C_(C), Primitive(mat) {
	}

	bool intersect(Ray& ray) const {
		float3 O = ray.O_ - C_;
		float3 D = ray.D_;

		float sum_d_sqrd = D.x * D.x + D.y * D.y + D.z * D.z;
		float e = O.x * O.x + O.y * O.y + O.z * O.z - R_ * R_ - r_ * r_;
		float f = O.x * D.x + O.y * D.y + O.z * D.z;
		float four_a_sqrd = 4.0 * R_ * R_;

		float C0 = e * e - four_a_sqrd * (r_ * r_ - O.y * O.y);
		float C1 = 4.0 * f * e + 2 * four_a_sqrd * O.y * D.y;
		float C2 = 2.0 * sum_d_sqrd * e + 4.0 * f * f + four_a_sqrd * D.y * D.y;
		float C3 = 4.0 * sum_d_sqrd * f;
		float C4 = sum_d_sqrd * sum_d_sqrd;

		double C[5] = { C0,C1,C2,C3,C4 };
		double S[4] = { -1,-1,-1,-1 };
		int nIntersections = SolveQuartic(C, S);

		if (nIntersections == 0)
			return false;

		for (int i = 0; i < nIntersections; ++i) {
			float t = S[i];
			if ((t < ray.t_) && t > 0) {
				ray.t_ = t; ray.mat_ = mat_;
				float3  p = ray.O_ + ray.D_ * t;
				ray.N_ = computeNormal(O + D * t);
			};
		}
		return true;
	}

	float3 computeNormal(float3 p) const {
		float paramSquared = R_ * R_ + r_ * r_;

		float x = p.x;
		float y = p.y;
		float z = p.z;
		float sumSquared = x * x + y * y + z * z;

		float3 tmp = float3{
			(float)(4.0 * x * (sumSquared - paramSquared)),
			(float)(4.0 * y * (sumSquared - paramSquared + 2.0 * R_ * R_)),
			(float)(4.0 * z * (sumSquared - paramSquared)) };

		return normalize(tmp);
	}
};

class Sphere : public Primitive {
public:
	float3 C_;
	float r2_;

	Sphere(const float3 C, const float r2, Material* mat) : C_(C), r2_(r2), Primitive(mat) {
	}

	bool intersect(Ray& ray) const {
		float t0, t1;

		float3 L = ray.O_ - C_;
		float a = dot(ray.D_, ray.D_);
		float b = 2 * dot(ray.D_, L);
		float c = dot(L, L) - r2_;

		if (!solveQuadratic(a, b, c, t0, t1)) return false;
		if (t0 > t1) std::swap(t0, t1);

		if (t0 < 0) {
			t0 = t1; // if t0 is negative, let's use t1 instead 
			if (t0 < 0) return false; // both t0 and t1 are negative 
		}
		float t = t0;
		float3 N = normalize(C_ - (ray.O_ + t * ray.D_));
		if ((t < ray.t_) && t > 0) {
			ray.t_ = t; ray.mat_ = mat_; ray.N_ = -N;
			float u, v;
			get_sphere_uv(ray.O_ + t * ray.D_, C_, u, v);
			ray.u = u; ray.v = v;
		};
		//normal is the vector from the center of the circle, taking into the other direction so it points outwards
		return true;
	}
};

class Plane : public Primitive {
public:
	float3 N_;
	float3 P0_;

	Plane(const float3 N, const float3 P0, Material* mat) : N_(normalize(N)), P0_(P0), Primitive(mat) {
	}

	bool intersect(Ray& ray) const {
		float DdotN = dot(N_, ray.D_);//if dot is 0 the ray is paralel with the plane in thus no intersection
		if (DdotN > eps) return false;
		float3 P0O = P0_ - ray.O_;
		float t = dot(P0O, N_) / DdotN;
		if ((t < ray.t_) && t > (0)) { ray.t_ = t; ray.mat_ = mat_; ray.N_ = N_; }
		return true;
	}
};

class Triangle : public Primitive {
public:
	union { __m128 V04_; float3 V0_; };
	union { __m128 V14_; float3 V1_; };
	union { __m128 V24_; float3 V2_; };

	float3 N_;

	Triangle(const float3& V0, const float3& V1, const float3& V2, const float3& N, Material* mat) :
		V0_(V0), V1_(V1), V2_(V2), N_(N), Primitive(mat) {
	}

	bool intersect(Ray& ray) const {
		float3 V0V1 = V1_ - V0_;
		float3 V0V2 = V2_ - V0_;
		float3 Pvec = cross(ray.D_, V0V2);
		float det = dot(V0V1, Pvec);

		if (fabs(det) < eps) return false;
		float invDet = 1 / det;

		float3 Tvec = ray.O_ - V0_;
		float u = dot(Tvec, Pvec) * invDet;
		if (u < 0 || u > 1) return false;

		float3 Qvec = cross(Tvec, V0V1);
		float v = dot(ray.D_, Qvec) * invDet;
		if (v < 0 || u + v > 1) return false;

		float t = dot(V0V2, Qvec) * invDet;
		//float3 N = normalize(cross(V0V1, V0V2));
		if ((t < ray.t_) && t > 0) { ray.t_ = t; ray.mat_ = mat_; ray.N_ = N_; };

		return true;
	}
};