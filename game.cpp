#include "precomp.h" // include (only) this in every .cpp file

//https://github.com/tinyobjloader/tinyobjloader
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


std::ostream& operator<<(std::ostream& os, const BVHNode& dt) {
	os << "--------------BVH--------------\n";
	for (int i = 0; i < 3; ++i)
		os << dt.bounds_.bmin_[i] << " ";
	os << "\n";
	for (int i = 0; i < 3; ++i)
		os << dt.bounds_.bmax_[i] << " ";
	os << "\n";
	os << dt.isLeaf_ << " ";
	os << dt.left_ << " " << dt.right_ << "\n";
	os << dt.first_ << " " << dt.count_;
	os << "----------------------------------\n";
	return os;
}

//Some utility functionst
inline float Radians(float deg) {
	return (M_PI / 180) * deg;
}

union minsmaxs {
	float3 f;
	__m128 f4_;
};

ImageTexture::ImageTexture(const std::string& filename) {
	int nn = 0;
	data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
}

void showControls() {
	std::string outputString = "w/a/s/d: moving the camera\n1/2: zoom in/out\nLeft/right arrow: rotate left/right\nup/down arrow: move camera up/down\nf: use fish eye lense\ng: use path tracer\nh: use gpu for path tracing\n";
	std::cout << outputString << std::endl;
}

ostream& operator<<(ostream& os, const float3& dir)
{
	os << dir.x << " " << dir.y << " " << dir.z << "\n";
	return os;
}

float3 getRandomDirection(const float3& N) {
	float3 randVec = normalize(float3{ Rand(2) - 1,Rand(2) - 1,Rand(2) - 1 });
	float theDot = dot(randVec, N);

	while (theDot < 0) {
		randVec = normalize(float3{ Rand(2) - 1,Rand(2) - 1,Rand(2) - 1 });
		theDot = dot(randVec, N);
	}

	return randVec;
}

float3 getRandomDirection2(float3& N, float r1, float r2) {
	float3 W;
	if (fabs(N.x) > 0.99)
		W = float3{ 0,1,0 };
	else
		W = float3{ 1,0,0 };
	float3 T = normalize(cross(N, W));
	float3 B = normalize(cross(T, N));

	//x,y, and z are random points on the unit hemisphere
	float x = cos(2 * M_PI * r1) * sqrtf(1.0f - r2 * r2);
	float y = sin(2 * M_PI * r1) * sqrtf(1.0f - r2 * r2);
	float z = r2;
	float3 P = normalize(float3{ x, y, z });

	float Px = P.x * T.x + P.y * B.x + P.z * N.x;
	float Py = P.x * T.y + P.y * B.y + P.z * N.y;
	float Pz = P.x * T.z + P.y * B.z + P.z * N.z;

	return normalize(float3{ Px,Py,Pz });
}

void get_sphere_uv(const float3& p, const float3 c, float& u, float& v) {
	float3 tp = normalize(p - c);
	float phi = atan2(tp.z, tp.x);
	float theta = asin(tp.y);
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / M_PI;
}


bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) x0 = x1 = -0.5 * b / a;
	else {
		float q = (b > 0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1) std::swap(x0, x1);

	return true;
}

class Light {
public:
	float3 loc_;
	float3 color_;
	float intensity_;

	Light(float3 loc, float3 color, float intensity) : loc_(loc), color_(color), intensity_(intensity) {}
	virtual float getFallof(const float3& rayOrigin) const = 0;
};

class NormalLight : public Light {
public:
	NormalLight(float3 loc, float3 color, float intensity) : Light(loc, color, intensity) {}
	float getFallof(const float3& rayOrigin) const {
		return 1.0f;
	}
};

class SpotLight : public Light {
public:
	float3 direction_;
	float totalWidth_, fallofWidth_;
	SpotLight(float3 loc, float3 color, float intensity, float totalWidth, float fallofWidth, float3 direction) :
		Light(loc, color, intensity), totalWidth_(std::cos(Radians(totalWidth))),
		fallofWidth_(std::cos(Radians(fallofWidth))), direction_(normalize(direction)) {

	}

	float getFallof(const float3& rayOrigin) const {
		float cosTheta = dot(normalize(rayOrigin - loc_), direction_);
		if (cosTheta < totalWidth_) return 0;
		if (cosTheta > fallofWidth_) return 1;
		float delta = (cosTheta - totalWidth_) / (fallofWidth_ - totalWidth_);
		return (delta * delta) * (delta * delta);
	}
};



//Camera
class Camera {
public:
	float3 E_;
	float3 V_;
	float3 C_;
	float d_;
	float3 P0_;
	float3 P1_;
	float3 P2_;

	//the float3 Pt_ are the initial positions transformed with the camera matrix.
	float3 Et_;
	float3 P0t_;
	float3 P1t_;
	float3 P2t_;
	mat4 cameraMatrix_{};

	Camera(const float3& E, const float3& V, const float d) : E_(E), V_(V), d_(d) {
		setInitialPoints();
	}

	void setInitialPoints() {
		C_ = E_ + d_ * V_;
		P0_ = C_ + float3{ -1,-1,0 };
		P1_ = C_ + float3{ 1,-1,0 };
		P2_ = C_ + float3{ -1,1,0 };
	}

	void zoom(const float increase) {
		d_ += increase;
		setInitialPoints();
	}

	void transformInitialPoints() {
		Et_ = cameraMatrix_.TransformPoint(E_);
		P0t_ = cameraMatrix_.TransformPoint(P0_);
		P1t_ = cameraMatrix_.TransformPoint(P1_);
		P2t_ = cameraMatrix_.TransformPoint(P2_);
	}

	Ray getRay(const float u, const float v) {
		float3 Puv = P0t_ + u * (P1t_ - P0t_) + v * (P2t_ - P0t_);
		Puv -= Et_;//ray from camera origin to point on world screen
		Puv /= length(Puv);//normalized
		return { Et_,Puv };
	}

	Ray getFishRay(const float i, const float j, const int width, const int height, const float aperture) {
		float x = 2 * i / (float)width - 1;
		float y = 2 * j / (float)height - 1;
		float r = std::sqrt(x * x + y * y);
		if (r > 1)return { Et_,float3{-1,-1,-1} };
		float phi;
		if (r == 0) phi = 0;
		else if (x < 0) phi = M_PI - std::asin(y / r);
		else phi = asin(y / r);
		float theta = r * Radians(aperture) / 2.0f;
		float xdir = std::sin(theta) * std::cos(phi);
		float ydir = std::sin(theta) * std::sin(phi);
		float zdir = std::cos(theta);
		return { Et_,cameraMatrix_.TransformVector(normalize(float3{xdir,ydir,zdir})) };
	}

	void transformPoints() {
		E_ = cameraMatrix_.TransformPoint(E_);
		P0_ = cameraMatrix_.TransformPoint(P0_);
		P1_ = cameraMatrix_.TransformPoint(P1_);
		P2_ = cameraMatrix_.TransformPoint(P2_);
	}

	void translate(const float x, const float y, const float z) {
		const auto translation = mat4::Translate(x, y, z);
		cameraMatrix_ = cameraMatrix_ * translation;

	}

	void rotateY(const float a) {
		const auto rotation = mat4::RotateY(a);
		cameraMatrix_ = cameraMatrix_ * rotation;
	}

	void rotateX(const float a) {
		const auto rotation = mat4::RotateX(a);
		cameraMatrix_ = cameraMatrix_ * rotation;
	}
};

//returns the refracted direction and the ratio of refrected light in Fr.
float3 refract(const float3& I, const float3& N, const float& ior, float& Fr) {
	float cosi = clamp(-1.0f, 1.0f, dot(I, N));
	float etai = 1, etat = ior;
	float3 n = N;
	if (cosi < 0) { cosi = -cosi; }
	else { std::swap(etai, etat); n = -1 * N; }

	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	if (sint >= 1) {
		Fr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		Fr = (Rs * Rs + Rp * Rp) / 2;
	}

	float eta = etai / etat;
	float k = 1 - eta * eta * (1 - cosi * cosi);

	return k < 0 ? float3{ 0, 0, 0 } : eta * I + (eta * cosi - sqrtf(k)) * n;
}

struct bvhobject {
	BVHNodeConstructor* BVH;
	std::vector<Triangle> triangles_;
};

class Scene {
public:
	std::vector<Primitive*> primitives_;
	std::vector<Triangle> roomTriangles;
	trianglec* trianglesc_;//for cuda

	std::vector<bvhobject> objects;

	BVHNodeConstructor* BVHroom;
	std::vector<Light*> lights_;//used for whitted style ray tracing, for path tracing the light is a primitive with emittance

	Camera camera_;

	float3* frameBuffer_;//Holds all the pixel values accumelated over several frames, destroyed by Game.
	float nFrames_ = 0;

	bool fish = false;
	bool path = false;
	bool useGpu = false;

	//Camera always initialized in the center of the screen
	Scene(float3* frameBuffer) : frameBuffer_(frameBuffer), camera_(Camera({ 0,0,0 }, { 0,0,1 }, 1)) {
	}

	~Scene() {
		for (auto p : primitives_)
			delete p;
		for (auto p : lights_)
			delete p;
	}

	void resetBuffer() {
		memset(frameBuffer_, 0, SCRWIDTH * SCRHEIGHT * 3 * 4);
		nFrames_ = 0;
	}

	void addLight(Light* l) {
		lights_.push_back(l);
	}
	void addPrimitive(Primitive* p) {
		primitives_.push_back(p);
	}

	//Reads object file and add all the triangles to triangles_
	void readObjFile(std::string& filename, float3 trans, float rotX, float rotY, float scale, Material* objectMat, int matIndex) {
		std::string inputfile = filename;
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string warn, err;
		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str())) {
			throw std::runtime_error(warn + err);
		}

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str())) {
			throw std::runtime_error(warn + err);
		}

		std::vector<Material*> listOfMats;

		for (const auto& mat : materials) {
			if (!mat.diffuse_texname.empty()) {
				ImageTexture* imgText = new ImageTexture(mat.diffuse_texname);
				Material* newMat = new Material(imgText, MaterialType::diffuse, 0, 1);
				listOfMats.push_back(newMat);
			}
		}

		bvhobject newObject;
		newObject.BVH = new BVHNodeConstructor();

		const auto tMat = mat4::Translate(trans);
		const auto rMat = mat4::RotateX(rotX) * mat4::RotateY(rotY);
		const auto sMat = mat4::Scale(scale);

		for (size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces(polygon)
			size_t index_offset = 0;

			for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				int fv = shapes[s].mesh.num_face_vertices[f];
				int materialId = shapes[s].mesh.material_ids[f];
				//std::cout << triangles_.size() << std::endl;
				std::vector<float3> V(3);

				// Loop over vertices in the face.
				for (size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
					tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
					tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];

					//tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
					//tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
					//tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
					//N = float3{ nx,ny,nz };
					//V[v] = tMat.TransformPoint(sMat.TransformPoint({ vx,vy,vz }));
					V[v] = tMat.TransformPoint(rMat.TransformPoint(sMat.TransformPoint({ vx,vy,vz })));

				}
				index_offset += fv;

				float3 N = normalize(cross(V[0] - V[1], V[0] - V[2]));
				Triangle t = Triangle(V[0], V[1], V[2], N, objectMat);
				newObject.BVH->triangles_.push_back(std::move(t));
			}
		}

		std::cout << "-------------------------------------------------\n";
		Timer timer;
		newObject.BVH->initializePool();
		float elapsedTime1 = timer.elapsed();
		std::cout << "Initialized BVH pool in time: " << elapsedTime1 << "s\n";
		timer.reset();
		newObject.BVH->ConstructBVH();
		float elapsedTime2 = timer.elapsed();

		std::cout << "Created BVH in time: " << elapsedTime2 << "s\n";
		std::cout << "total time: " << elapsedTime2 + elapsedTime1 << "s\n";
		std::cout << "-------------------------------------------------\n";

		int N = newObject.BVH->triangles_.size();
		trianglesc_ = new trianglec[N];

		for (int i = 0; i < N; ++i) {
			trianglec temp;

			temp.V0_[0] = newObject.BVH->triangles_[i].V0_.x;
			temp.V0_[1] = newObject.BVH->triangles_[i].V0_.y;
			temp.V0_[2] = newObject.BVH->triangles_[i].V0_.z;

			temp.V1_[0] = newObject.BVH->triangles_[i].V1_.x;
			temp.V1_[1] = newObject.BVH->triangles_[i].V1_.y;
			temp.V1_[2] = newObject.BVH->triangles_[i].V1_.z;

			temp.V2_[0] = newObject.BVH->triangles_[i].V2_.x;
			temp.V2_[1] = newObject.BVH->triangles_[i].V2_.y;
			temp.V2_[2] = newObject.BVH->triangles_[i].V2_.z;

			temp.N_[0] = newObject.BVH->triangles_[i].N_.x;
			temp.N_[1] = newObject.BVH->triangles_[i].N_.y;
			temp.N_[2] = newObject.BVH->triangles_[i].N_.z;

			trianglesc_[i] = temp;
		}

		initializeBVH(convertBVH(newObject.BVH->pool_), newObject.BVH->pool_.size());
		initializeTriangles(trianglesc_, N, matIndex);

		objects.push_back(newObject);

		initializeRayBuffer(SCRWIDTH, SCRHEIGHT);
	}

	void gpuTrace() {
		int N = objects[0].triangles_.size();
		float P0t[3]{ camera_.P0t_.x,camera_.P0t_.y,camera_.P0t_.z };
		float P1t[3]{ camera_.P1t_.x,camera_.P1t_.y,camera_.P1t_.z };
		float P2t[3]{ camera_.P2t_.x, camera_.P2t_.y, camera_.P2t_.z };
		float Et[3]{ camera_.Et_.x, camera_.Et_.y, camera_.Et_.z };


		generateRays(P0t, P1t, P2t, Et);
		traceRays(N);
		float* fb = getFrameBuffer();
		for (int i = 0, j = 0; i < SCRWIDTH * SCRHEIGHT; ++i, j += 3) {
			frameBuffer_[i].x += fb[j];
			frameBuffer_[i].y += fb[j + 1];
			frameBuffer_[i].z += fb[j + 2];
		}
		delete fb;
	}

	//For whitted style ray tracing
	float3 directIllumination(Ray& ray) {
		float3 totalLight{ 0,0,0 };

		for (const auto& light : lights_) {
			float3 intersectionPoint = (ray.O_ + (ray.t_) * ray.D_) + ray.N_ * eps;
			float3 rayDirection = light->loc_ - intersectionPoint;
			float rayLength = length(rayDirection);
			rayDirection = normalize(rayDirection);

			Ray shadowRay{ intersectionPoint,rayDirection };

			for (const auto& primitive : primitives_) {
				primitive->intersect(shadowRay);
			}

			for (const auto& obj : objects)
				obj.BVH->Intersect(ray);

			if (shadowRay.t_ < 0 || shadowRay.t_ > rayLength) {//not occluded
				const auto newColor = light->color_ * light->getFallof(intersectionPoint) * light->intensity_ * (dot(ray.N_, rayDirection));
				if (newColor.x >= 0 && newColor.y >= 0 && newColor.z >= 0)
					totalLight += newColor;
			}
		}
		return totalLight;
	}

	void traceRay(Ray& ray) {
		for (const auto& primitive : primitives_) primitive->intersect(ray);
		for (const auto& obj : objects) obj.BVH->Intersect(ray);
	}

	//Skydome
	Primitive* skyDome = new Sphere({ 0,0,0 }, 10000, new Material{ new ImageTexture("assets/forest.ppm"),MaterialType::diffuse,0,1.3 });
	float3 getSkyDomeColor(Ray& ray) {
		skyDome->intersect(ray);
		return skyDome->mat_->getColor(ray.u, ray.v, { 0,0,0 });
	}

	//Path tracing
	float3 pathTrace(Ray& ray, int recursionDepth) {
		if (recursionDepth >= 15)  return { 0,0,0 };

		traceRay(ray);
		if (ray.t_ == std::numeric_limits<float>::max()) return getSkyDomeColor(ray);// getSkyDomeColor(ray);

		if (ray.mat_->matType_ == MaterialType::light)return ray.mat_->emittance_;

		float3 I = (ray.O_ + (ray.t_) * ray.D_);//intersectionPoint + ray.N_*eps

		if (ray.mat_->matType_ == MaterialType::dialectric) {
			float Fr = -1;//amount of light reflected
			float3 refractionDirection = refract(ray.D_, ray.N_, ray.mat_->index_, Fr);

			bool outside = dot(ray.D_, ray.N_) < 0;
			float3 n = ray.N_;
			if (!outside)
				n = -1 * ray.N_;

			float randomNumber = Rand(1.0f);
			if (randomNumber < (1 - Fr)) {
				Ray refrectiveRay{ I + n * -eps,normalize(refractionDirection) };
				return pathTrace(refrectiveRay, recursionDepth + 1);
			}
			else {
				Ray reflectiveRay{ I + n * eps, reflect(ray.D_, n) };
				return pathTrace(reflectiveRay, recursionDepth + 1);
			}
		}

		if (ray.mat_->reflectivity_ > 0) {
			float randomNumber = Rand(1.0f);
			if (randomNumber < ray.mat_->reflectivity_) {
				Ray reflectedRay{ I + ray.N_ * eps,reflect(ray.D_, ray.N_) };
				return ray.mat_->getColor(ray.u, ray.v, I) * pathTrace(reflectedRay, recursionDepth + 1);
			}
		}


		float3 R = getRandomDirection2(ray.N_, Rand(1), Rand(1)) + ray.N_ * eps;//get a random direction on the hemisphere
		Ray newRay(I, R);

		float3 BRDF = ray.mat_->getColor(ray.u, ray.v, I) / M_PI;
		float3 Ei = pathTrace(newRay, recursionDepth + 1) * dot(ray.N_, R);

		return (M_PI * 2.0f * BRDF * Ei);
	}

	//Whitted style ray tracing
	float3 whittedTrace(Ray& ray, int recursionDepth) {
		if (recursionDepth >= 5) {
			return{ 0,0,0 };
		}
		for (const auto& primitive : primitives_) {
			primitive->intersect(ray);
		}

		for (auto& obj : objects) {
			obj.BVH->Intersect(ray);
		}

		if (!(ray.t_ == std::numeric_limits<float>::max())) {
			float3 directColor{ 0,0,0 };
			float3 reflectiveColor{ 0,0,0 };
			float3 intersectionPoint = (ray.O_ + (ray.t_) * ray.D_);
			if (ray.mat_->matType_ == MaterialType::dialectric) {
				float Fr = -1;//amount of light reflected
				float3 refractionDirection = refract(ray.D_, ray.N_, ray.mat_->index_, Fr);

				bool outside = dot(ray.D_, ray.N_) < 0;
				float3 n = ray.N_;
				if (!outside)
					n = -1 * ray.N_;

				float3 newDirectColor{ 0,0,0 };
				float3 newReflectiveColor{ 0,0,0 };
				if (Fr > 0)
					newReflectiveColor = Fr * whittedTrace(Ray{ intersectionPoint + n * eps, reflect(ray.D_, n) }, recursionDepth + 1);
				if (Fr < 1)
					newDirectColor = (1 - Fr) * whittedTrace(Ray{ intersectionPoint + n * -eps,normalize(refractionDirection) }, recursionDepth + 1);

				//Beer's law
				if (!outside) {//this ray travelled through an object, the recursive call gives the color obtained from outside the object
					//which is what needs to be scaled. The mat color is used as the absorbtion values
					newDirectColor.x *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).x * ray.t_);
					newDirectColor.y *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).y * ray.t_);
					newDirectColor.z *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).z * ray.t_);
					newReflectiveColor.x *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).x * ray.t_);
					newReflectiveColor.y *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).y * ray.t_);
					newReflectiveColor.z *= exp(-ray.mat_->getColor(0, 0, intersectionPoint).z * ray.t_);
				}

				directColor += newDirectColor;
				reflectiveColor += newReflectiveColor;
			}
			else {
				if (ray.mat_->reflectivity_ > 0) {
					float3 reflectedDirection = reflect(ray.D_, ray.N_);
					reflectiveColor = ray.mat_->reflectivity_ * whittedTrace(Ray{ intersectionPoint + ray.N_ * eps, reflectedDirection }, recursionDepth + 1);
				}
				if (ray.mat_->reflectivity_ < 1) {
					directColor = (1 - ray.mat_->reflectivity_) * ray.mat_->getColor(ray.u, ray.v, intersectionPoint) * directIllumination(ray);
				}
			}
			return directColor + reflectiveColor;
		}
		return getSkyDomeColor(ray);
	}
};

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
Scene* scene;

void Game::Init()
{
	// allocate memory for the floating point render target
	floatBuffer = new float3[SCRWIDTH * SCRHEIGHT];
	// clear the render target
	memset(floatBuffer, 0, SCRWIDTH * SCRHEIGHT * 4);
	scene = new Scene(floatBuffer);
	scene->camera_.translate(0, -15, 0);
	//nef.obj has 2m triangles, dragon.obj about 250k
	Material* newMat = new Material{ {0.8, 0.5, 0.7},MaterialType::diffuse,0,1.3 };
	Material* newMat2 = new Material{ { 0.5, 0.6, 0.8 },MaterialType::reflective,0.8,1.3 };

	//texture used for skydome in gpu application
	ImageTexture* skydometexture = new ImageTexture("assets/forest.ppm");
	float skydomecenter[3]{ 0,0,0 };
	intializeSkydome(skydometexture->data, skydometexture->nx, skydometexture->ny, skydomecenter, 1000000);
	
	initializeObjs();

	float alb[3] = { 0.8, 0.5, 0.7 };
	addMat(alb, 0, false, false, 1.0f, 0);
	scene->readObjFile(std::string("assets/dragon.obj"), { 25, -5, 75 }, 3.0, 3.7, 0.6, newMat,0);

	float alb2[3] = { 0.5, 0.6, 0.8 };
	addMat(alb2,0.8, false, false, 1.0f, 1);
	scene->readObjFile(std::string("assets/CrumpledDevelopable.obj"), { 10, 20, 75 }, 3.0, 3.0, 65, newMat2,1);
	
	allocateMats();

	//lights for whitted style
	scene->addLight(new NormalLight({ 0,-25,3 }, { 1,1,1 }, 1.1));
	scene->addLight(new NormalLight({ 0,0,0 }, { 1,1,1 }, 1.1));
	scene->addLight(new SpotLight({ 17,-10,8 }, { 1,1,1 }, 1.4, 45, 20, { 0,1,0 }));
	scene->addLight(new SpotLight({ 21,-8,18 }, { 1,1,1 }, 1.4, 45, 20, { 0,1,0 }));
	scene->addLight(new SpotLight({ 27,-10,8 }, { 1,1,1 }, 1.4, 45, 20, { 0,1,0 }));

	//light for path tracing scene
	scene->addPrimitive(new Sphere({ 0,-100,-100 }, 10000, new Material{ {1,1,1},MaterialType::light,0,1.3,{3.5,3.5,3.5} }));

	showControls();
}



// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
	delete[] floatBuffer;
}

// -----------------------------------------------------------
// Convert the floating point pixel buffer to integer pixels
// -----------------------------------------------------------
void Game::PresentFloatBuffer(float nFrames)
{
	int pixelIndex = 0;
	uint* screenBuffer = screen->GetBuffer();
	for (int y = 0; y < SCRHEIGHT; y++) for (int x = 0; x < SCRWIDTH; x++)
	{
		float3 pixel = floatBuffer[pixelIndex] / nFrames;
		int red = clamp((int)(pixel.x * 256.0f), 0, 255);
		int green = clamp((int)(pixel.y * 256.0f), 0, 255);
		int blue = clamp((int)(pixel.z * 256.0f), 0, 255);
		screenBuffer[pixelIndex++] = (red << 16) + (green << 8) + blue;
	}
}

// Handle input
void Game::KeyDown(int key) {
	switch (key) {
	case 65: scene->camera_.translate(-1, 0, 0); scene->resetBuffer(); break;
	case 68: scene->camera_.translate(1, 0, 0); scene->resetBuffer(); break;
	case 87: scene->camera_.translate(0, 0, 1); scene->resetBuffer(); break;
	case 83: scene->camera_.translate(0, 0, -1); scene->resetBuffer(); break;
	case 263: scene->camera_.rotateY(-0.1); scene->resetBuffer(); break;
	case 262: scene->camera_.rotateY(0.1); scene->resetBuffer(); break;
	case 49: scene->camera_.zoom(0.1); scene->resetBuffer(); break;
	case 50: scene->camera_.zoom(-0.1); scene->resetBuffer(); break;
	case 265: scene->camera_.translate(0, -1, 0); scene->resetBuffer(); break;
	case 264: scene->camera_.translate(0, 1, 0); scene->resetBuffer(); break;
	case 70: scene->fish = !scene->fish; scene->resetBuffer(); break;
	case 71: scene->path = !scene->path; scene->resetBuffer(); break;
	case 72: scene->useGpu = !scene->useGpu; scene->resetBuffer(); break;
	}
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick(float deltaTime)
{
	Timer timer;

	scene->camera_.transformInitialPoints();

	float aspectRatio = (float)SCRWIDTH / (float)SCRHEIGHT;
	scene->nFrames_++;

	float xStepSize = aspectRatio / SCRWIDTH;
	float yStepSize = 1 / SCRHEIGHT;

	bool fish = scene->fish;
	bool path = scene->path;
	bool useGpu = scene->useGpu;

	if(useGpu)
		scene->gpuTrace();
	else {
		Ray newRay(float3{ -1,-1,-1 }, float3{ -1,-1,-1 });

		for (int y = 0; y < SCRHEIGHT; y++) for (int x = 0; x < SCRWIDTH; x++)
		{
			float u = (float)x / (float)SCRHEIGHT + Rand(yStepSize);// *aspectRatio;
			float v = (float)y / (float)SCRWIDTH * aspectRatio + Rand(xStepSize);

			if (fish) newRay = scene->camera_.getFishRay(x, y, SCRWIDTH, SCRHEIGHT, 180);
			else newRay = scene->camera_.getRay(u, v);

			float3 color;
			if (newRay.D_ == float3{ -1,-1,-1 })color = { 0,0,0 };
			else {
				if (path) color = scene->pathTrace(newRay, 0);
				else color = scene->whittedTrace(newRay, 0);
			}

			floatBuffer[x + y * SCRWIDTH] += color;
		}
	}

	PresentFloatBuffer(scene->nFrames_);

	std::string elapsedTime = "Frame time: " + to_string(timer.elapsed());
	std::string nFrames = "Frame samples: " + to_string((int)scene->nFrames_);
	const char* timestring = elapsedTime.c_str();
	const char* frametring = nFrames.c_str();
	screen->Print(timestring, 4, 4, 0xffffff);
	screen->Print(frametring, 4,16, 0xffffff);
}