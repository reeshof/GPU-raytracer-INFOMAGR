#pragma once
//BVG
struct AABB {
	union { float3 mins_; __m128 mins4_; float bmin_[4]; };
	union { float3 maxs_; __m128 maxs4_; float bmax_[4]; };

	AABB() {
	}

	AABB(Triangle& T) {
		setMinMax();
		adjustBounds(T);
	}
	void setMinMax() {
		for (int i = 0; i < 3; ++i)bmin_[i] = 99999999.0f;
		for (int i = 0; i < 3; ++i)bmax_[i] = -99999999.0f;
	}
	void adjustBounds(const Triangle& T) {
		mins4_ = _mm_min_ps(mins4_, T.V04_);//
		mins4_ = _mm_min_ps(mins4_, T.V14_);//
		mins4_ = _mm_min_ps(mins4_, T.V24_);//

		maxs4_ = _mm_max_ps(maxs4_, T.V04_);
		maxs4_ = _mm_max_ps(maxs4_, T.V14_);
		maxs4_ = _mm_max_ps(maxs4_, T.V24_);
	}

	void adjustBounds(AABB& bb) {
		mins4_ = _mm_min_ps(mins4_, bb.mins4_);//
		maxs4_ = _mm_max_ps(maxs4_, bb.maxs4_);//
	}

	float calculateSA() {
		union { __m128 e4; float e[4]; };
		e4 = _mm_sub_ps(maxs4_, mins4_);
		return (e[0] * e[1] + e[0] * e[2] + e[1] * e[2]);//actual surface area * 2 but its a constant so doesnt matter
	}
};

struct TriangleBound {
	AABB bb_;
	float C_[3];
	TriangleBound() {}
	TriangleBound(Triangle& T) {
		bb_ = AABB(T);

		C_[0] = (T.V0_.x + T.V1_.x + T.V2_.x) / 3.0f;
		C_[1] = (T.V0_.y + T.V1_.y + T.V2_.y) / 3.0f;
		C_[2] = (T.V0_.z + T.V1_.z + T.V2_.z) / 3.0f;
	}
};

struct bin {
	AABB bounds_;
	int N = 0;//n of triangles in bin

	bin() {
		bounds_ = AABB();
		bounds_.setMinMax();
	}

	void addTriangle(TriangleBound& T) {
		bounds_.adjustBounds(T.bb_);
		++N;
	}
};

struct BVHNode {
	AABB bounds_;
	bool isLeaf_ = true;
	int left_, right_;//child indices
	int first_, count_;//indice of first triangle, and number of triangles

	void SubDivide(int& poolptr, std::vector<BVHNode>& pool, std::vector<Triangle>& triangles);
	int Partition(std::vector<Triangle>& triangles, int axis, float partitionValue);
	void BVHNode::traverse(Ray& ray, std::vector<BVHNode>& pool, std::vector<Triangle>& triangles);
};

class BVHNodeConstructor {
public:
	std::vector<Triangle> triangles_;
	std::vector<BVHNode> pool_;
	BVHNodeConstructor() {
		
	};

	void initializePool() {
		int N = triangles_.size();
		pool_ = std::vector<BVHNode>(N * 2 - 1);
	}

	void ConstructBVH();
	void Intersect(Ray& ray);
};



std::pair<int, float> determineSplitPlane(AABB& boundingBox);
int getBinNr(const float C, const float& k0, const float& k1);
bool intersectAABB(AABB& b, Ray& r);
std::pair<int, float> determineSplitPlane(AABB& boundingBox, std::vector<Triangle>& triangles, int first, int count);
void showBounds(AABB& b);
AABB CalculateBounds(std::vector<Triangle>& triangles, int first, int count);
BVHutil* convertBVH(std::vector<BVHNode> BVH);
