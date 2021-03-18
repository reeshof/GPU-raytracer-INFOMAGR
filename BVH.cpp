#include "precomp.h"

//BVG
bool intersectAABB(AABB& b, Ray& r) {
	float tmin = -INFINITY, tmax = INFINITY;

	float tx1 = (b.bmin_[0] - r.O_.x) / r.D_.x;
	float tx2 = (b.bmax_[0] - r.O_.x) / r.D_.x;

	tmin = max(tmin, min(tx1, tx2));
	tmax = min(tmax, max(tx1, tx2));

	float ty1 = (b.bmin_[1] - r.O_.y) / r.D_.y;
	float ty2 = (b.bmax_[1] - r.O_.y) / r.D_.y;

	tmin = max(tmin, min(ty1, ty2));
	tmax = min(tmax, max(ty1, ty2));

	float tz1 = (b.bmin_[2] - r.O_.z) / r.D_.z;
	float tz2 = (b.bmax_[2] - r.O_.z) / r.D_.z;

	tmin = max(tmin, min(tz1, tz2));
	tmax = min(tmax, max(tz1, tz2));

	return (tmax >= tmin && tmax > 0);
}

std::pair<int, float> determineSplitPlane(AABB& boundingBox) {
	float xdiff = std::abs(boundingBox.bmax_[0] - boundingBox.bmin_[0]);
	float ydiff = std::abs(boundingBox.bmax_[1] - boundingBox.bmin_[1]);
	float zdiff = std::abs(boundingBox.bmax_[2] - boundingBox.bmin_[2]);

	//the split plane is parallel to the xy-plane if the zdiff is the highest, and the same in the other cases
	if (xdiff > ydiff && xdiff > zdiff)return { 0, (boundingBox.bmax_[0] + boundingBox.bmin_[0]) / 2.0f };
	if (ydiff > zdiff) return { 1,(boundingBox.bmax_[1] + boundingBox.bmin_[1]) / 2.0f };
	return { 2,(boundingBox.bmax_[2] + boundingBox.bmin_[2]) / 2.0f };
}

int getBinNr(const float C, const float& k0, const float& k1) {
	float binId = k1 * (C - k0);
	return binId;
}

std::pair<int, float> determineSplitPlane(AABB& boundingBox, std::vector<Triangle>& triangles, int first, int count) {
	auto result = determineSplitPlane(boundingBox);//binning on the dominant axis

	const int nbins = 8;//using 8 bins
	const int axis = result.first;
	bin bins[nbins];

	const float k1 = (nbins * (1 - eps)) / (boundingBox.bmax_[axis] - boundingBox.bmin_[axis]);
	const float k0 = boundingBox.bmin_[axis];

	for (int i = first; i < first + count; ++i) {
		TriangleBound tempt = TriangleBound(triangles[i]);
		int binNr = getBinNr(tempt.C_[axis], k0, k1);
		bins[binNr].addTriangle(tempt);
	}

	float Nl[nbins - 1] = { 0 };
	float Nr[nbins - 1] = { 0 };
	float Al[nbins - 1] = { 0 };
	float Ar[nbins - 1] = { 0 };

	Nl[0] = bins[0].N;
	Al[0] = bins[0].bounds_.calculateSA();
	AABB accBb = bins[0].bounds_;
	for (int i = 1; i < nbins - 1; ++i) {
		Nl[i] = bins[i].N + Nl[i - 1];
		accBb.adjustBounds(bins[i].bounds_);
		Al[i] = accBb.calculateSA();
	}

	Nr[6] = bins[7].N;
	Ar[6] = bins[7].bounds_.calculateSA();
	accBb = bins[7].bounds_;
	for (int i = 5; i >= 0; --i) {
		Nr[i] = bins[i + 1].N + Nr[i + 1];
		accBb.adjustBounds(bins[i + 1].bounds_);
		Ar[i] = accBb.calculateSA();
	}

	float C[nbins - 1] = { 0 };
	float minCost = std::numeric_limits<float>::max();
	int minIndex = -1;
	for (int i = 0; i < nbins - 1; ++i) {
		float cost = Al[i] * Nl[i] + Ar[i] * Nr[i];
		C[i] = cost;
		if (cost < minCost) {
			minCost = cost;
			minIndex = i;
		}
	}

	float stepSize = (boundingBox.bmax_[axis] - boundingBox.bmin_[axis]) / nbins;
	float splitValue = boundingBox.bmin_[axis] + stepSize * (minIndex + 1);

	return { axis,splitValue };
}

void showBounds(AABB& b) {
	for (const auto& min : b.bmin_)std::cout << min << " ";
	for (const auto& max : b.bmax_)std::cout << max << " ";
	std::cout << std::endl;
}

AABB CalculateBounds(std::vector<Triangle>& triangles, int first, int count) {
	AABB boundingbox;
	boundingbox.setMinMax();

	for (int i = first; i < first + count; ++i) {
		boundingbox.adjustBounds(triangles[i]);
	}

	return boundingbox;
}

void BVHNodeConstructor::Intersect(Ray& ray) {
	pool_[0].traverse(ray, pool_, triangles_);
}

void BVHNodeConstructor::ConstructBVH() {

	int rootIndex = 0;
	int poolPtr = 2;

	pool_[rootIndex].first_ = 0;
	pool_[rootIndex].count_ = triangles_.size();
	pool_[rootIndex].bounds_ = CalculateBounds(triangles_, pool_[rootIndex].first_, pool_[rootIndex].count_);
	pool_[rootIndex].SubDivide(poolPtr, pool_, triangles_);

	std::cout << " - n of triangles: " << triangles_.size() << "\n";
	std::cout << " - n of nodes: " << poolPtr - 1 << "\n";
}

int BVHNode::Partition(std::vector<Triangle>& triangles, int axis, float partitionValue) {
	auto it = std::partition(triangles.begin() + first_, triangles.begin() + first_ + count_, [&](Triangle& t) {
		if (axis == 0) { float centroidValue = (t.V0_.x + t.V1_.x + t.V2_.x) / 3.0f; return centroidValue < partitionValue; }
		if (axis == 1) { float centroidValue = (t.V0_.y + t.V1_.y + t.V2_.y) / 3.0f; return centroidValue < partitionValue; }
		if (axis == 2) { float centroidValue = (t.V0_.z + t.V1_.z + t.V2_.z) / 3.0f; return centroidValue < partitionValue; }
		return true;
		});
	int dist = std::distance(triangles.begin() + first_, it);
	return dist;
}

void BVHNode::SubDivide(int& poolptr, std::vector<BVHNode>& pool, std::vector<Triangle>& triangles) {
	if (count_ < 3) {
		return;
	}
	left_ = poolptr++;
	right_ = poolptr++;

	auto [splitplane, splitvalue] = determineSplitPlane(bounds_, triangles, first_, count_);
	//auto [splitplane, splitvalue] = determineSplitPlane(bounds_); use this if no binning required
	int countFirst = Partition(triangles, splitplane, splitvalue);

	//to prevent infinite loops
	if (countFirst == 0 || count_ - countFirst == 0) {
		bounds_ = CalculateBounds(triangles, first_, count_);
		return;
	}

	pool[left_].first_ = first_;
	pool[left_].count_ = countFirst;
	pool[left_].bounds_ = CalculateBounds(triangles, pool[left_].first_, pool[left_].count_);

	pool[right_].first_ = first_ + countFirst;
	pool[right_].count_ = count_ - countFirst;
	pool[right_].bounds_ = CalculateBounds(triangles, pool[right_].first_, pool[right_].count_);

	pool[left_].SubDivide(poolptr, pool, triangles);
	pool[right_].SubDivide(poolptr, pool, triangles);

	isLeaf_ = false;
}

void BVHNode::traverse(Ray& ray, std::vector<BVHNode>& pool, std::vector<Triangle>& triangles) {
	if (!intersectAABB(bounds_, ray))return;
	if (isLeaf_) {
		for (int i = first_; i < first_ + count_; ++i) {
			triangles[i].intersect(ray);
		}
	}
	else {
		pool[left_].traverse(ray, pool, triangles);
		pool[left_ + 1].traverse(ray, pool, triangles);
	}
}

BVHutil* convertBVH(std::vector<BVHNode> BVH)
{
	BVHutil* newBVH = new BVHutil[BVH.size()];
	for (int i = 0; i < BVH.size(); ++i) {
		newBVH[i].minBounds_[0] = BVH[i].bounds_.bmin_[0];
		newBVH[i].minBounds_[1] = BVH[i].bounds_.bmin_[1];
		newBVH[i].minBounds_[2] = BVH[i].bounds_.bmin_[2];

		newBVH[i].maxBounds_[0] = BVH[i].bounds_.bmax_[0];
		newBVH[i].maxBounds_[1] = BVH[i].bounds_.bmax_[1];
		newBVH[i].maxBounds_[2] = BVH[i].bounds_.bmax_[2];

		newBVH[i].left_ = BVH[i].left_;
		newBVH[i].right_ = BVH[i].right_;
		newBVH[i].first_ = BVH[i].first_;
		newBVH[i].count_ = BVH[i].count_;
		newBVH[i].isLeaf_ = BVH[i].isLeaf_;
	}

	return newBVH;
}