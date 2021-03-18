#pragma once

//need this to get the BVH from cpu to gpu.
//couldnt use the BVH used their because of the datatypes that conflicted with the cuda namespace
struct BVHutil
{
	float minBounds_ [3];
	float maxBounds_ [3];
	bool isLeaf_ = true;
	int left_, right_;
	int first_, count_;
};

