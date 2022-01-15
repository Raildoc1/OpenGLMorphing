#include "HarmonicMapper.h"

HarmonicMapper::HarmonicMapper(MeshData &source, MeshData &target)
{
	this->source = &source;
	this->target = &target;
}

void HarmonicMapper::init()
{
	initialized = true;
}
