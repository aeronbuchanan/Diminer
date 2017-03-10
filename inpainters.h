/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of Diminer, a digital inpainting resource.
 *
 * Diminer is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Diminer is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Diminer.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <math.h>
#include <limits>
#include <tuple>

#include "diminer.h"
#include "patch.h"

// DEBUG
#include "CImg.h" 

namespace Diminer
{

template <class M, class T>
class SpatialGridSquare
{
public:
	SpatialGridSquare() : xmax(0), xmin(0), ymax(0), ymin(0), contents(0) {};

	M sqrdDistanceTo(M x, M y)
	{
		double dx = 0;
		double dy = 0;

		if ( x < xmin ) dx = x - xmin;
		else if ( x > xmax ) dx = x - xmax;

		if ( y < ymin ) dy = y - ymin;
		else if ( y > ymax ) dy = y - ymax;

		return (M)(dx * dx + dy * dy);
	}

	M xmax, xmin, ymax, ymin;
	T* contents;
};

template <class M, class T>
using SpatialGrid = std::vector<SpatialGridSquare<M, T>>;

typedef std::pair<BoundaryColors::iterator, BoundaryColors::iterator> BoundaryColorRange;
typedef std::vector<BoundaryColorRange> BoundaryColorRanges;

typedef SpatialGridSquare<double, BoundaryColorRanges> BoundaryColorGridSq;
typedef SpatialGrid<double, BoundaryColorRanges> BoundaryColorGrid;

class Inpainter
{
public:
	Inpainter(BoundaryColors const * const);
	Inpainter(BoundaryColors const * const _b, CImg<uchar> * _img) : m_boundary(_b), debug(_img) {}
	~Inpainter() {}

	BoundaryColorRanges withinRange(CoordPtr);
	virtual Color pixelColor(CoordPtr);

protected:
	BoundaryColors const * m_boundary;
	BoundaryColorGrid m_grid;

	CImg<uchar> * debug;
};

class BleedInpainter : public Inpainter
{
public:
	BleedInpainter(BoundaryColors const * const _b) : Inpainter(_b) {}
	BleedInpainter(BoundaryColors const * const _b, CImg<uchar> * _img) : Inpainter(_b, _img) {}

	Color pixelColor(CoordPtr _c);
};

class WeightedInpainter : public Inpainter
{
public:
	WeightedInpainter(BoundaryColors const * const _b, float _pow = 2.f) : Inpainter(_b), m_pow(_pow) {}

	Color pixelColor(CoordPtr _c);

private:
	float m_pow;
};


typedef std::tuple<CoordPtr, float, float, float> CoordFFF;
typedef std::vector<CoordFFF> BoundaryGrads;

class GradientWeightedInpainter : public Inpainter
{
public:
	GradientWeightedInpainter(BoundaryColors const * const _b, CImg<uchar> const * const _img, CImg<uchar> const * const _mask, float _pow = 2.f, float _jitter = 0.35, int _dilation = 2) : Inpainter(_b), m_pow(_pow), m_jitter(100 * _jitter) { init(_img, _mask); }

	Color pixelColor(CoordPtr _c);

private:
	void init(CImg<uchar> const * const _img, CImg<uchar> const * const _mask);

	std::vector<std::vector<int> > m_boundarySides;
	BoundaryGrads m_boundaryGrad;
	float m_maxGrad;

	float m_pow;
	int m_jitter;
	CImg<float> m_attenuation;
};

} // end namespace Diminer

