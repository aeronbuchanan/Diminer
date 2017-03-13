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
#include "boundaryRegions.h"

namespace Diminer
{

// TODO: Inpainter factory
enum Inpainters {BLEED, BLEND, GRADS};

template <class T>
class SpatialGridSquare
{
public:
	SpatialGridSquare() : xmax(0), xmin(0), ymax(0), ymin(0) { contents = std::make_shared<T>(); };
	SpatialGridSquare(double _xn, double _yn, double _xx, double _yx) : xmax(_xx), xmin(_xn), ymax(_yx), ymin(_yn) { contents = std::make_shared<T>(); };
	SpatialGridSquare(double _xn, double _yn, double _xx, double _yx, T* _c) : xmax(_xx), xmin(_xn), ymax(_yx), ymin(_yn), contents(_c) {};

	double sqrdDistanceTo(double x, double y)
	{
		double dx = 0;
		double dy = 0; 

		if ( x < xmin ) dx = xmin - x;
		else if ( x >= xmax ) dx = x - xmax + 1;

		if ( y < ymin ) dy = ymin - y;
		else if ( y > ymax ) dy = y - ymax + 1;

		return dx * dx + dy * dy;
	}

	double hasInside(double x, double y)
	{
		return x >= xmin && x < xmax && y >= ymin && y < ymax;
	}

	// min pixels are included, max pixels are not
	double xmax, xmin, ymax, ymin;
	std::shared_ptr<T> contents;
};

template <class T>
using SpatialGrid = std::vector<SpatialGridSquare<T>>;

typedef std::pair<BoundaryColors::const_iterator, BoundaryColors::const_iterator> BoundaryColorRange;
typedef std::vector<BoundaryColorRange> BoundaryColorRanges;

typedef SpatialGridSquare<BoundaryColorRanges> BoundaryColorGridSq;
typedef SpatialGrid<BoundaryColorRanges> BoundaryColorGrid;

class Inpainter
{
public:
	Inpainter(BoundaryColors const * const);
	~Inpainter() {}

	virtual Color pixelColor(CoordPtr const &);

protected:
	BoundaryColorRanges withinRange(CoordPtr const &);

	BoundaryColors const * m_boundary;
	BoundaryColorGrid m_grid;
	double m_gridSize;
};

class BleedInpainter : public Inpainter
{
public:
	BleedInpainter(BoundaryColors const * const _b) : Inpainter(_b) {}

	Color pixelColor(CoordPtr const & _c);
};

class WeightedInpainter : public Inpainter
{
public:
	WeightedInpainter(BoundaryColors const * const _b, float _pow = 2.f) : Inpainter(_b), m_pow(_pow) {}

	Color pixelColor(CoordPtr const & _c);

private:
	float m_pow;
};


typedef std::tuple<CoordPtr, float, float, float> CoordFFF;
typedef std::vector<CoordFFF> BoundaryGrads;

class GradientWeightedInpainter : public Inpainter
{
public:
	GradientWeightedInpainter(
		BoundaryColors const * const _b, 
		SourceImage const * const _img, 
		MaskImage const * const _mask, 
		MaskImage const * const _regions,
		GradImage const * const _grads,
		int _regionID,
		float _pow = 5.f, 
		float _jitter = 0.35, 
		int _dilation = 2) 
	: Inpainter(_b), m_gradImg(_grads), m_pow(_pow), m_jitter(100 * _jitter) 
	{ 
		init(_img, _mask, _regions, _grads,  _regionID);
	}

	Color pixelColor(CoordPtr const & _c);

private:
	void init(SourceImage const * const _img, MaskImage const * const _mask, MaskImage const * const _regions, GradImage const * const _grads, int _regionID);
	int side(CoordPtr const &, CoordPtr const &, GradImage const * const);

	BoundaryRegions m_boundaryRegions;
	SalientPointRanges m_maxGradPoints;
	GradImage const * const m_gradImg;

	float m_pow;
	int m_jitter;
};

} // end namespace Diminer

