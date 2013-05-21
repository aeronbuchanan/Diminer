/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of Diminer, an digital inpainting resource.
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

#include "common.h"

class Inpainter
{
public:
	Inpainter(BoundaryColors const * const _b) : m_boundary(_b) {}
	~Inpainter() {}

	virtual Color pixelColor(Coord) { return Color(); }

protected:
	BoundaryColors const * m_boundary;
};

class BleedInpainter : public Inpainter
{
public:
	BleedInpainter(BoundaryColors const * const _b) : Inpainter(_b) {}

	Color pixelColor(Coord _c);
};

class WeightedInpainter : public Inpainter
{
public:
	WeightedInpainter(BoundaryColors const * const _b, float _pow = 2.f) : Inpainter(_b), m_pow(_pow) {}

	Color pixelColor(Coord _c);

private:
	float m_pow;
};

class GradientWeightedInpainter : public Inpainter
{
public:
	GradientWeightedInpainter(BoundaryColors const * const _b, CImg<uchar> const * const _img, float _pow = 2.f) : Inpainter(_b), m_pow(_pow) { init(_img); }
	GradientWeightedInpainter(BoundaryColors const * const _b, CImg<uchar> const * const _img, AnimDisp * const _a) : Inpainter(_b), debugDisp(_a), m_pow(2.f) { init(_img); }

	Color pixelColor(Coord _c);

private:
	void init(CImg<uchar> const * const _img);

	std::vector<std::vector<int> > m_boundarySides;
	BoundaryGrads m_boundaryGrad;
	float m_pow;

	AnimDisp * debugDisp; // DEBUG

};

