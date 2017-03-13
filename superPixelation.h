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

#include <stdlib.h>
#include <memory>

#include "diminer.h"

namespace Diminer
{

struct KPoint
{ 
	KPoint() : c1(0), c2(0), c3(0), x(0), y(0) {}; 
	KPoint(double a1, double a2, double a3, double a4, double a5) : c1(a1), c2(a2), c3(a3), x(a4), y(a5) {}; 
	double c1, c2, c3, x, y; 

	static double KPointDistance(KPoint const & a, KPoint const & b, double spatialBias)
	{
		double d1 = a.c1 - b.c1;
		double d2 = a.c2 - b.c2;
		double d3 = a.c3 - b.c3;
		double cd_sqrd = d1*d1 + d2*d2 + d3*d3;
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		double dd_sqrd = dx*dx + dy*dy;
		return sqrt( cd_sqrd + dd_sqrd * spatialBias);
	};
}

typedef std::shared_ptr<KPoint> KPointPtr;
typedef std::vector<KPointPtr> KPoints;


class KMeans
{
public:
	KMeans(KPoints const * const _cs) : centres(_cs) {};

	void iterate(int numSteps);

private:
	

	KPoints const * const centres;
}


		double f = 100 / (spacing * spacing);

class SuperPixelator
{
public:
	SuperPixelator(SourceImage const * const _img) : m_img(_img) {};

private:
	void tweekCentre(KPointPtr const &);

	SourceImage const * const m_img;
}


} // end namespace Diminer




