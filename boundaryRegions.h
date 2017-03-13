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

struct BoundaryRegion
{
	CImg<float> colData;
	int xmin, xmax, ymin, ymax;
};


typedef std::shared_ptr<BoundaryRegion> BoundaryRegionPtr;
typedef std::vector<BoundaryRegionPtr> BoundaryRegions;
// TODO: provide macro for creating shared_ptrs

typedef Coords::const_iterator SalientPointItr;
typedef std::vector<SalientPointItr> SalientPointRanges;

double findSalientPoints(Coords const * const boundaryCoords, GradImage const * const gradImag, SalientPointRanges & output);

} // end namespace Diminer




