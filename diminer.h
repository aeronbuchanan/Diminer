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

// TODO: clean up include hierarchy 

#pragma once

#include <vector>
#include <utility>

#include "vecn.h"
#include "CImg.h"
using namespace cimg_library;

#define IS_NOT_MASKED 0
#define IS_MASKED 255

typedef unsigned char uchar;
typedef unsigned int uint;

namespace Diminer
{

class Color
{
public:
   Color() : r(0), g(0), b(0), a(255) {}

   uint r, g, b, a;
};

// TODO: img mask test should be switchable
bool imgMaskTest(Color c);

typedef VecN<int, 2> Coord;
typedef std::vector<Coord> Coords;
typedef std::pair<Coord, Color> CoordCol;
typedef std::vector<CoordCol> BoundaryColors;

} // end namespace Diminer



