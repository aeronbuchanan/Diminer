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
#include <memory>

#include "vecn.h"
#include "CImg.h"
using namespace cimg_library;

#define IS_NOT_MASKED 0
#define IS_MASKED 255

typedef unsigned char uchar;

// TODO: un-lazy these structs...

namespace Diminer
{

class Color
{
public:
	Color() : r(0), g(0), b(0), a(255) {}
	Color(uint _c) : r(_c), g(_c), b(_c), a(255) {}
	Color(uint _c, uint _a) : r(_c), g(_c), b(_c), a(_a) {}
	Color(uint _r, uint _g, uint _b) : r(_r), g(_g), b(_b), a(255) {}
	Color(uint _r, uint _g, uint _b, uint _a) : r(_r), g(_g), b(_b), a(_a) {}

	uint r, g, b, a;
};

class CoordBase
{
public:
	CoordBase() : x(0), y(0) {};
	CoordBase(int _x, int _y) : x(_x), y(_y) {};
	int x, y;
};

// TODO: proper inheritance!
class Coord : public CoordBase
{
public:
	Coord() : col(Color()) {};
	Coord(int _x, int _y) : CoordBase(_x, _y), col(Color()) {}
	Coord(int _x, int _y, Color _col) : CoordBase(_x, _y), col(_col) {}
	Color col;
};

typedef std::shared_ptr<Coord> CoordPtr;

typedef std::vector<CoordPtr> Coords;
typedef std::vector<CoordPtr> BoundaryColors;

// TODO: img mask test should be switchable
bool imgMaskTest(Color c);


} // end namespace Diminer



