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

#include <vector>
#include <tuple>

#include "CImg.h"
using namespace cimg_library;

// TODO: Diminer namespace

typedef unsigned char uchar;
typedef unsigned int uint;

typedef std::pair<int, int> Coord;
typedef std::vector<Coord> Coords;

struct Color { uchar r; uchar g; uchar b; };

typedef std::pair<Coord, float> CoordF;
typedef std::vector<CoordF> BoundaryInfo;

typedef std::pair<Coord, Color> CoordCol;
typedef std::vector<CoordCol> BoundaryColors;

typedef std::tuple<Coord, float, float, float> CoordFFF; // FFF = {mag, nx, ny}
typedef std::vector<CoordFFF> BoundaryGrads;


class AnimDisp
{
public:
	AnimDisp(CImgDisplay* const _disp) : m_disp(_disp) {}

	void update(CImg<uchar> * const _img) { m_disp->operator=(*_img); m_disp->wait(); }

	CImgDisplay* m_disp;
};

static bool maskTest(uchar _v) { return _v > 127; } // Humph
static bool imgMaskTest(Color _c) { return _c.r <= 1 && _c.g <= 1 && _c.b >= 254; }
//static bool imgMaskTest(Color _c) { return _c.r >= 254 && _c.g >= 87 && _c.g <= 98 && _c.b <= 1; }

class FillHelper
{
public:
	FillHelper(CImg<uchar> const * const _src, CImg<uchar> * const _tgt) : m_src(_src), m_tgt(_tgt) {}

	void fill( Coord _xy, uchar col)
	{
		m_stack.clear();

		// start fill
		m_stack.push_back(_xy);
		while ( !m_stack.empty() )
		{
			Coord foc = m_stack.back();
			m_stack.pop_back();
			int fx = foc.first;
			int fy = foc.second;

			// mark
			(*m_tgt)(fx,fy) = col;

			// scan
			bool ffup = lookUp(fx, fy, false);
			bool ffwn = lookDown(fx, fy, false);

			// scan back
			bool fup = ffup;
			bool fwn = ffwn;
			int nx = fx;
			while ( --nx > 0 && maskTest((*m_src)(nx,fy)) && (*m_tgt)(nx,fy) == 0 )
			{
				(*m_tgt)(nx,fy) = col;

				fup = lookUp(nx, fy, fup);
				fwn = lookDown(nx, fy, fwn);
			}

			// scan forward
			fup = ffup;
			fwn = ffwn;
			nx = fx;
			while ( ++nx < m_src->width() && maskTest((*m_src)(nx,fy)) && (*m_tgt)(nx,fy) == 0 )
			{
				(*m_tgt)(nx,fy) = col;

				fup = lookUp(nx, fy, fup);
				fwn = lookDown(nx, fy, fwn);
			}
		}
	}

private:
	bool lookUp(int _x, int _y, bool _b) { if ( _y > 0 ) _b = inspect(_x, --_y, _b); return _b;	}
	bool lookDown(int _x, int _y, bool _b) { if ( _y < m_src->height() - 1 ) _b = inspect(_x, ++_y, _b); return _b; }

	bool inspect(int _x, int _y, bool _b)
	{
		if ( maskTest((*m_src)(_x, _y)) && (*m_tgt)(_x, _y) == 0 )
		{
			if ( !_b )
			{
				_b = true;
				m_stack.push_back(Coord(_x, _y));
			}
		}
		else
		{
			_b = false;
		}
		return _b;
	}

	CImg<uchar> const * m_src;
	CImg<uchar> * m_tgt;
	Coords m_stack;
};
