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

#include <iostream>
#include <algorithm>
#include <limits>

#include "inpainters.h"

using namespace Diminer;

Inpainter::Inpainter(BoundaryColors const * const _b) : m_boundary(_b), m_gridSize(100)
{
	std::cout << "Gridding... ";
	std::flush(std::cout); // TODO: flush

	// populate spatial grid
	double xstart = std::numeric_limits<double>::max();
	double xend = 0;
	double ystart = std::numeric_limits<double>::max();
	double yend = 0;
	for ( auto ci = _b->begin(); ci != _b->end(); ci++ )
	{
		double x = double((*ci)->x);
		double y = double((*ci)->y);
		if ( x >= xend ) xend = x + 1;
		if ( x <= xstart ) xstart = x - 1;
		if ( y >= yend ) yend = y + 1;
		if ( y <= ystart ) ystart = y - 1;
	}

	for ( double x = xstart; x < xend; x += m_gridSize )
		for ( double y = ystart; y < yend; y += m_gridSize )
			m_grid.push_back(BoundaryColorGridSq(x, y, x + m_gridSize, y + m_gridSize));

	if ( ! m_grid.size() ) return;

	// find grid sq for first boundary point
	BoundaryColorGridSq * curr = &m_grid[0];
	for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
		if ( sq->hasInside((*_b)[0]->x, (*_b)[0]->y) ) { curr = &(*sq); break; }

	auto start = _b->begin();

	for ( auto ci = _b->begin(); ci != _b->end(); ci++ )
	{
		CoordPtr c = (*ci);

		if ( ! curr->hasInside(c->x, c->y) )
		{
			curr->contents->push_back(BoundaryColorRange(start, ci));
			start = ci;
			for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
				if ( sq->hasInside(c->x, c->y) ) { curr = &(*sq); break; }
		}
	}

	curr->contents->push_back(BoundaryColorRange(start, _b->end()));

	// remove empty grid squares
	while ( m_grid.begin()->contents->size() == 0 )
		m_grid.erase(m_grid.begin());
	auto last = m_grid.begin();
	for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
	{
		if ( sq->contents->size() == 0 )
		{
			m_grid.erase(sq);
			sq = last;
		}
		else
		{
			last = sq;
		}
	}
	
	std::cout << "complete." << std::endl;	

	// TODO: check full coverage
	/* DEBUG
	uint count = 0;
	for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
	{
		std::cout << "GridSq (" << sq->xmin << ", " << sq->ymin << ")--(" << sq->xmax << ", " << sq->ymax << ") :" << std::endl;
		for ( auto pi = sq->contents->begin(); pi != sq->contents->end(); pi++ )
			for ( auto ci = pi->first; ci != pi->second; ci++ )
				std::cout << ++count << "(" << (*ci)->x << ", " << (*ci)->y << ") ";
		std::cout << std::endl;
	}
	if ( count != m_boundary->size() )
		std::cout << "ERROR: wrong number of coords in grid - expecting " << m_boundary->size() << ", but found " << count << std::endl;
	// DEBUG END */
}

typedef std::pair<double, std::shared_ptr<BoundaryColorRanges> > dBCR;
bool compare_dBCRs(dBCR const & a, dBCR const & b) { return (a.first < b.first); }

// TODO: better use of pass by reference
BoundaryColorRanges Inpainter::withinRange(CoordPtr const & _c)
{
	BoundaryColorRanges rs;

	if ( m_grid.size() >= 9)
	{
		std::vector<dBCR> v;
		for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
			v.push_back(dBCR(sq->sqrdDistanceTo(_c->x, _c->y), sq->contents));

		// TODO: pre-bake 8-neighbour connectivity so only nearest quadrant need be found
		std::sort(v.begin(), v.end(), compare_dBCRs);

		for ( auto ri = v.begin(); ri != v.end() && ri < v.begin() + 9; ri++)
			rs.insert(rs.end(), ri->second->begin(), ri->second->end());
	}
	else
	{
		for ( auto sq = m_grid.begin(); sq != m_grid.end(); sq++ )
			rs.insert(rs.end(), sq->contents->begin(), sq->contents->end());
	}

	return rs;
}

Color Inpainter::pixelColor(CoordPtr const &)
{
	return Color();
}

Color BleedInpainter::pixelColor(CoordPtr const & _c)
{
	double bestDist = std::numeric_limits<double>::max();
	Color bestCol = Color();
	double bestCount = 0;

	BoundaryColorRanges rs = withinRange(_c);

	for ( auto ri = rs.begin() ; ri != rs.end(); ri++ )
	{
		for ( auto bc = ri->first; bc != ri->second; bc++ )
		{
			CoordPtr coord = (*bc);
			Color color = coord->col;
			double dx = double(_c->x) - double(coord->x);
			double dy = double(_c->y) - double(coord->y);
			double dist = dx * dx + dy * dy;

			if ( dist < bestDist )
			{
				bestDist = dist;
				bestCol = color;
				bestCount = 1;
			}
			else if ( dist == bestDist )
			{
				double r = bestCol.r * bestCount;
				double g = bestCol.g * bestCount;
				double b = bestCol.b * bestCount;
				++bestCount;
				bestCol.r = (r + color.r) / bestCount;
				bestCol.g = (g + color.g) / bestCount;
				bestCol.b = (b + color.b) / bestCount;
			}
		}
	}
	
	/* DEBUG
	double DEBUGbestDist = std::numeric_limits<double>::max();
	Coord DEBUG2 = Coord();
	for ( auto bc = m_boundary->begin(); bc != m_boundary->end(); bc++ )
	{
			CoordPtr coord = (*bc);
			double dx = double(_c->x) - double(coord->x);
			double dy = double(_c->y) - double(coord->y);
			double dist = dx * dx + dy * dy;

			if ( dist < DEBUGbestDist )
			{
				DEBUGbestDist = dist;
				DEBUG2 = (*coord);
			}
	}
	if ( DEBUGbestDist != bestDist && (DEBUG1.x != DEBUG2.x || DEBUG1.y != DEBUG2.y) )
		std::cout << "ERROR! at (" << _c->x << ", " << _c->y << ") expected (" << DEBUG2.x << ", " << DEBUG2.y << ")@" << DEBUGbestDist << ", but got (" << DEBUG1.x << ", " << DEBUG1.y << ")@" << bestDist << std::endl;
	// DEBUG END */

	return bestCol;
}

Color WeightedInpainter::pixelColor(CoordPtr const & _c)
{
	double weightSum = 0;
	double redSum = 0;
	double greSum = 0;
	double bluSum = 0;

	BoundaryColorRanges rs = withinRange(_c);

	for ( auto ri = rs.begin() ; ri != rs.end(); ri++ )
	{
		for ( auto bc = ri->first; bc != ri->second; bc++ )
		{
			CoordPtr coord = (*bc);
			Color color = coord->col;
			double dx = double(_c->x) - double(coord->x);
			double dy = double(_c->y) - double(coord->y);
			double distSqrd = dx * dx + dy * dy;
			double weight = 1.f / pow(distSqrd, m_pow); // closer pixels more heavily weighted
			// TODO: check pow() is faster than pre-processed look-up

			redSum += double(color.r) * weight;
			greSum += double(color.g) * weight;
			bluSum += double(color.b) * weight;
			weightSum += weight;
		}
	}

	Color c;
	c.r = uchar(redSum / weightSum);
	c.g = uchar(greSum / weightSum);
	c.b = uchar(bluSum / weightSum);

	return c;
}

void GradientWeightedInpainter::init(CImg<uchar> const * const _img, CImg<uchar> const * const _mask)
{
	if ( m_boundary->size() < 3 ) return; // not worth it

	CImg<uchar> img(*_img);
	CImg<uchar> mask(*_mask);

	// greyscale
	CImg<double> g(img.width(), img.height(), 1, 1, 0.f);
	cimg_forXYC(img, x,y,k)
		g(x,y) += img(x,y,k);

	// gradient
	m_gradImg.assign(img.width(), img.height(), 3, 1, 0.f);
	//CImg<uchar> debugGrad(img.width(), img.height(), 1, 3, 0);

	// pow series look-up
	int d = ceil( log(1e8) / log(m_pow) ); // find out how far we need to go...
	m_attenuation.assign(d, d, 1, 1, 1.0);
	cimg_forXY(m_attenuation,x,y)
		if ( x > 0 || y > 0 ) m_attenuation(x,y) = 1.f / pow(x * x + y * y, m_pow); // closer pixels are more heavily weighted

	// DEBUG
	cimg_forXY(m_attenuation,x,y)
	{
		if ( x == 0 ) std::cout << "\n";
		std::cout << m_attenuation(x,y) << ", ";
	}
	std::cout << std::endl;

/*
	@misc{NRIGO,
	Author = {Pavel Holoborodko},
	Title = {Noise Robust Gradient Operators.},
	howpublished = {\url{http://www.holoborodko.com/pavel/image-processing/edge-detection/}}
	year = {2009}
	}

	dI/dx = 1/32 [ -1 -2 0 2 1; -2 -4 0 4 2; -1 -2 0 2 1 ]

	(like an extended Sobel kernel)
*/

	// TODO: only calculate for boundary
	// TODO: optimize recalculations
	for ( int y = 0; y < g.height(); ++y )
	{
		int y_pp = std::max( 0, y - 2);
		int y_p = std::max( 0, y - 1);
		int y_n = std::min( g.height() - 1, y + 1);
		int y_nn = std::min( g.height() - 1, y + 2);

		for ( int x = 0; x < g.width(); ++x )
		{
			int x_pp = std::max( 0, x - 2);
			int x_p = std::max( 0, x - 1);
			int x_n = std::min( g.width() - 1, x + 1);
			int x_nn = std::min( g.width() - 1, x + 2);

			double gxppyp = mask(x_pp, y_p) ? g(x,y) : g(x_pp, y_p);
			double gxpyp = mask(x_p, y_p) ? g(x,y) : g(x_p, y_p);
			double gxnyp = mask(x_n, y_p) ? g(x,y) : g(x_n, y_p);
			double gxnnyp = mask(x_nn, y_p) ? g(x,y) : g(x_nn, y_p);

			double gxppy = mask(x_pp, y) ? g(x,y) : g(x_pp, y);
			double gxpy = mask(x_p, y) ? g(x,y) : g(x_p, y);
			double gxny = mask(x_n, y) ? g(x,y) : g(x_n, y);
			double gxnny = mask(x_nn, y) ? g(x,y) : g(x_nn, y);

			double gxppyn = mask(x_pp, y_n) ? g(x,y) : g(x_pp, y_n);
			double gxpyn = mask(x_p, y_n) ? g(x,y) : g(x_p, y_n);
			double gxnyn = mask(x_n, y_n) ? g(x,y) : g(x_n, y_n);
			double gxnnyn = mask(x_nn, y_n) ? g(x,y) : g(x_nn, y_n);

			double gx = ( gxnnyn + gxnnyp + 2 * (gxnyn + gxnny + gxnyp) + 4 * gxny ) - ( gxppyp + gxppyn + 2 * (gxppy + gxpyp + gxpyn) + 4 * gxpy );
			gx /= 32;

			double gxpypp = mask(x_p, y_pp) ? g(x,y) : g(x_p, y_pp, 0);
			//double gxpyp;
			//double gxpyn;
			double gxpynn = mask(x_p, y_nn) ? g(x,y) : g(x_p, y_nn, 0);

			double gxypp = mask(x, y_pp) ? g(x,y) : g(x, y_pp, 0);
			double gxyp = mask(x, y_p) ? g(x,y) : g(x, y_p, 0);
			double gxyn = mask(x, y_n) ? g(x,y) : g(x, y_n, 0);
			double gxynn = mask(x, y_nn) ? g(x,y) : g(x, y_nn, 0);

			double gxnypp = mask(x_n, y_pp) ? g(x,y) : g(x_n, y_pp, 0);
			//double gxnyp;
			//double gxnyn;
			double gxnynn = mask(x_n, y_nn) ? g(x,y) : g(x_n, y_nn, 0);

			double gy = ( gxnynn + gxpynn + 2 * (gxnyn + gxynn + gxpyn) + 4 * gxyn ) - ( gxpypp + gxnypp + 2 * (gxypp + gxpyp + gxnyp) + 4 * gxyp );
			gy /= 32;

			double mag = sqrt( gx * gx + gy * gy );
			m_gradImg(x,y,0) = mag;
			m_gradImg(x,y,1) = gx / mag;
			m_gradImg(x,y,2) = gy / mag;

			//debugGrad(x,y,0,0) = mag;
			//debugGrad(x,y,0,1) = abs(gx);
			//debugGrad(x,y,0,2) = abs(gy);
		}
	}

	//debugGrad.save("gradient.bmp");

	// non-max suppress
	// TODO: check for and deal with looped boundaries
	// TODO: measure "peakiness" and also supress insufficiently peaky maximums
	double maxGrad = 0.f;
	auto prev = m_boundary->begin();
	auto curr = prev + 1;
	auto next = curr + 1;
	for ( ; next != m_boundary->end(); prev++, curr++, next++ )
	{
		uint x = (*curr)->x;
		uint y = (*curr)->y;
		double mag = m_gradImg(x, y, 0);
		if ( mag >= m_gradImg((*prev)->x, (*prev)->y, 0) && mag > m_gradImg((*next)->x, (*next)->y, 0) )
		{
			m_maxGradPoints.push_back( *curr );
			if ( mag > maxGrad ) maxGrad = mag;
		}
	}

	std::cout << "Found " << m_maxGradPoints.size() << " worthy gradient points" << std::endl;

	for ( auto gi = m_maxGradPoints.begin(); gi != m_maxGradPoints.end(); gi++ )
		m_gradImg((*gi)->x, (*gi)->y, 0) = pow(1.f - m_gradImg((*gi)->x, (*gi)->y, 0) / maxGrad, 2.f);

	/* note which side each boundary coord is on
	for ( auto gc = m_maxGradPoints.begin(); gc != m_maxGradPoints.end(); gc++ )
	{
		CoordPtr coord = *gc;
		double mag = m_gradImg(coord->x, coord->y, 0);
		double nx = m_gradImg(coord->x, coord->y, 1);
		double ny = m_gradImg(coord->x, coord->y, 2);

		std::vector<int> info(m_boundary->size());
		for ( auto gc = m_boundary->begin(); gc != m_boundary->end(); gc++ )
		{
			CoordPtr c = *gc;
			double dx = double(c->x) - double(coord->x);
			double dy = double(c->y) - double(coord->y);
			double dot = dx * nx + dy * ny;
			info[j] = dot < 0 ? -1 : dot > 0 ? 1 : 0; // function
		}
		m_boundarySides.push_back(info);
	}
	*/
}

// Determine on which side point 'c' is of line of gradient at boundary point 'b'
int GradientWeightedInpainter::dot(CoordPtr const & c, CoordPtr const & b)
{
	double nx = m_gradImg(b->x, b->y, 1);
	double ny = m_gradImg(b->x, b->y, 2);
	double dx = double(c->x) - double(b->x);
	double dy = double(c->y) - double(b->y);
	double dot = dx * nx + dy * ny;
	return dot < 0 ? -1 : dot > 0 ? 1 : 0;
}

Color GradientWeightedInpainter::pixelColor(CoordPtr const & _c)
{
	if ( m_boundary->size() < 3 ) return Color(); // don;t bother

	// determine which side of all the major grad 
	std::vector<int> sides(m_maxGradPoints.size());
	for ( uint i = 0 ; i < m_maxGradPoints.size(); ++i )
		sides[i] = dot(_c, m_maxGradPoints[i]);

	uint maxd = m_attenuation.width();
	BoundaryColorRanges rs = withinRange(_c);

	double r = 0;
	double g = 0;
	double b = 0;
	double k = 0;

	// distance weightings, attenuated by gradient segments
	for ( auto ri = rs.begin() ; ri != rs.end(); ri++ )
	{
		for ( auto bc = ri->first; bc != ri->second; bc++ )
		{
			if ( rand() % 100 >= m_jitter )
			{
				CoordPtr coord = (*bc);

				uint dx = abs(double(_c->x) - double(coord->x));
				uint dy = abs(double(_c->y) - double(coord->y));

				double ki = 0;

				// closer pixels more heavily weighted
				if ( dx < maxd && dy < maxd ) ki = m_attenuation(dx, dy); // TODO: check look-up is faster than pow()
				else ki = m_attenuation(maxd - 1, maxd - 1);

				// influence of pixels on different side of gradient lines attenuated
				for ( uint i = 0 ; i < m_maxGradPoints.size(); ++i ) // TODO: reduce the search space on grad lines
					if ( sides[i] != dot(coord, m_maxGradPoints[i]) ) ki *= m_gradImg(m_maxGradPoints[i]->x, m_maxGradPoints[i]->y, 0) ;

				Color c = coord->col;
				r += ki * double(c.r);
				g += ki * double(c.g);
				b += ki * double(c.b);
				k += ki;
			}
		}
	}
	
	Color c;
	c.r = uchar(r / k);
	c.g = uchar(g / k);
	c.b = uchar(b / k);

	return c;
}

