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
#include "boundaryChains.h"

using namespace Diminer;

Inpainter::Inpainter(BoundaryColors const * const _b) : m_boundary(_b), m_gridSize(200)
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
	int count = 0;
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

// Determine on which side point 'c' is of line of gradient at boundary point 'b'
int GradientWeightedInpainter::side(CoordPtr const & c, CoordPtr const & b, GradImage const * const gradImg)
{
	double nx = (*gradImg)(b->x, b->y, 0, 1);
	double ny = (*gradImg)(b->x, b->y, 0, 2);
	double dx = double(c->x) - double(b->x);
	double dy = double(c->y) - double(b->y);
	double dot = dx * nx + dy * ny;
	return dot < 0 ? -1 : dot > 0 ? 1 : 0;
}

void GradientWeightedInpainter::init(SourceImage const * const _img, MaskImage const * const _mask, MaskImage const * const _regions, GradImage const * const _grads, int regionID)
{
	if ( m_boundary->size() < 3 ) return; // not worth it

	// just for typing convenience
	SourceImage img(*_img);
	MaskImage mask(*_mask);
	MaskImage regions(*_regions);
	GradImage gradImg(*_grads);

	int W = img.width();
	int H = img.height();

	// pow series look-up
	int extent = ceil( log(1e8) / log(m_pow) ); // find out how far we need to go...
	CImg<float> attenuation(extent, extent, 1, 1, 1.0);
	cimg_forXY(attenuation,x,y)
		if ( x > 0 || y > 0 ) attenuation(x,y) = 1.f / pow(x * x + y * y, m_pow); // closer pixels are more heavily weighted

	// non-max suppress
	// TODO: check for and deal with looped boundaries
	// TODO: measure "peakiness" and also supress insufficiently peaky maximums
	double maxGrad = 0.f;
	m_maxGradPoints.push_back(m_boundary->begin());
	auto prev = m_boundary->begin();
	auto curr = prev + 1;
	auto next = curr + 1;
	for ( ; next != m_boundary->end(); prev++, curr++, next++ )
	{
		int x = (*curr)->x;
		int y = (*curr)->y;
		double mag = gradImg(x, y, 0, 0);
		if ( mag >= gradImg((*prev)->x, (*prev)->y, 0, 0) && mag > gradImg((*next)->x, (*next)->y, 0, 0) )
		{
			m_maxGradPoints.push_back( curr );
			if ( mag > maxGrad ) maxGrad = mag;
		}
	}
	if ( m_boundary->size() > 1 )
		m_maxGradPoints.push_back(m_boundary->end() - 1);

	std::cout << "Found " << m_maxGradPoints.size() << " worthy gradient points" << std::endl;

	if ( m_maxGradPoints.size() < 3 ) return;

	/* DEBUG
	std::cout << "[" << std::endl;
	for ( auto cmi = cms.begin(); cmi != cms.end(); cmi++ )
	{
		std::cout << "[" << std::endl;
		(*cmi).printChains();
		std::cout << "]," << std::endl;
	}
	std::cout << "[] ]" << std::endl;
	// DEBUG END */

	// store attenuation factors
	for ( auto gi = m_maxGradPoints.begin(); gi != m_maxGradPoints.end(); gi++ )
		gradImg((**gi)->x, (**gi)->y, 0, 0) = pow(1.f - gradImg((**gi)->x, (**gi)->y, 0, 0) / maxGrad, 2.f);

	// create regions for each "superpixel" boundary segment
	auto gCurr = m_maxGradPoints.begin();
	auto gNext = gCurr + 1;
	//int count = 0; // DEBUG
	// TODO: shift secnd to first
	for ( ; gNext != m_maxGradPoints.end(); gCurr++, gNext++ )
	{
		CoordPtr first = **gCurr;
		CoordPtr secnd = **gNext;

		double first_mag = gradImg(first->x, first->y, 0, 0);
		double secnd_mag = gradImg(secnd->x, secnd->y, 0, 0);

		CoordPtr inside = std::make_shared<Coord>( (first->x + secnd->x) / 2.0, (first->y + secnd->y) / 2.0);
		int inside_first = side(inside, first, _grads);
		int inside_secnd = side(inside, secnd, _grads);

		std::shared_ptr<BoundaryRegion> patch = std::make_shared<BoundaryRegion>();

		patch->xmin = std::max(0, std::min(int(first->x), int(secnd->x)) - extent);
		patch->xmax = std::min(W - 1, std::max(int(first->x), int(secnd->x)) + extent);
		patch->ymin = std::max(0, std::min(int(first->y), int(secnd->y)) - extent);
		patch->ymax = std::min(H - 1, std::max(int(first->y), int(secnd->y)) + extent);
		patch->colData.assign(patch->xmax - patch->xmin + 1, patch->ymax - patch->ymin + 1, 1, 4, 0);

		//std::cout << count << ": Boundary points (" << first->x << ", " << first->y << ") & (" << secnd->x << ", " << secnd->y << ") ==> box (" << patch->xmin << ", " << patch->ymin << ")--(" << patch->xmax << ", " << patch->ymax << ")" << std::endl; // DEBUG
		//std::cout << "Inside point (" << inside->x << ", " << inside->y << ") sidedness: " << inside_first << "; " << inside_secnd << std::endl;

		cimg_forXY(patch->colData,xx,yy)
		{
			int x = patch->xmin + xx;
			int y = patch->ymin + yy;

			//if ( xx == 0 && yy != 0 ) std::cout << std::endl;
			//std::cout << int(regions(x,y)) << ", ";

			if ( mask(x,y) && regions(x,y) == regionID )
			{
				CoordPtr cc = std::make_shared<Coord>(x, y);

				double r = 0;
				double g = 0;
				double b = 0;
				double k = 0;

				for ( auto bi = *gCurr; bi != *gNext; bi++ )
				{
					CoordPtr bc = (*bi);

					int dx = abs(double(x) - double(bc->x));
					int dy = abs(double(y) - double(bc->y));

					double ki = 0;

					// closer pixels more heavily weighted
					if ( dx < extent && dy < extent ) ki = attenuation(dx, dy); // TODO: check look-up is faster than pow()
					else ki = attenuation(extent - 1, extent - 1);

					// same side as gradient lines?
					if ( side(cc, first, _grads) != inside_first ) ki *= first_mag;
					if ( side(cc, secnd, _grads) != inside_secnd ) ki *= secnd_mag;

					Color c = bc->col;
					r += ki * double(c.r);
					g += ki * double(c.g);
					b += ki * double(c.b);
					k += ki;
				}

				patch->colData(xx, yy, 0, 0) = r;
				patch->colData(xx, yy, 0, 1) = g;
				patch->colData(xx, yy, 0, 2) = b;
				patch->colData(xx, yy, 0, 3) = k;
			}
		}

		m_boundaryRegions.push_back(patch);

		/* DEBUG
		count++;
		char * name = (char*)malloc(250);
		sprintf(name, "debug_regions_%03d.png", count);
		patch->colData.save(name);
		free(name);
		// DEBUG END */
	}

	/* note which side each boundary coord is on
	for ( auto gc = m_maxGradPoints.begin(); gc != m_maxGradPoints.end(); gc++ )
	{
		CoordPtr coord = *gc;
		double mag = gradImg(coord->x, coord->y, 0, 0);
		double nx = gradImg(coord->x, coord->y, 0, 1);
		double ny = gradImg(coord->x, coord->y, 0, 2);

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

Color GradientWeightedInpainter::pixelColor(CoordPtr const & _c)
{
	if ( m_boundary->size() < 3 ) return Color(); // don;t bother

	double r = 0;
	double g = 0;
	double b = 0;
	double k = 0;

	//std::cout << "(" << _c->x << ", " << _c->y << "):" << std::endl;

	for ( auto pi = m_boundaryRegions.begin(); pi != m_boundaryRegions.end(); pi++ )
	{
		auto patch = (*pi);
		if ( _c->x >= patch->xmin && _c->x <= patch->xmax && _c->y >= patch->ymin && _c->y <= patch->ymax )
		{
			int xx = _c->x - patch->xmin;
			int yy = _c->y - patch->ymin;
			r += patch->colData(xx, yy, 0, 0);
			g += patch->colData(xx, yy, 0, 1);
			b += patch->colData(xx, yy, 0, 2);
			k += patch->colData(xx, yy, 0, 3);

			//std::cout << "+ [" << patch->colData(xx, yy, 0, 0) << ", " << patch->colData(xx, yy, 0, 1) << ", " <<patch->colData(xx, yy, 0, 2) << ", " << patch->colData(xx, yy, 0, 3) << "] = [" << r << ", " << g << ", " << b << ", " << k << "] = [" << int(uchar( r / k )) << "; " << int(uchar( g / k )) << "; " << int(uchar( b / k )) << "]" << std::endl;
		}
	}

	return Color( uchar(r / k), uchar(g / k), uchar(g / k) );
}

