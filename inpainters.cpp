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
#include <limits>

#include "inpainters.h"
#include "boundaryChains.h"

using namespace Diminer;

Color Inpainter::pixelColor(Coord) { return Color(); }

Color BleedInpainter::pixelColor(Coord _c)
{
	float bestDist = std::numeric_limits<float>::max();
	uchar bestR = 0;
	uchar bestG = 0;
	uchar bestB = 0;
	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		Coord coord = (*m_boundary)[i].first;
		Color color = (*m_boundary)[i].second;
		float dx = _c.x() - coord.x();
		float dy = _c.y() - coord.y();

		float dist = dx * dx + dy * dy;
		if ( dist < bestDist )
		{
			bestDist = dist;
			bestR = color.r;
			bestG = color.g;
			bestB = color.b;
		}
	}

	Color c;
	c.r = bestR;
	c.g = bestG;
	c.b = bestB;

	return c;
}

Color WeightedInpainter::pixelColor(Coord _c)
{
	float weightSum = 0;
	float redSum = 0;
	float greSum = 0;
	float bluSum = 0;

	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		Coord coord = (*m_boundary)[i].first;
		Color color = (*m_boundary)[i].second;
		float dx = _c.x() - coord.x();
		float dy = _c.y() - coord.y();
		float distSqrd = dx * dx + dy * dy;
		float weight = 1.f / pow(distSqrd, m_pow); // closer pixels more heavily weighted

		redSum += float(color.r) * weight;
		greSum += float(color.g) * weight;
		bluSum += float(color.b) * weight;
		weightSum += weight;
	}

	Color c;
	c.r = uchar(redSum / weightSum);
	c.g = uchar(greSum / weightSum);
	c.b = uchar(bluSum / weightSum);

	return c;
}

void GradientWeightedInpainter::init(CImg<uchar> const * const _img, CImg<uchar> const * const _mask)
{
	CImg<uchar> img(*_img);
	CImg<uchar> mask(*_mask);

	// greyscale
	CImg<float> g(img.width(), img.height(), 1, 1, 0.f);
	cimg_forXYC(img, x,y,k)
		g(x,y) += img(x,y,k);

	// gradient
	CImg<float> grad(img.width(), img.height(), 3, 1, 0.f);
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

			float gxppyp = mask(x_pp, y_p) ? g(x,y) : g(x_pp, y_p);
			float gxpyp = mask(x_p, y_p) ? g(x,y) : g(x_p, y_p);
			float gxnyp = mask(x_n, y_p) ? g(x,y) : g(x_n, y_p);
			float gxnnyp = mask(x_nn, y_p) ? g(x,y) : g(x_nn, y_p);

			float gxppy = mask(x_pp, y) ? g(x,y) : g(x_pp, y);
			float gxpy = mask(x_p, y) ? g(x,y) : g(x_p, y);
			float gxny = mask(x_n, y) ? g(x,y) : g(x_n, y);
			float gxnny = mask(x_nn, y) ? g(x,y) : g(x_nn, y);

			float gxppyn = mask(x_pp, y_n) ? g(x,y) : g(x_pp, y_n);
			float gxpyn = mask(x_p, y_n) ? g(x,y) : g(x_p, y_n);
			float gxnyn = mask(x_n, y_n) ? g(x,y) : g(x_n, y_n);
			float gxnnyn = mask(x_nn, y_n) ? g(x,y) : g(x_nn, y_n);

			float gx = ( gxnnyn + gxnnyp + 2 * (gxnyn + gxnny + gxnyp) + 4 * gxny ) - ( gxppyp + gxppyn + 2 * (gxppy + gxpyp + gxpyn) + 4 * gxpy );
			gx /= 32;

			float gxpypp = mask(x_p, y_pp) ? g(x,y) : g(x_p, y_pp, 0);
			//float gxpyp;
			//float gxpyn;
			float gxpynn = mask(x_p, y_nn) ? g(x,y) : g(x_p, y_nn, 0);

			float gxypp = mask(x, y_pp) ? g(x,y) : g(x, y_pp, 0);
			float gxyp = mask(x, y_p) ? g(x,y) : g(x, y_p, 0);
			float gxyn = mask(x, y_n) ? g(x,y) : g(x, y_n, 0);
			float gxynn = mask(x, y_nn) ? g(x,y) : g(x, y_nn, 0);

			float gxnypp = mask(x_n, y_pp) ? g(x,y) : g(x_n, y_pp, 0);
			//float gxnyp;
			//float gxnyn;
			float gxnynn = mask(x_n, y_nn) ? g(x,y) : g(x_n, y_nn, 0);

			float gy = ( gxnynn + gxpynn + 2 * (gxnyn + gxynn + gxpyn) + 4 * gxyn ) - ( gxpypp + gxnypp + 2 * (gxypp + gxpyp + gxnyp) + 4 * gxyp );
			gy /= 32;

			float mag = sqrt( gx * gx + gy * gy );
			grad(x,y,0) = mag;
			grad(x,y,1) = gx / mag;
			grad(x,y,2) = gy / mag;

			//debugGrad(x,y,0,0) = mag;
			//debugGrad(x,y,0,1) = abs(gx);
			//debugGrad(x,y,0,2) = abs(gy);
		}
	}

	//debugGrad.save("gradient.bmp");

	std::cout << "DEBUG: length of boundary = " << m_boundary->size() << std::endl;

	// get chain boundary coords in order
	ChainManager cm;
	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		cm.addCoord((*m_boundary)[i].first);
	}
	if ( ! cm.isGood(img.width(), img.height()) )
	{
		std::cout << "ERROR: boundary ordering failure (fragmented or neither a loop nor spanning the image) - simplify mask and try again?" << std::endl;
	}

	// walk along chainGrouping
	Coords chainedBoundary = cm.getOrderedCoords();

	std::cout << "DEBUG: " << (m_boundary->size() != chainedBoundary.size() ? "FAIL" : "success") << " :: size(m_boundary) = " << m_boundary->size() << ", size(chainedBoundary) = " << chainedBoundary.size() << std::endl;

	// non-max suppress
	Coord next = chainedBoundary.front();
	Coord curr = chainedBoundary.back();
	Coord prev;
	for ( uint i = 1; i < chainedBoundary.size(); ++i )
	{
		prev = curr;
		curr = next;
		next = chainedBoundary[i];
		if ( grad(curr.x(),curr.y(),0) >= grad(prev.x(), prev.y(),0) && grad(curr.x(),curr.y(),0) > grad(next.x(), next.y(),0) )
			m_boundaryGrad.push_back(CoordFFF(curr, grad(curr.x(),curr.y(),0), grad(curr.x(),curr.y(),1), grad(curr.x(),curr.y(),2)));
	}
	// last
	prev = curr;
	curr = next;
	next = chainedBoundary.front();
	if ( grad(curr.x(),curr.y(),0) >= grad(prev.x(), prev.y(),0) && grad(curr.x(),curr.y(),0) > grad(next.x(), next.y(),0) )
		m_boundaryGrad.push_back(CoordFFF(curr, grad(curr.x(),curr.y(),0), grad(curr.x(),curr.y(),1), grad(curr.x(),curr.y(),2)));

	//std::cout << "DEBUG: " << (chainedBoundary.size() != m_boundaryGrad.size() ? "FAIL" : "success") << " - size(chainedBoundary) = " << chainedBoundary.size() << ", size(m_boundaryGrad) = " << m_boundaryGrad.size() << std::endl;
	std::cout << "Number of significant gradients = " << m_boundaryGrad.size() << std::endl;

	// note which side each boundary coord is on
	m_maxGrad = 0.f;
	for ( uint i = 0; i < m_boundaryGrad.size(); ++i )
	{
		Coord coord = std::get<0>(m_boundaryGrad[i]);
		float mag = std::get<1>(m_boundaryGrad[i]);
		float nx = std::get<2>(m_boundaryGrad[i]);
		float ny = std::get<3>(m_boundaryGrad[i]);

		if ( mag > m_maxGrad ) m_maxGrad = mag;

		std::vector<int> info(m_boundary->size());
		for ( uint j = 0; j < m_boundary->size(); ++j )
		{
			Coord c = (*m_boundary)[j].first;
			float dx = c.x() - coord.x();
			float dy = c.y() - coord.y();
			float dot = dx * nx + dy * ny;
			info[j] = dot < 0 ? -1 : dot > 0 ? 1 : 0;
		}
		m_boundarySides.push_back(info);
	}

}

Color GradientWeightedInpainter::pixelColor(Coord _c)
{
	std::vector<double> ks(m_boundary->size(), 0.f);

	uint maxd = m_attenuation.width();

	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		Coord coord = (*m_boundary)[i].first;

		uint dx = abs(_c.x() - coord.x());
		uint dy = abs(_c.y() - coord.y());

		// closer pixels more heavily weighted
		if ( dx < maxd && dy < maxd ) ks[i] = m_attenuation(dx, dy); 
		else ks[i] = m_attenuation(maxd - 1, maxd - 1);
	}

	for ( uint j = 0; j < m_boundaryGrad.size(); ++j )
	{
		Coord coord = std::get<0>(m_boundaryGrad[j]);

		float dx = _c.x() - coord.x();
		float dy = _c.y() - coord.y();

		float mag = std::get<1>(m_boundaryGrad[j]);
		float nx = std::get<2>(m_boundaryGrad[j]);
		float ny = std::get<3>(m_boundaryGrad[j]);

		float k_atnuation = pow(1.f - mag / m_maxGrad, 2.f);

		float dot = dx * nx + dy * ny;
		int s = (dot < 0 ? -1 : (dot > 0 ? 1 : 0));

		for ( uint i = 0; i < m_boundary->size(); ++i )
		{
			float s_i = m_boundarySides[j][i];

			ks[i] *= ( s_i == s ) ? 1.f : k_atnuation;
		}
	}

	double r = 0;
	double g = 0;
	double b = 0;
	double k = 0;

	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		Color c_i = (*m_boundary)[i].second;

		float k_i = std::max(ks[i], double(std::numeric_limits<float>::min()));

		if ( rand() % 100 >= m_jitter )
		{
			r += k_i * float(c_i.r);
			g += k_i * float(c_i.g);
			b += k_i * float(c_i.b);
			k += k_i;
		}
	}

	Color c;
	c.r = uchar(r / k);
	c.g = uchar(g / k);
	c.b = uchar(b / k);

	return c;
}

