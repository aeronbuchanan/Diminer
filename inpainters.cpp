/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of Diminutive, an digital inpainting resource.
 *
 * Diminutive is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Diminutive is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Diminutive.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "inpainters.h"

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
		float dx = _c.first - coord.first;
		float dy = _c.second - coord.second;

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
		float dx = _c.first - coord.first;
		float dy = _c.second - coord.second;
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

void GradientWeightedInpainter::init(CImg<uchar> const * const _img)
{
	CImg<uchar> img(*_img);
	img.blur(1.5);

	// greyscale
	CImg<float> g(img.width(), img.height(), 1, 1, 0.f);
	cimg_forXYC(img, x,y,k)
		g(x,y) += img(x,y,k);

	// gradient
	// TODO: take note of mask
	CImg<float> grad(img.width(), img.height(), 3, 1, 0.f);

	CImg_3x3(I, float);
	cimg_for3x3(g, x, y, 0, 0, I, float)
	{
		float gx = (Inc - Ipc);
		float gy = (Icn - Icp);
		float mag = sqrt( gx * gx + gy * gy );
		grad(x,y,0) = mag;
		grad(x,y,1) = gx / mag;
		grad(x,y,2) = gy / mag;
	}

	std::cout << "DEBUG: length of boundary = " << m_boundary->size() << std::endl;

#ifdef DEBUGIMG
	CImg<uchar> debugImg(img.width(), img.height(), 1, 3, 0.f);
	for ( uint i = 0; i < m_boundary->size(); ++i )
	{
		Coord c = (*m_boundary)[i].first;
		float x = c.first;
		float y = c.second;
		debugImg(x,y,0,0) = 255;
		debugImg(x,y,0,1) = 255;
		debugImg(x,y,0,2) = 255;
	}
	debugDisp->update(&debugImg);
#endif

	// chain boundary coords in order
	std::vector<bool> chained(m_boundary->size(), false);
	Coords chainedBoundary;
	chainedBoundary.push_back((*m_boundary)[0].first);
	chained[0] = true;
	uint last_j = 0;
	int j = 1;
	int jDirection = ( j > last_j ? 1 : -1 );
	for ( uint i = 1; i < m_boundary->size(); ++i )
	{
		Coord c = chainedBoundary.back();
		int x = c.first;
		int y = c.second;
		bool nextFound = false;
		uint first_j = j;

		//std::cout << "SET(" << first_j << ") -> ";

		uint diag_j = 0;
		bool firstIter = true;
		while ( !nextFound )
		{
			//std::cout << j;

			if ( j == first_j && !firstIter )
			{
				if ( diag_j )
				{
					chainedBoundary.push_back((*m_boundary)[diag_j].first);
					chained[diag_j] = true;
					nextFound = true;
					jDirection = ( j > last_j ? 1 : -1 );
					last_j = j;

#ifdef DEBUGIMG
					int xx = chainedBoundary.back().first;
					int yy = chainedBoundary.back().second;
					debugImg(xx,yy,0,0) = 0;
					debugImg(xx,yy,0,1) = 0;
					debugImg(xx,yy,0,2) = 255;
					//debugDisp->update(&debugImg);
#endif

					//std::cout << "added(" << diag_j << ") " << std::endl;
				}
				else
				{
					// There are sections of boundary with mask on both sides
					std::cout << "ERROR: chain ended prematurely" << std::endl;
					i = m_boundary->size();
					break;
				}
			}
			else
			{
				if ( !chained[j] )
				{
					int dx = abs( x - (*m_boundary)[j].first.first );
					int dy = abs( y - (*m_boundary)[j].first.second );

					//std::cout << "(" << dx << "," << dy << ") ";

					// region is 4-connected, therefore next boundary will be
					// \Delta = (1,0) or (0,1) or (1,1) only.
					if ( (dx == 1 && dy == 0) || (dx == 0 && dy == 1) )
					{
						chainedBoundary.push_back((*m_boundary)[j].first);
						chained[j] = true;
						nextFound = true;
						jDirection = ( j > last_j ? 1 : -1 );
						last_j = j;

#ifdef DEBUGIMG
						int xx = chainedBoundary.back().first;
						int yy = chainedBoundary.back().second;
						debugImg(xx,yy,0,0) = 0;
						debugImg(xx,yy,0,1) = 0;
						debugImg(xx,yy,0,2) = 255;
						//debugDisp->update(&debugImg);
#endif

						//std::cout << "added(" << j << ")" << std::endl;
					}
					else if ( (dx == 1) && (dy == 1) )
					{
						diag_j = j;
						//std::cout << "marked(" << j << ") ";
					}
				}

				else
				{
					//std::cout << "* ";
				}
			}

			if ( firstIter ) firstIter = false;

			j += jDirection;
			if ( j < 0 ) j = m_boundary->size() - 1;
			else if ( j >= m_boundary->size() ) j = 0;
		}
	}
	std::cout << "DEBUG: " << (m_boundary->size() != chainedBoundary.size() ? "FAIL" : "success") << " - size(m_boundary) = " << m_boundary->size() << ", size(chainedBoundary) = " << chainedBoundary.size() << std::endl;

	// non-max suppress
	Coord next = chainedBoundary.front();
	Coord curr = chainedBoundary.back();
	Coord prev;
	for ( uint i = 1; i < chainedBoundary.size(); ++i )
	{
		prev = curr;
		curr = next;
		next = chainedBoundary[i];
		if ( grad(curr.first,curr.second,0) >= grad(prev.first, prev.second,0) && grad(curr.first,curr.second,0) > grad(next.first, next.second,0) )
			m_boundaryGrad.push_back(CoordFFF(curr, grad(curr.first,curr.second,0), grad(curr.first,curr.second,1), grad(curr.first,curr.second,2)));
	}
	// last
	prev = curr;
	curr = next;
	next = chainedBoundary.front();
	if ( grad(curr.first,curr.second,0) >= grad(prev.first, prev.second,0) && grad(curr.first,curr.second,0) > grad(next.first, next.second,0) )
		m_boundaryGrad.push_back(CoordFFF(curr, grad(curr.first,curr.second,0), grad(curr.first,curr.second,1), grad(curr.first,curr.second,2)));

	//std::cout << "DEBUG: " << (chainedBoundary.size() != m_boundaryGrad.size() ? "FAIL" : "success") << " - size(chainedBoundary) = " << chainedBoundary.size() << ", size(m_boundaryGrad) = " << m_boundaryGrad.size() << std::endl;
	std::cout << "Number of significant gradients = " << m_boundaryGrad.size() << std::endl;

	// note which side each boundary coord is on
	for ( uint i = 0; i < m_boundaryGrad.size(); ++i )
	{
		Coord coord = std::get<0>(m_boundaryGrad[i]);
		float nx = std::get<2>(m_boundaryGrad[i]);
		float ny = std::get<3>(m_boundaryGrad[i]);

		std::vector<int> info(m_boundary->size());
		for ( uint j = 0; j < m_boundary->size(); ++j )
		{
			Coord c = (*m_boundary)[j].first;
			float dx = c.first - coord.first;
			float dy = c.second - coord.second;
			float dot = dx * nx + dy * ny;
			info[j] = dot < 0 ? -1 : dot > 0 ? 1 : 0;
		}
		m_boundarySides.push_back(info);
	}

}

Color GradientWeightedInpainter::pixelColor(Coord _c)
{
	double r = 0;
	double g = 0;
	double b = 0;
	double k = 0;

	for ( uint j = 0; j < m_boundaryGrad.size(); ++j )
	{
		Coord coord = std::get<0>(m_boundaryGrad[j]);

		float dx = _c.first - coord.first;
		float dy = _c.second - coord.second;

		float mag = std::get<1>(m_boundaryGrad[j]);
		float nx = std::get<2>(m_boundaryGrad[j]);
		float ny = std::get<3>(m_boundaryGrad[j]);

		float k_mag = pow(mag / 255.f, 0.1); // TODO: make this a parameter
		float k_revMag = pow(1.f - mag / 255.f, 0.35); // TODO: make this a parameter

		float dot = dx * nx + dy * ny;
		int s = (dot < 0 ? -1 : (dot > 0 ? 1 : 0));

		for ( uint i = 0; i < m_boundary->size(); ++i )
		{
			Coord coord = (*m_boundary)[i].first;
			Color c_i = (*m_boundary)[i].second;

			float dx = _c.first - coord.first;
			float dy = _c.second - coord.second;
			float distSqrd = dx * dx + dy * dy;
			float k_i = 1.f / pow(distSqrd, m_pow); // closer pixels more heavily weighted
			k_i *= k_mag;

			float s_i = m_boundarySides[j][i];

			if ( s_i == s )
			{
				//r += k_i * float(c_i.r);
				//g += k_i * float(c_i.g);
				//b += k_i * float(c_i.b);
				//k += k_i;
				r += 255;
			}
			else
			{
				b += 255;
			}
			k += 1;
		}
	}

	Color c;
	c.r = uchar(r / k);
	c.g = uchar(g / k);
	c.b = uchar(b / k);

	return c;
}

