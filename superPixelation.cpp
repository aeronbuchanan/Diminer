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


#include "superPixelation.h"

using namespace Diminer;


SuperPixelator::SuperPixelator(
	SourceImage const * const _img,
	GradImage const * const _grads,
	
	) : m_img(SourceImage(*_img)), m_gradImg(GradImage(*_grads))
{
	// prepping for k-means super-pixel clustering
	std::cout << "K-means" << std::endl;

	double GRAD_THRESHOLD = 50;
#define NEXT_SIGNIFICANT_GRAD_POINT(itr) while(itr != m_maxGradPoints.end() && m_gradImg((**itr)->x,(**itr)->y,0,0) < GRAD_THRESHOLD) itr++;
	std::vector<std::shared_ptr<KPoint> > cs;
	auto gCurr = m_maxGradPoints.begin();
	auto gNext = gCurr + 1;
	NEXT_SIGNIFICANT_GRAD_POINT(gNext)
	// TODO: shift secnd to first
	while ( gNext != m_maxGradPoints.end() )
	{
		auto midpnt = *gCurr;
		bool flip = false;
		for ( auto skpptr = *gCurr; skpptr != *gNext; skpptr++, flip = !flip )
			if ( flip ) midpnt++;

		// get mid-point
		int mx = (*midpnt)->x;
		int my = (*midpnt)->y;

		tweekCoords(mx, my);

		double c1 = m_img(mx, my, 0, 0);
		double c2 = m_img(mx, my, 0, 1);
		double c3 = m_img(mx, my, 0, 2);

		cs.push_back(std::make_shared<KPoint>(c1, c2, c3, mx, my));

		gCurr = gNext;
		gNext++;
		if ( *gNext != m_maxGradPoints.back() )
			NEXT_SIGNIFICANT_GRAD_POINT(gNext)
	}

	int spacing = ceil(double(m_boundary->size()) / double(cs.size()));
	int spacingSqrd = spacing * spacing;
	double spatialBias = 100.f / double(spacingSqrd);

	std::cout << "Initialized boundary with " << cs.size() << " seed points (spacing = " << spacing << ")" << std::endl;

	std::vector<std::shared_ptr<KPoint> > xs;
	for ( int x = spacing / 2; x < W; x += spacing )
	{
		for ( int y = spacing / 2; y < H; y += spacing )
		{
			if ( ! mask(x, y) )
			{
				bool notClose = true;
				for ( auto ci = cs.begin(); notClose && ci != cs.end(); ci++ )
				{
					double dx = x - (*ci)->x;
					double dy = y - (*ci)->y;
					notClose = (dx * dx + dy * dy) >= spacingSqrd;
				}
				if ( notClose )
				{
					int nx = x;
					int ny = y;
					tweekCoords(nx, ny);

					xs.push_back(std::make_shared<KPoint>(m_img(nx, ny, 0, 0), m_img(nx, ny, 0, 1), m_img(nx, ny, 0, 2), nx, ny));
				}
			}
		}
	}

	cs.insert(cs.end(), xs.begin(), xs.end());

	std::cout << "Added " << xs.size() << " interior grid seed points" << std::endl;

	CImg<float> labeling(W, H, 1, 2, std::numeric_limits<float>::max());


	// perform 10 steps of the k-means algo
	std::cout << "Starting Algorithm" << std::endl;


};



void SuperPixelator::tweekCenter(int & x, int & y)
{
	// slip to a pixel of lower gradient if possible
	double minGrad = m_gradImg(x, y, 0, 0);
	int nx = x;
	int ny = y;
	for ( int dx = -1; dx <= 1; dx++ )
	{
		for ( int dy = -1; dy <= 1; dy++ )
		{
			if ( dx != 0 || dy != 0 )
			{
				int gx = x + dx;
				int gy = y + dy;
				if ( gx >= 0 && gx < W && gy >= 0 && gy < H && ! mask(gx, gy) && (*m_gradImg)(gx, gy, 0, 0) < minGrad )
				{
					minGrad = (*m_gradImg)(gx, gy, 0, 0);
					nx = gx;
					ny = gy;
				}
			}
		}
	}
	x = nx;
	y = ny;
};

KMeans::iterate(int numSteps)
{

	for ( int j = 0; j < NUM_STEPS; j++ )
	{
		std::cout << "Step " << j << std::endl;
		// update points
		cimg_forXY(m_labeling, x, y)
		{
			if ( ! mask(x, y) )
			{
				std::shared_ptr<KPoint> kp = std::make_shared<KPoint>(m_img(x, y, 0, 0), m_img(x, y, 0, 1), m_img(x, y, 0, 2), x, y);
				int i = 0;
				for ( auto ci = cs.begin(); ci != cs.end(); ci++, i++ )
				{
					float d = KPointDistance(**ci, *kp, m_spatialBias);
					if ( d < labeling(x, y, 0, 1) )
					{
						m_labeling(x, y, 0, 0) = i;
						m_labeling(x, y, 0, 1) = d;
					}
				}
			}
		}
		// update means
		if ( j < NUM_STEPS - 1 )
		{
			std::vector<int> counts(cs.size(), 0);
			for ( auto ci = cs.begin(); ci != cs.end(); ci++ ) { (*ci)->c1 = 0; (*ci)->c2 = 0; (*ci)->c3 = 0; (*ci)->x = 0; (*ci)->y = 0; }
			cimg_forXY(labeling, x, y)
			{
				if ( ! mask(x, y) )
				{
					int k = labeling(x, y, 0, 0);
					cs[k]->c1 += img(x, y, 0, 0);
					cs[k]->c2 += img(x, y, 0, 1);
					cs[k]->c3 += img(x, y, 0, 2);
					cs[k]->x += x;
					cs[k]->y += y;
					counts[k]++;
				}
			}
			int i = 0;
			for ( auto ci = cs.begin(); ci != cs.end(); ci++, i++ )
			{
				double count = counts[i];
				(*ci)->c1 /= count;
				(*ci)->c2 /= count;
				(*ci)->c3 /= count;
				(*ci)->x /= count;
				(*ci)->y /= count;
			}
		}
	}

	std::cout << "Getting boundary chains" << std::endl;



	std::cout << "Finished" << std::endl;



}
