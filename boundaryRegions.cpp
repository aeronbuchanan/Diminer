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


#include "boundaryRegions.h"

using namespace Diminer;

double Diminer::findSalientPoints(Coords const * const boundaryCoords, GradImage const * const gradImag, SalientPointRanges & output)
{
	// non-max suppress
	// TODO: check for and deal with looped boundaries
	// TODO: measure "peakiness" and also supress insufficiently peaky maximums

	// if not a loop, include ends
	output.push_back(boundaryCoords->begin());

	double maxGrad = 0.f;
	auto prev = boundaryCoords->begin();
	auto curr = prev + 1;
	auto next = curr + 1;
	for ( ; next != boundaryCoords->end(); prev++, curr++, next++ )
	{
		int x = (*curr)->x;
		int y = (*curr)->y;
		double mag = (*gradImag)(x, y, 0, 0);
		if ( mag >= (*gradImag)((*prev)->x, (*prev)->y, 0, 0) && mag > (*gradImag)((*next)->x, (*next)->y, 0, 0) )
		{
			output.push_back( curr );
			if ( mag > maxGrad ) maxGrad = mag;
		}
	}

	if ( boundaryCoords->size() > 1 )
		output.push_back(boundaryCoords->end() - 1);

	std::cout << "Found " << output.size() << " worthy gradient points" << std::endl;

	return maxGrad;
}




