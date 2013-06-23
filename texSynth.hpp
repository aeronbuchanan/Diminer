/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of TexSynth, a digital inpainting resource.
 *
 * TexSynth is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TexSynth is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TexSynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <limits>
#include <algorithm>

#include "circularSeam.h"

template<uint N>
uint TexSynth::TexSynther<N>::selectPatchFor(cimg_library::CImg<float> & _I, cimg_library::CImg<float> & _M, uint _x, uint _y) const
{
	// printf("selectPatchFor\n");

	if ( m_patches.size() == 0 ) return 0; // TODO: throw, etc.

	typedef std::pair<double, uint> DUIP;
	struct compareDUIP { bool operator()(DUIP const & a, DUIP const & b) { return a.first > b.first; } };

	Patch p(_I, _x, _y);
	Patch m(_M, _x, _y);

	std::vector<DUIP> heap;
	double threshold = std::numeric_limits<double>::max();

	for ( uint i = 0; i < m_patches.size(); ++i )
	{
		double d = Patch::maskedSqrdError(p, m_patches[i], m);

		if ( d < threshold )
		{
			heap.push_back(std::make_pair(d, i));
			std::push_heap(heap.begin(), heap.end(), compareDUIP());
			if ( d != 0 )
			{
				// set threshold as multiple of second highest, non zero error
				uint j = 1;
				while ( j < heap.size() && heap[j].first == 0 )
					j++;
				if ( j < heap.size() )
					threshold = heap[j].first * 1.1;
			}
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compareDUIP());

	int i = heap.size(); // heap is in reverse order
	uint count = 0;
	uint CANDIDATE_MAX = heap.size();

	while (	count < CANDIDATE_MAX && i > 0 && heap[--i].first <= threshold	)
		++count;

	printf("DEBUG stats: %lu in total; %lu in selection; %d under threshold (%f); ", m_patches.size(), heap.size(), count, threshold);

	if ( count == 0 ) return printf("ACHTUNG!\n"); // TODO: throw, etc.

	uint s = heap.size() - 1 - (rand() % count); // selection

	printf("final selection = %i\n", heap[s].second);

	return heap[s].second;
}

template<uint N>
void TexSynth::TexSynther<N>::extendTextureIn(cimg_library::CImg<float> & _img, cimg_library::CImg<float> & _msk) const
{
	//printf("extendTextureIn\n");

	Patch ones(255.f);

	CImg<float> debug(_img);

	uint shift = m_patchWidth * 0.6;

	printf("I(%dx%d) and %d\n", _img.width(), _img.height(), m_patchWidth);

	uint iii = 0;
	char name[100];

	for ( int y = 0; y < _img.height(); y += shift )
	{
		uint yy = y;
		if ( yy >= _img.height() - m_patchWidth )
		{
			yy = _img.height() - m_patchWidth;
			y = _img.height();
		}

		for ( int x = 0; x < _img.width(); x += shift )
		{
			uint xx = x;
			if ( xx >= _img.width() - m_patchWidth )
			{
				xx = _img.width() - m_patchWidth;
				x = _img.width();
			}

			Patch p = m_patches[selectPatchFor(_img, _msk, xx, yy)];

			// DEBUG
			p.insertInto(debug, xx, yy);

			// TODO: put the following somewhere more general
			Patch image(_img, xx, yy);
			Patch mask(_msk, xx, yy);
			Patch diffs = Patch::diffSqrd(p, image);

			Table<float> tMask = CircSeams::findMin( diffs.asTable(), mask.channelAsTable(0) );

			Patch blendMask(tMask);

			//sprintf(name, "debug_%s_%d.png", "blendMask", iii);
			//blendMask.save(name);

			// Merge patches
			p = p.blendedWith(image, blendMask);
			p.insertInto(_img, xx, yy);

			ones.insertInto(_msk, xx, yy);

			//sprintf(name, "debug_%s_%d.png", "fullImage", iii);
			//_img.save(name);

			iii++;
		}
	}

	debug.save("debug.png");
}

