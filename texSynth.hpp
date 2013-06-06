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

#pragma once

#include <limits>
#include <algorithm>

template<uint N>
template<class T>
uint TexSynth::TexSynther<N>::selectPatchFor(cimg_library::CImg<T> & _I,cimg_library:: CImg<T> & _M, uint _x, uint _y) const
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
		if ( heap.size() > 0) threshold = heap.front().first * 1.3; // (1.0 * m.magnitudeSqrd());
		double d = Patch::maskedSqrdError(p, m_patches[i], m);
		if ( d < threshold )
		{
			heap.push_back(std::make_pair(d, i));
			std::push_heap(heap.begin(), heap.end(), compareDUIP());
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compareDUIP());

	int i = heap.size(); // heap is in reverse order
	uint count = 0;
	uint CANDIDATE_MAX = 100;

	while (	++count < CANDIDATE_MAX && i > 0 &&	heap[--i].first < threshold	);

	--count;
	if ( count == 0 ) return 0; // TODO: throw, etc.

	uint s = heap.size() - 1 - (rand() % count); // selection

	return heap[s].second;
}


template<uint N>
template<class T>
void TexSynth::TexSynther<N>::extendTextureIn(cimg_library::CImg<T> & _I,cimg_library:: CImg<T> & _M) const
{
	//printf("extendTextureIn\n");

	Patch ones(1);

	uint shift = m_patchWidth * 0.6;

	printf("I(%dx%d) and %d\n", _I.width(), _I.height(), m_patchWidth);

	for ( uint y = 0; y < _I.height() - m_patchWidth; y += shift )
		for ( uint x = 0; x < _I.width() - m_patchWidth; x += shift )
		{
			Patch p = 	m_patches[selectPatchFor(_I, _M, x, y)];
			p.insertInto(_I, x, y);
			ones.insertInto(_M, x, y);
		}
}

