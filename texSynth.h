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

#include <utility>
#include <vector>
#include <unordered_map>

// TODO: make Diminer image-library agnostic
using cimg_library::CImg;

#include "patch.h"

namespace Diminer
{

template<uint N>
class Diminerer
{
public:
	typedef ImagePatch<N, 3> Patch;
	typedef std::vector<Patch> PatchLibrary;

	Diminerer() : m_patchWidth(N) { }

	void extendTextureIn(cimg_library::CImg<float> & _image, cimg_library::CImg<float> & _mask) const;

	uint selectPatchFor(cimg_library::CImg<float> & _image, cimg_library:: CImg<float> & _mask, uint _x, uint _y) const;

	void addPatches(PatchLibrary const & _pl) { m_patches.insert( m_patches.end(), _pl.begin(), _pl.end() ); }
	void clearPatches() { m_patches.clear(); }

private:

	PatchLibrary m_patches;
	uint m_patchWidth;
};

} // end namespace Diminer

#include "Diminer.hpp"
