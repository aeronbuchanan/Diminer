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

#include <iostream>
#include <string>

// TODO: make TexSynth image-library agnostic
#include "CImg.h"
using cimg_library::CImg;

#include "table.h"

namespace TexSynth
{

template<uint N, uint M = 3>
class ImagePatch
{	
protected:
	typedef VecN<float, N * N * M> ThisVecN;

	// ordering (inner -> outer) is (c,i,j)
	// optimizes conversion to gray, etc
	// bad for conversion to/from CImg (I think)
	static inline uint kCoord(uint _i, uint _j, uint _c) { return (((_j * N) + _i) * M) + _c; } //!< convert (i,j,c) coord to linear index
	static inline uint iCoord(uint _k) { return ( _k / M) % N; } //!< convert linear index to horizontal patch coord
	static inline uint jCoord(uint _k) { return ( _k / M) / N; } //!< convert linear index to vertical patch coord
	static inline uint cCoord(uint _k) { return _k % M; } //!< convert linear index to color channel
	static inline uint xCoord(uint _x, uint, uint _k) { return _x + iCoord(_k); }
	static inline uint yCoord(uint, uint _y, uint _k) { return _y + jCoord(_k); }

private:
	ThisVecN m_vec;

public:
	ImagePatch() {}
	ImagePatch(float _v) : m_vec(_v) {}
	ImagePatch(ThisVecN _v) : m_vec(_v) {}
	ImagePatch(Table<float> _t) { this->setFromTable(_t); }
	ImagePatch(CImg<float> const & _img, uint _x, uint _y) { this->extractFrom(_img, _x, _y); }

	// TODO: use size_t everywhere else too
	size_t size() const { return N * N * M; }
	size_t width() const { return N; }
	size_t height() const { return N; }
	size_t depth() const { return M; }

	float & operator[](uint _k) { return m_vec[_k]; }
	float operator[](uint _k) const { return m_vec[_k]; }

	float & operator()(uint _i, uint _j, uint _c) { return (*this)[kCoord(_i, _j, _c)]; }
	float operator()(uint _i, uint _j, uint _c) const { return (*this)[kCoord(_i, _j, _c)]; }

	ImagePatch<N, M> & setFromTable(Table<float> const & _t)
	{
		// TODO
		if ( _t.width() != N || _t.height() != N )
			printf("WARNING: Table -> ImagePatch dimension mismatch.\n");

		/*
		ImagePatch<N, M> p;
		for ( uint j = 0; j < N; ++j )
			for ( uint i = 0; i < N; ++i )
				for ( uint c = 0; c < M; ++c )
					(*this)(i, j, c) = _t(i, j);
		*/

		// assuming both are row-wise
		uint i = 0;
		uint k = 0;
		while ( i < this->size() )
		{
			(*this)[i] = _t[k];
			++i;
			if ( i % M == 0 )
				++k;
		}
		return *this;
	}

	Table<float> asTable() const { return this->grayVersion().channelAsTable(0); }
	Table<float> channelAsTable(uint _ch) const
	{
		// TODO: silently returns zero table if _ch > depth()
		Table<float> t(width(), height(), 0.f);
		// assuming both are row-wise
		uint k = 0;
		for ( uint i = 0; i < this->size(); ++i )
			if ( i % M == _ch )
				t[k++] = (*this)[i];

		return t;
	}


	bool extractFrom(CImg<float> const & _img, uint _x, uint _y)
	{
		bool badCoords = false;

		if ( _x > _img.width() - N ) { printf(" [extract: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N ) { printf(" [extract: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < ImagePatch::size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					(*this)[k] = _img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k));
		}

		return badCoords;
	}

	bool insertInto(CImg<float> & _img, uint _x, uint _y) const
	{
		bool badCoords = false;

		if ( _x > _img.width() - N ) { printf(" [insert: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N ) { printf(" [insert: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < this->size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					_img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k)) = (*this)[k];
		}

		return badCoords;
	}

	ImagePatch<N, M> & convertToGray()
	{
		if ( M > 1 )
		{
			// TODO: seems to be coming out too dark...
			for ( uint j = 0; j < N; ++j )
			{
				for ( uint i = 0; i < N; ++i )
				{
					float v = 0;
					for ( uint c = 0; c < M; ++c )
						v += (*this)[kCoord(i, j, c)];
					for ( uint c = 0; c < M; ++c )
						(*this)[kCoord(i, j, c)] = v / M;
				}
			}
		}
		return *this;
	}

	ImagePatch<N, M> grayVersion() const { return ImagePatch<N, M>(*this).convertToGray(); }

	void print() const { m_vec.print(); }
	void printn() const { m_vec.printn(); }
	void save(std::string const & _name) const { this->save(_name.c_str()); }
	void save(char * const _name) const { CImg<float> t(N, N, 1, M); insertInto(t, 0, 0); t.save(_name); }

	// TODO: refactor below code for inter-type compatibility

	ImagePatch<N, M> blendedWith(ImagePatch<N, M> const & _q, ImagePatch<N, M> const & _m) const
	{
		ImagePatch<N, M> p(*this);
		for ( uint i = 0; i < this->size(); ++i )
		{
			float a = _m[i] / 255.f;
			p[i] = (1.f - a) * _q[i] +  a * p[i];
		}
		return p;
	}

	ImagePatch<N, M> & maskWith(ImagePatch<N, M> _m) { m_vec.maskWith(_m.m_vec); return *this; }

	static float sqrdError(ImagePatch<N, M> const & _a, ImagePatch<N, M> const & _b) { return diff(_a, _b).m_v.magnitudeSqrd(); }
	static float maskedSqrdError(ImagePatch<N, M> const & _a, ImagePatch<N, M> const & _b, ImagePatch<N, M> const & _m) { return diff(_a, _b).maskWith(_m).m_vec.magnitudeSqrd(); }

	static ImagePatch<N, M> diff(ImagePatch<N, M> const & a, ImagePatch<N, M> const & b) { return ImagePatch<N, M>(a.m_vec - b.m_vec); }
	static ImagePatch<N, M> diffSqrd(ImagePatch<N, M> const & a, ImagePatch<N, M> const & b) { return ImagePatch<N, M>( (a.m_vec - b.m_vec).maskWith(a.m_vec - b.m_vec) ); }
	static ImagePatch<N, M> abs(ImagePatch<N, M> const & a) { return ImagePatch<N, M>(a.m_vec.asAbsolute()); }

/*
bool test()
{
	bool good = true;
	for ( uint j = 0; j < N; ++j )
		for ( uint i = 0; i < N; ++i )
			for ( uint c = 0; c < M; ++c )
			{
				uint k = kCoord(i,j,c);
				uint ic = iCoord(k);
				uint jc = jCoord(k);
				uint cc = cCoord(k);
				if ( i != ic || j != jc || c != cc )
				{
					good = false;
					printf("ERROR: (%u,%u,%u) -> [%u] -> (%u,%u,%u)\n", i, j, c, k, iCoord(k), jCoord(k), cCoord(k));
				}
			}
	return good;
}
*/

};

} // end namespace TexSynth
