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

#include <iostream>
#include "CImg.h"
using cimg_library::CImg;

typedef unsigned char uchar;

namespace TexSynth
{

//! Another fixed vector class
template<class T, uint N>
class VecN
{

#define FOREACH_(I,N) for( uint (I) = 0; (I) < (N); ++(I) )

public:
	VecN() { setAll(T()); }
	explicit VecN(T const _v) { setAll(_v); }
	explicit VecN(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } }
	VecN(VecN<T,N> const & _v) { *this = _v; }
	template<class Y, uint M>
	VecN(VecN<Y, M> const & _v) { setAll(0); FOREACH_(i, std::min(N,M)) { m_v[i] = _v[i]; } }

	VecN<T,N> & operator=(T const _v) { return setAll(_v); }
	VecN<T,N> & operator=(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }
	VecN<T,N> & operator=(VecN<T,N> const & _v) { FOREACH_(i,N) { m_v[i] = _v.m_v[i]; } return *this; }
	template<class Y>
	VecN<T,N> & operator=(VecN<Y,N> const & _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }

	T & operator[](uint _i) { return m_v[std::min(_i, N - 1)]; }
	T operator[](uint _i) const { return m_v[std::min(_i, N - 1)]; }

	template<class Y>
	VecN<T,N> & operator+=(VecN<Y,N> const & _v) { FOREACH_(i,N) { m_v[i] += _v[i]; } return *this; }
	VecN<T,N> & operator+=(T _v) { return *this += VecN<T,N>(_v); }

	template<class Y>
	VecN<T,N> & operator-=(VecN<Y,N> const & _v) { return *this += (-_v); }
	VecN<T,N> & operator-=(T _v) { return *this -= VecN<T,N>(_v); }

	VecN<T,N> & operator*=(T _v) { FOREACH_(i,N) { m_v[i] *= _v; } return *this; }
	VecN<T,N> & operator/=(T _v) { FOREACH_(i,N) { m_v[i] /= _v; } return *this; }

	template<class Y>
	VecN<T,N> operator+(VecN<Y,N> const & _v) const { return VecN<Y,N>(*this) += _v; }
	VecN<T,N> operator+(T _v) const { return VecN<T,N>(*this) += _v; }

	template<class Y>
	VecN<T,N> operator-(VecN<Y,N> const & _v) const { return VecN<Y,N>(*this) -= _v; }
	VecN<T,N> operator-(T _v) const { return VecN<T,N>(*this) -= _v; }

	VecN<T,N> operator*(T _v) const { return VecN<T,N>(*this) *= _v; }
	VecN<T,N> operator/(T _v) const { return VecN<T,N>(*this) /= _v; }

	VecN<T,N> operator-() const { return VecN<T,N>(*this) *= -1; }

	VecN<T,N> & operator++() { return *this += 1; }
	VecN<T,N> & operator--() { return *this -= 1; }

	VecN<T,N> operator++(int) { VecN<T,N> pre(*this); ++(*this); return pre; }
	VecN<T,N> operator--(int) { VecN<T,N> pre(*this); --(*this); return pre; }

	template<class Y>
	T dot(VecN<Y,N> const & _v) const { T k = T(); FOREACH_(i,N) { k += m_v[i] * _v[i]; } return k; }

	//! Element-wise multiplication (dot-multiply of matlab)
	template<class Y>
	VecN<T,N> maskWith(VecN<Y,N> const & _m) { VecN<T,N> v(*this); FOREACH_(i,N) { v[i] *= _m[i]; } return v; }

	T magnitudeSqrd() const { return this->dot(*this); }
	T magnitude() const { return sqrt(magnitudeSqrd()); }
	template<class Y>
	T maskedMagSqrd(VecN<Y,N> const & _m) const { return this->maskWith(_m).dot(*this); }

	VecN<T,N> abs() const { VecN<T,N> v; FOREACH_(i,N) { v[i] = std::abs(m_v[i]); } return v; }

	void print() const { std::cout << "("; FOREACH_(i, (N - 1)) { std::cout << m_v[i] << ", "; } std::cout << m_v[N - 1] << ")"; }
	void printn() const { print(); std::cout << std::endl; }

private:
	VecN<T,N> & setAll(T _v) { FOREACH_(i,N) { m_v[i] = _v; } return *this; }

	T m_v[N];

#undef FOREACH_
};

template<uint N>
class ImagePatch
{
	typedef VecN<float, N * N * 3> VEC;

	static inline uint xCoord(uint _x, uint, uint _i, uint) { return _x + _i; }
	static inline uint yCoord(uint, uint _y, uint, uint _j) { return _y + _j; }

	VEC m_v;

public:
	ImagePatch() {}
	explicit ImagePatch(float _v) : m_v(_v) {}
	template<class U>
	ImagePatch(CImg<U> const & _img, uint _x, uint _y) { this->extractFrom(_img, _x, _y); }

	template<class U>
	bool extractFrom(CImg<U> const & _img, uint _x, uint _y)
	{
		bool badCoords = false;

		if ( _x > _img.width() - N - 1 ) { printf(" [extract: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N - 1 ) { printf(" [extract: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			uint k = 0;
			for ( uint j = 0; j < N; ++j )
				for ( uint i = 0; i < N; ++i )
					for ( uint c = 0; c < 3; ++c )
						m_v[k++] = _img(xCoord(_x,_y,i,j), yCoord(_x,_y,i,j), 0, c);
		}

		return badCoords;
	}

	template<class U>
	bool insertInto(CImg<U> & _img, uint _x, uint _y) const
	{
		bool badCoords = false;

		if ( _x > _img.width() - N - 1 ) { printf(" [insert: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N - 1 ) { printf(" [insert: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			uint k = 0;
			for ( uint j = 0; j < N; ++j )
				for ( uint i = 0; i < N; ++i )
					for ( uint c = 0; c < 3; ++c )
						_img(xCoord(_x,_y,i,j), yCoord(_x,_y,i,j), 0, c) = m_v[k++];
		}

		return badCoords;
	}

	ImagePatch<N> & maskWith(ImagePatch<N> _m) { m_v.maskWith(_m.m_v); return *this; }

	template<uint M>
	static float sqrdError(ImagePatch<M> const & a, ImagePatch<M> const & b) { return (a.m_v - b.m_v).magnitudeSqrd(); }

	template<uint M>
	static float maskedSqrdError(ImagePatch<M> const & a, ImagePatch<M> const & b, ImagePatch<M> const & m) { return (a.m_v - b.m_v).maskWith(m.m_v).magnitudeSqrd(); }

};

} // end namespace TexSynth
