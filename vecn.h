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

namespace TexSynth
{

//! Another fixed vector class
template<class T, uint N>
class VecN
{

#define FOREACH_(I,N) for( uint (I) = 0; (I) < (N); ++(I) )

public:
	VecN() { setAll(T()); }
	explicit VecN(T const & _v) { setAll(_v); }
	explicit VecN(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } }

	template<class Y, uint M>
	explicit VecN(VecN<Y, M> const & _v) { setAll(0); FOREACH_(i, std::min(N,M)) { m_v[i] = _v[i]; } }

	// TODO: accept {} list constructors
	VecN(T const & _v0, T const & _v1)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
	}
	VecN(T const & _v0, T const & _v1, T const & _v2)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
		getElement<2>() = _v2;
	}
	VecN(T const & _v0, T const & _v1, T const & _v2, T const & _v3)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
		getElement<2>() = _v2;
		getElement<3>() = _v3;
	}

	VecN<T, N> & operator=(T const _v) { return setAll(_v); }
	VecN<T, N> & operator=(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }
	template<class Y>
	VecN<T, N> & operator=(VecN<Y,N> const & _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }

	// TODO: WARNING: silently returns last element if index out-of-bounds
	T & operator[](uint _i) { return m_v[std::min(_i, N - 1)]; }
	T operator[](uint _i) const { return m_v[std::min(_i, N - 1)]; }

private:
	template<uint I>
	T & getElement() { static_assert(I < N, "Index out-of-bounds."); return m_v[I]; }
	template<uint I>
	T getElement() const { static_assert(I < N, "Index out-of-bounds."); return m_v[I]; }

public:
	T & x() { return getElement<0>(); }
	T x() const { return getElement<0>(); }
	T & y() { return getElement<1>(); }
	T y() const { return getElement<1>(); }
	T & z() { return getElement<2>(); }
	T z() const { return getElement<2>(); }
	T & w() { return getElement<3>(); }
	T w() const { return getElement<3>(); }

	template<class Y>
	VecN<T, N> & operator+=(VecN<Y, N> const & _v) { FOREACH_(i,N) { m_v[i] += _v[i]; } return *this; }
	VecN<T, N> & operator+=(T _v) { return *this += VecN<T, N>(_v); }

	template<class Y>
	VecN<T, N> & operator-=(VecN<Y, N> const & _v) { return *this += (-_v); }
	VecN<T, N> & operator-=(T _v) { return *this -= VecN<T, N>(_v); }

	VecN<T, N> & operator*=(T _v) { FOREACH_(i,N) { m_v[i] *= _v; } return *this; }
	VecN<T, N> & operator/=(T _v) { FOREACH_(i,N) { m_v[i] /= _v; } return *this; }

	template<class Y>
	VecN<T, N> operator+(VecN<Y,N> const & _v) const { return VecN<Y,N>(*this) += _v; }
	VecN<T, N> operator+(T _v) const { return VecN<T, N>(*this) += _v; }

	template<class Y>
	VecN<T, N> operator-(VecN<Y, N> const & _v) const { return VecN<Y,N>(*this) -= _v; }
	VecN<T, N> operator-(T _v) const { return VecN<T, N>(*this) -= _v; }

	VecN<T, N> operator*(T _v) const { return VecN<T, N>(*this) *= _v; }
	VecN<T, N> operator/(T _v) const { return VecN<T, N>(*this) /= _v; }

	VecN<T, N> operator-() const { return VecN<T, N>(*this) *= -1; }

	VecN<T, N> & operator++() { return *this += 1; }
	VecN<T, N> & operator--() { return *this -= 1; }

	VecN<T, N> operator++(int) { VecN<T, N> pre(*this); ++(*this); return pre; }
	VecN<T, N> operator--(int) { VecN<T, N> pre(*this); --(*this); return pre; }

	template<class Y>
	T dot(VecN<Y,N> const & _v) const { T k = T(); FOREACH_(i,N) { k += m_v[i] * _v[i]; } return k; }

	//! Element-wise multiplication (dot-multiply of matlab)
	template<class Y>
	VecN<T, N> maskedWith(VecN<Y, N> const & _m) const { VecN<T, N> v(*this); FOREACH_(i,N) { v[i] *= _m[i]; } return v; }
	//! In-place element-wise multiplication (dot-multiply of matlab)
	template<class Y>
	VecN<T, N> & maskWith(VecN<Y, N> const & _m) { FOREACH_(i,N) { m_v[i] *= _m[i]; } return *this; }

	T magnitudeSqrd() const { return this->dot(*this); }
	T magnitude() const { return sqrt(magnitudeSqrd()); }
	template<class Y>
	T maskedMagSqrd(VecN<Y, N> const & _m) const { return this->maskedWith(_m).dot(*this); }

	VecN<T, N> & makeAbsolute() { FOREACH_(i,N) { m_v[i] = std::abs(m_v[i]); } return *this; }
	VecN<T, N> asAbsolute() const { return VecN<T, N>(*this).makeAbsolute(); }

	void print() const { std::cout << "("; FOREACH_(i, (N - 1)) { std::cout << m_v[i] << ", "; } std::cout << m_v[N - 1] << ")"; }
	void printn() const { print(); std::cout << std::endl; }

	static inline uint size() { return N; }

	typedef T Type;

protected:
	VecN<T, N> & setAll(T _v) { FOREACH_(i,N) { m_v[i] = _v; } return *this; }

	T m_v[N];

#undef FOREACH_
};

typedef VecN<int, 2> Coord;


} // end namespace TexSynth
