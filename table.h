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
#include <iomanip>
#include "vecn.h"

namespace TexSynth
{

//! Another matrix class
template<typename T>
class Table
{
private:
	// ordering is row wise (0,0) (1,0) (2,0) (3,0) ...
	uint kCoord(uint _i, uint _j) const { return (_j * m_width) + _i; } //!< convert (i, j) coord to linear index
	//uint iCoord(uint _k) const { return ( _k / m_width ) % m_height; } //!< convert linear index to horizontal patch coord
	//uint jCoord(uint _k) const { return ( _k / m_width ) / m_height; } //!< convert linear index to vertical patch coord

	size_t m_width;
	size_t m_height;
	std::vector<T> m_vec;

	void init(T const & _v) { m_vec.resize(m_width * m_height, _v); }
	void init(uint _width, uint _height) { m_width = _width; m_height = _height; init(T()); }

public:
	Table(uint _width, uint _height) : m_width(_width), m_height(_height) { init(T()); }
	Table(uint _width, uint _height, T const & _v) : m_width(_width), m_height(_height) { init(_v); }

	template<typename U>
	explicit Table(Table<U> const & _t)
	{
		init(_t.width(), _t.height());
		for ( uint i = 0; i < m_vec.size(); ++i )
			m_vec[i] = _t[i];
	}

	// TODO: use size_t everywhere else too
	size_t size() const { return m_width * m_height; }
	size_t width() const { return m_width; }
	size_t height() const { return m_height; }

	T & operator[](uint _k) { return m_vec[_k]; }
	T operator[](uint _k) const { return m_vec[_k]; }

	T & operator()(uint _i, uint _j) { return m_vec[kCoord(_i, _j)]; }
	T operator()(uint _i, uint _j) const { return m_vec[kCoord(_i, _j)]; }

	Table<T> operator*(float _v) const { Table<T> t(*this); for ( uint i = 0; i < t.size(); ++i ) t[i] *= _v; return t; }

	void display()
	{
		display(0, width() - 1, 0, height() - 1);
	}

	void display(uint si, uint ei, uint sj, uint ej)
	{
		//std::cout.precision(2);
		printf("Table<%s>: %lux%lu\n", identifyType(T()).c_str(), m_width, m_height);
		for ( uint j = sj; j <= ej; ++j )
		{
			for ( uint i = si; i <= ei - 1; ++i )
				std::cout << std::setw(10) << (*this)(i, j) << ", ";
			std::cout << std::setw(10) << (*this)(ei, j) << std::endl;
		}
		//std::cout.unsetf( std::ios::floatfield );
	}
};


} // end namespace TexSynth
