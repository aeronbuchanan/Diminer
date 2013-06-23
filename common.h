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

typedef unsigned char uchar;

namespace TexSynth
{

//! Useful debug typename identifiers
template<typename T> std::string identifyType(T const &) { return "unknown"; }
template<> std::string identifyType(int const &) { return "int"; }
template<> std::string identifyType(uint const &) { return "uint"; }
template<> std::string identifyType(char const &) { return "char"; }
template<> std::string identifyType(uchar const &) { return "uchar"; }
template<> std::string identifyType(float const &) { return "float"; }
template<> std::string identifyType(double const &) { return "double"; }

} // end namespace TexSynth
