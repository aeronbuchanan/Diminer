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

#include <vector>
#include <limits>
#include <algorithm>
#include <utility>

#include "CImg.h"
using cimg_library::CImg;

#include "table.h"
#include "patch.h"


namespace TexSynth
{

namespace CircSeams
{

enum Direction { NODIRECTION = 0, Downwards, Rightwards, Upwards, Leftwards };
std::ostream & operator<<(std::ostream & _stream, Direction _dir)
{
	switch ( _dir )
	{
	case ( Downwards ):_stream << "Downwards"; break;
	case ( Upwards ):_stream << "Upwards"; break;
	case ( Rightwards ):_stream << "Rightwards"; break;
	case ( Leftwards ):_stream << "Leftwards"; break;
	default:_stream << "Unknown"; break;
	}
	return _stream;
}

bool isVertical(Direction d) { return d == Downwards || d == Upwards; }
bool isHorizontal(Direction d) { return !isVertical(d); }
bool isForward(Direction d) { return d == Downwards || d == Rightwards; }
bool isReverse(Direction d) { return !isForward(d); }

// CircularSeams are traversed in an anticlockwise direction
// (with the inside on the left)
//   SO | SI   SI   EI   EI ^ EO   EO   SO
//      |F     ------>F     |R     <------R
//   EO V EI   SO   EO   SI | SO   EI   SI
// S at start; E at end in 'dir' Direction
// I on left; O on right in 'orth' Direction

struct DirCoord : public Coord
{
	DirCoord() : Coord(0, 0), direction(NODIRECTION) {}
	explicit DirCoord(Direction _d) : Coord(0, 0), direction(_d) {}
	DirCoord(Coord::Type _x, Coord::Type _y, Direction _d) : Coord(_x, _y), direction(_d) {}

	DirCoord & advIn() { orthComp() += orthStep(); return *this; }
	DirCoord & advOut() { orthComp() -= orthStep(); return *this; }
	DirCoord & advBack() { dirComp() -= dirStep(); return *this; }
	DirCoord & advForward() { dirComp() += dirStep(); return *this; }

	DirCoord & setOrth(Coord::Type _i) { orthComp() = _i; return *this; }
	DirCoord & setDir(Coord::Type _i) { dirComp() = _i; return *this; }

	Coord::Type & dirComp() { return isVertical(direction) ? y() : x(); }
	Coord::Type dirComp() const { return isVertical(direction) ? y() : x(); }
	Coord::Type & orthComp() { return isVertical(direction) ? x() : y(); }
	Coord::Type orthComp() const { return isVertical(direction) ? x() : y(); }

	int dirStep() const { return isForward(direction) ? +1 : -1; } //!< Outer to Inner
	int orthStep() const { return isVertical(direction) ? dirStep() : -dirStep(); } //!< Start to End

	DirCoord & setToInner(Coord::Type _max) { orthComp() = orthStep() < 0 ? 0 : _max; return *this; }
	DirCoord & setToOuter(Coord::Type _max) { orthComp() = orthStep() > 0 ? 0 : _max; return *this; }
	DirCoord & setToStart(Coord::Type _max) { dirComp() = dirStep() > 0 ? 0 : _max; return *this; }
	DirCoord & setToEnd  (Coord::Type _max) { dirComp() = dirStep() < 0 ? 0 : _max; return *this; }

	//! Sets 'that' in same direction context as 'this'!
	bool moreInnerThan(DirCoord const & _that) const { DirCoord that(_that); that.direction = this->direction; return this->orthStep() > 0 ? this->orthComp() > that.orthComp() : this->orthComp() < that.orthComp(); }
	bool moreOuterThan(DirCoord const & _that) const { DirCoord that(_that); that.direction = this->direction; return this->orthStep() > 0 ? this->orthComp() < that.orthComp() : this->orthComp() > that.orthComp(); }
	bool moreForwardThan(DirCoord const & _that) const { DirCoord that(_that); that.direction = this->direction; return this->dirStep() > 0 ? this->dirComp() > that.dirComp() : this->dirComp() < that.dirComp(); }
	bool moreBackwardThan(DirCoord const & _that) const { DirCoord that(_that); that.direction = this->direction; return this->dirStep() > 0 ? this->dirComp() < that.dirComp() : this->dirComp() > that.dirComp(); }

	Direction direction; // Could be templated for efficiency (see Seam::Coord)
};

std::ostream & operator<<(std::ostream & _stream, DirCoord const & _coord)
{
	_stream << "[" << _coord.direction << ": " << _coord.x() << ", " << _coord.y() << "]";
	return _stream;
}

struct Linf
{
	static int dist(DirCoord const & s, DirCoord const & e)
	{
		return std::max( s.x() > e.x() ? s.x() - e.x() : e.x() - s.x(), s.y() > e.y() ? s.y() - e.y() : e.y() - s.y() );
	}
};

struct L1
{
	static int dist(DirCoord const & s, DirCoord const & e)
	{
		return (s.x() > e.x() ? s.x() - e.x() : e.x() - s.x()) + (s.y() > e.y() ? s.y() - e.y() : e.y() - s.y());
	}
};

struct CircSeam
{
	// TODO: this is a very lazy implementation

	std::vector<DirCoord> m_coords;

	DirCoord & operator[](uint i) { return m_coords[i]; }
	DirCoord operator[](uint i) const { return m_coords[i]; }

	// TODO: incorporate corner-cutting into push_back
	// add member variable m_cost
	// add a while ( L1::dist(_c, m_coords(m_coords.size() - 2)) == 1 ) { <remove corner> }
	void push_back(DirCoord const & _c) { m_coords.push_back(_c); }
	void pop_back() { m_coords.pop_back(); }
	void pop_front() { m_coords.erase(m_coords.begin()); }

	std::vector<DirCoord>::iterator begin() { return m_coords.begin(); }
	std::vector<DirCoord>::iterator end() { return m_coords.end(); }
	std::vector<DirCoord>::const_iterator begin() const { return m_coords.begin(); }
	std::vector<DirCoord>::const_iterator end() const { return m_coords.end(); }

	DirCoord & front() { return m_coords.front(); }
	DirCoord const & front() const { return m_coords.front(); }
	DirCoord & back() { return m_coords.back(); }
	DirCoord const & back() const { return m_coords.back(); }

	size_t size() const { return m_coords.size(); }

	std::vector<DirCoord> & coords() { return m_coords; }

	void cutCorners(double & _totalCost, Table<float> const & _costs)
	{
		for ( int i = 0; i < m_coords.size(); ++i )
		{
			int j = i - 2;
			if ( j < 0 )
				j += m_coords.size();

			//std::cout << "comparing " << i << ":" << thisSeam[i] << " <> " << j << ":" << thisSeam[j] << " (" << L1::dist(thisSeam[i], thisSeam[j]) << ")" << std::endl;

			if ( m_coords[i].direction != m_coords[j].direction )
			{
				if ( Linf::dist(m_coords[i], m_coords[j]) == 1 )
				{
					int k = i - 1;
					if ( k < 0 )
						k += m_coords.size();

					auto ptr = m_coords.begin() + k;

					//std::cout << "REMOVING " << k << ":" << *ptr << std::endl;

					_totalCost -= _costs((*ptr).x(), (*ptr).y());
					m_coords.erase(ptr);
					i = std::max(i - 2, 0) - 1; // correction
				}
			}
		}
	}

	void display()
	{
		for ( auto itr = m_coords.rbegin(); itr != m_coords.rend(); ++itr )
			std::cout << *itr << std::endl;
	}

};

template<typename T>
class Helper
{
public:
	Helper(Table<T> const * _sc) : m_scores(_sc) {}

	bool coordIsValid(DirCoord const & c) const { return c.x() >= 0 && c.y() >= 0 && c.x() < m_scores->width() && c.y() < m_scores->height(); }

	T lookup(DirCoord const & c) const { return (*m_scores)(c.x(), c.y()); }

	Coord::Type dirSize(Direction _d) const { return DirCoord(m_scores->width(), m_scores->height(), _d).dirComp(); }
	Coord::Type orthSize(Direction _d) const { return DirCoord(m_scores->width(), m_scores->height(), _d).orthComp(); }

	DirCoord minCoord(Direction _d) const { return DirCoord(0, 0, _d); }
	DirCoord maxCoord(Direction _d) const { return DirCoord(m_scores->width() - 1, m_scores->height() - 1, _d); }
	Coord::Type maxOrth(Direction _d) const { return maxCoord(_d).orthComp(); }
	Coord::Type maxDir(Direction _d) const { return maxCoord(_d).dirComp(); }

	DirCoord startOutsideCoord(Direction _d) const { return DirCoord(_d).setToStart(maxDir(_d)).setToOuter(maxOrth(_d)); }
	DirCoord startInsideCoord(Direction _d) const { return DirCoord(_d).setToStart(maxDir(_d)).setToInner(maxOrth(_d)); }
	DirCoord endOuterCoord(Direction _d) const { return DirCoord(_d).setToEnd(maxDir(_d)).setToOuter(maxOrth(_d)); }
	DirCoord endInnerCoord(Direction _d) const { return DirCoord(_d).setToEnd(maxDir(_d)).setToInner(maxOrth(_d)); }

	Coord::Type distFromStart(DirCoord const & _c) const
	{
		Coord::Type s = startOutsideCoord(_c.direction).dirComp();
		Coord::Type c = _c.dirComp();
		return c > s ? c - s : s - c;
	}

	Coord::Type distFromOutside(DirCoord const & _c) const
	{
		Coord::Type o = startOutsideCoord(_c.direction).orthComp();
		Coord::Type c = _c.orthComp();
		return c > o ? c - o : o - c;
	}

	bool isWithinItsQuadrant(DirCoord const & _c) const
	{
		return distFromStart(_c) < dirSize(_c.direction) - distFromOutside(_c) && distFromStart(_c) > distFromOutside(_c) && coordIsValid(_c);
	}

	Direction directionOfQuadrantContaining(DirCoord const _c) const
	{
		Direction d = NODIRECTION;
		DirCoord c(_c);
		for ( Direction t : {Downwards, Upwards, Leftwards, Rightwards} )
		{
			c.direction = t;
			if ( isWithinItsQuadrant(c) )
			{
				d = t;
				break;
			}
		}
		return d;
	}

	DirCoord withQuadrantDirection(DirCoord const & _c) const { DirCoord c(_c); c.direction = directionOfQuadrantContaining(c); return c; }

	DirCoord prevIn(DirCoord const & _c) const { return withQuadrantDirection(DirCoord(_c).advBack().advIn()); }
	DirCoord prevOut(DirCoord const & _c) const { return withQuadrantDirection(DirCoord(_c).advBack().advOut()); }
	DirCoord prevStr(DirCoord const & _c) const { return withQuadrantDirection(DirCoord(_c).advBack()); }

	DirCoord bestPrevInStr(DirCoord const & c) const { return lookup(prevIn(c)) < lookup(prevStr(c)) ? prevIn(c) : prevStr(c); }
	DirCoord bestPrevOutStr(DirCoord const & c) const { return lookup(prevOut(c)) < lookup(prevStr(c)) ? prevOut(c) : prevStr(c); }

	DirCoord bestPrevBoth(DirCoord const & c) { return lookup(prevOut(c)) < lookup(bestPrevInStr(c)) ? prevOut(c) : bestPrevInStr(c); }

	DirCoord bestPrev(DirCoord const & _c)
	{
		DirCoord b;
		b = distFromOutside(_c) > 0 ? bestPrevBoth(_c) : bestPrevInStr(_c);
		// Must not go into central pixels
		if ( !directionOfQuadrantContaining(b) )
			b = distFromOutside(_c) > 0 ? bestPrevOutStr(_c) : prevStr(_c);
		return b;
	}

	DirCoord bestInDirectionOf(DirCoord const & src, DirCoord const & tgt, DirCoord const & avoid)
	{
		DirCoord s1 = src;
		DirCoord s2 = src;
		DirCoord s3 = src;
		if ( src.moreForwardThan(tgt) )
		{
			s1.advBack();
			s3 = s1;
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
				s3.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
				s3.advOut();
			}
			else
			{
				s2 = s1;
			}
		}
		else if ( src.moreBackwardThan(tgt) )
		{
			s1.advForward();
			s3 = s1;
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
				s3.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
				s3.advOut();
			}
			else
			{
				s2 = s1;
			}
		}
		else
		{
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
			}
			else
			{
				// src IS tgt
			}
			s1 = s2;
			s3 = s2;
		}

		struct equal { static bool coords(DirCoord const & a, DirCoord const & b){ return a.x() == b.x() && a.y() == b.y(); } };

		T d1 = equal::coords(s1, avoid) ? std::numeric_limits<T>::max() : lookup(s1);
		T d2 = equal::coords(s2, avoid) ? std::numeric_limits<T>::max() : lookup(s2);
		T d3 = equal::coords(s3, avoid) ? std::numeric_limits<T>::max() : lookup(s3);

		//std::cout << "1. " << s1 << ":" << d1 << std::endl;
		//std::cout << "2. " << s2 << ":" << d2 << std::endl;
		//std::cout << "3. " << s3 << ":" << d3 << std::endl;

		DirCoord best = src;

		if ( std::min(std::min(d1, d2), d3) < std::numeric_limits<T>::max() ) // if solution exists, find it...
			best = d1 < d2 ? (d1 < d3 ? s1 : (d3 < d2 ? s3 : s2)) : (d2 < d3 ? s2 : (d3 < d1 ? s3 : s1));

		return best;
	}

	void display(CircSeam const & cs)
	{
		display(cs, 0, m_scores.width() - 1, 0, m_scores->height() - 1);
	}

	void display(CircSeam const & cs, uint si, uint ei, uint sj, uint ej)
	{
		float total = 0;
		for ( auto c : cs )
			total += ( (*m_scores)(c.x(), c.y()) );
		printf("CircSeam: length = %u; cost = %f;\n", cs.size(), total);
		for ( uint j = sj; j <= ej; ++j )
		{
			for ( uint i = si; i <= ei; ++i )
			{
				bool found = false;
				for ( auto c : cs )
					found |= ( c.x() == i && c.y() == j );
				printf("%8.0f%s ", (*m_scores)(i,j), found ? "*" : " ");
			}
			printf("\n");
		}
	}

private:
	Table<T> const * m_scores;
};

//! Minimum cost seam all the way around, not using any entries where _mask < 255.f
Table<float> findMin(Table<float> const & _costs, Table<float> const & _mask)
{
	// Note: this is an approximation - the loop closure uses heuristics to get an answer rather than anything that can guarantee optimality
	// Note: current version necessarily includes the center pixel(s) in the output mask

	// TODO: use /mincut/ instead

	// increase size of the working table by 1 in every direction
	// to allow for prohibited areas going to the edge
	Table<float> costs(_costs.width() + 2, _costs.height() + 2, pow(2.f, 17.f)); // 256^2 * 2
	Table<double> cumulative(costs.width(), costs.height(), 0.f);
	Table<DirCoord> refs(costs.width(), costs.height(), DirCoord(Downwards));

	// TODO: should be some sort of global variable
	float epsilon = 0.1f;

	// fill costs
	for ( uint j = 0; j < _costs.height(); ++j )
		for ( uint i = 0; i < _costs.width(); ++i )
			costs(i + 1, j + 1) = _mask(i, j) < 255.f - epsilon ? std::numeric_limits<float>::infinity() : _costs(i, j);

	Helper<double> helper(&cumulative);

	// hack - algorithm fails if <inf>s don't go to center in at least one quadrant
	DirCoord hack;
	DirCoord center(cumulative.width() / 2, cumulative.height() / 2, NODIRECTION);
	for ( uint j = 1; j < costs.height() - 1; ++j )
		for ( uint i = 1; i < costs.width() - 1; ++i )
		{
			if ( costs(i, j) == std::numeric_limits<float>::infinity() )
			{
				DirCoord candidate(i, j, NODIRECTION);
				if ( L1::dist(candidate, center) < L1::dist(hack, center) )
					hack = candidate;
			}
			else
			{
				if ( !helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) ||
					 helper.directionOfQuadrantContaining( DirCoord(i + 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j + 1, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i + 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j + 1, NODIRECTION) ) ||
					 helper.directionOfQuadrantContaining( DirCoord(i + 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j - 1, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i + 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j - 1, NODIRECTION) ) ||
					 helper.directionOfQuadrantContaining( DirCoord(i - 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j + 1, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i - 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j + 1, NODIRECTION) ) ||
					 helper.directionOfQuadrantContaining( DirCoord(i - 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j - 1, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i - 1, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) && helper.directionOfQuadrantContaining( DirCoord(i, j, NODIRECTION) ) != helper.directionOfQuadrantContaining( DirCoord(i, j - 1, NODIRECTION) )
					)
				{
					costs(i, j) = std::numeric_limits<float>::infinity();
				}
			}
		}
	// create barrier
	hack.direction = helper.directionOfQuadrantContaining(hack);
	// std::cout << "HACK = " << hack <<std::endl;
	if ( hack.direction )
	{
		bool lastWasInwards = true;
		while ( helper.isWithinItsQuadrant(hack) )
		{
			if ( lastWasInwards && hack.moreBackwardThan(center) )
			{
				hack.advForward();
				lastWasInwards = false;
			}
			else if ( lastWasInwards && hack.moreForwardThan(center) )
			{
				hack.advBack();
				lastWasInwards = false;
			}
			else
			{
				hack.advIn();
				lastWasInwards = true;
			}

			costs(hack.x(), hack.y()) = std::numeric_limits<float>::infinity();
		}
	}

	// order of quadrants
	std::vector<Direction> dirs = { Downwards, Rightwards, Upwards, Leftwards };

	// initialize cumulatives
	for ( DirCoord i = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(i); i.advBack().advIn() )
		cumulative(i.x(), i.y()) = costs(i.x(), i.y()) == std::numeric_limits<float>::infinity() ? std::numeric_limits<float>::infinity() : 0;

	double bestCost = 0;
	CircSeam bestSeam;

	// dynamic programming
	for ( Direction d : dirs )
	{
		for ( DirCoord base = helper.startOutsideCoord(d); !base.moreForwardThan(helper.endOuterCoord(d)); base.advForward() )
		{
			// create a default 'around the outside' seam
			bestSeam.push_back(base);
			bestCost += costs(base.x(), base.y());

			for ( DirCoord curr = base; helper.isWithinItsQuadrant(curr); curr.advIn() )
			{
				DirCoord b = helper.bestPrev(curr);
				cumulative(curr.x(), curr.y()) = costs(curr.x(), curr.y()) + cumulative(b.x(), b.y());
				refs(curr.x(), curr.y()) = b;
			}
		}
	}

	// **** Backtracking *****

	// create a list the end points in order of total seam cost
	typedef std::pair<double, DirCoord> DDC;
	struct compDDC { bool operator()(DDC const & a, DDC const & b) { return a.first > b.first; } };
	std::vector<DDC> heap;

	for ( DirCoord e = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(e); e.advBack().advIn() )
	{
		double v = cumulative(e.x(), e.y());
		if ( v < std::numeric_limits<double>::infinity() )
		{
			heap.push_back(std::make_pair(v, e));
			std::push_heap(heap.begin(), heap.end(), compDDC());
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compDDC());

	Helper<float> costsHelper(&costs); // DEBUG

	// go through seams (by end point) in order of cost (lowest first)
	// Note: heap in reverse order
	for ( auto seamEndPtr = heap.rbegin(); seamEndPtr != heap.rend(); ++seamEndPtr )
	{
		//std::cout << "--" << std::endl;

		CircSeam thisSeam;
		bool leftLastQuadrant = false;
		double thisCost = (*seamEndPtr).first;

		// actual backtracking
		DirCoord curr = (*seamEndPtr).second;
		DirCoord test(curr.x(), curr.y(), dirs.back());
		while ( leftLastQuadrant != helper.isWithinItsQuadrant(test) )
		{
			thisSeam.push_back(curr);
			curr = refs(curr.x(), curr.y());
			test.x() = curr.x();
			test.y() = curr.y();
			if ( !leftLastQuadrant && !helper.isWithinItsQuadrant(test) )
				leftLastQuadrant = true;
		}

		// cut corners on quadrant transitions
		// (because of isWithinItsQuadrant(b) during dynamic programming)
		thisSeam.cutCorners(thisCost, costs);

		if (  thisCost < bestCost )
		{
			if ( Linf::dist(thisSeam.front(), thisSeam.back()) == 1 )
			{
				// done
				bestCost = thisCost;
				bestSeam = thisSeam;
			}
			else
			{
				//printf("This needs work...\n");

				// try connecting from back (greedy approach)
				while ( Linf::dist(thisSeam.front(), thisSeam.back()) > 1 )
				{
					DirCoord s = thisSeam.front();
					DirCoord e = thisSeam.back();
					DirCoord avoid = thisSeam[thisSeam.size() - 2]; // seam must be bigger than 1

					//std::cout << "Advancing " << e << " towards " << s << " avoiding " << avoid << std::endl;

					DirCoord u = costsHelper.bestInDirectionOf(e, s, avoid);

					if ( Linf::dist(u, e) == 0 ) // not going anywhere
						continue; // not joined up so on to the next seam

					// corner to be cut
					while ( thisSeam.size() > 2 && Linf::dist(u, thisSeam[1]) < 2)
					{
						//std::cout << "Removing start: " << thisSeam.front() << std::endl;

						s = thisSeam.front();
						thisCost -= costs(s.x(), s.y());
						thisSeam.pop_front();
					}

					// corner to be cut
					while ( thisSeam.size() > 1 && Linf::dist(u, thisSeam[thisSeam.size() - 2]) == 1 )
					{
						//std::cout << "Removing end: " << thisSeam.back() << std::endl;

						e = thisSeam.back();
						thisCost -= costs(e.x(), e.y());
						thisSeam.pop_back();

						// correct previous change to direction
						u.direction = thisSeam.back().direction;
					}

					// check for a non-continuation of direction
					if ( thisSeam.back().moreOuterThan(u) && thisSeam.size() > 1 )
						if ( thisSeam[thisSeam.size() - 2].moreBackwardThan(thisSeam.back()) )
						{
							//std::cout << thisSeam[thisSeam.size() - 2] << " -> " << thisSeam.back() << " is forwards" << std::endl;
							if ( !thisSeam.back().moreBackwardThan(u) )
							{
								//std::cout << u << " is not continuation" << std::endl;
								u.direction = s.direction;
							}
						}
						else
						{
							//std::cout << thisSeam[thisSeam.size() - 2] << " -> " << thisSeam.back() << " is backwards" << std::endl;
							if ( !thisSeam.back().moreForwardThan(u) )
							{
								//std::cout << u << " is not continuation" << std::endl;
								u.direction = s.direction;
							}
						}

					//std::cout << " ==> " << u << std::endl;

					thisSeam.push_back( u );
					thisCost += costs(u.x(), u.y());

					//std::cout << "end gap dist = " << L1::dist(thisSeam.front(), thisSeam.back()) << " & cost = " << thisCost << std::endl;
				}

				thisSeam.cutCorners(thisCost, costs);

				if ( thisCost < bestCost )
				{
					bestCost = thisCost;
					bestSeam = thisSeam;
				}

				// try other things?
			}

		}

	}

	// create mask
	Table<float> output(_costs.width(), _costs.height(), 1.f);

	struct l
	{
		float dummy; // TODO! tidy up this laziness
		Table<float> * ref;

		float & operator()(DirCoord const & c)
		{
			//std::cout<<c<<": "<<(check(c)?"set":"ignored")<<std::endl;
			return check(c) ? (*ref)(c.x() - 1, c.y() - 1) : dummy;
		}

		bool check(DirCoord const & c) { return c.x() > 0 && c.y() > 0 && c.x() <= ref->width() && c.y() <= ref->height(); }

		void fillToCorner(CircSeam const & _seam, uint _curr_i, uint _comp_i, Helper<float> const & _helper)
		{
			// fill to corner
			DirCoord curr = _seam[_curr_i];
			DirCoord comp = _seam[_comp_i];
			int i_dir = _comp_i - _curr_i;

			//std::cout<<"fill "<<(curr.moreForwardThan(comp)?"forward":"back")<<" from "<<curr<<std::endl;

			while ( !(curr.moreForwardThan(comp) || curr.moreBackwardThan(comp)) )
			{
				_comp_i += i_dir;
				comp = _seam[_comp_i];
			}

			curr.setToOuter(_helper.maxOrth(curr.direction));
			if ( curr.moreForwardThan(comp) )
				for ( curr.advForward(); _helper.isWithinItsQuadrant(curr); curr.advForward() )
					for ( DirCoord cc = curr; _helper.isWithinItsQuadrant(cc); cc.advIn() )
						(*this)(cc) = 0;
			else
				for ( curr.advBack(); _helper.isWithinItsQuadrant(curr); curr.advBack() )
					for ( DirCoord cc = curr; _helper.isWithinItsQuadrant(cc); cc.advIn() )
						(*this)(cc) = 0;
		}

	} output_;
	output_.ref = &output;

	for ( int i = 0; i < bestSeam.size(); ++i )
	{
		int j = i > 0 ? i - 1 : bestSeam.size() - 1;
		int k = i < bestSeam.size() - 1 ? i + 1 : 0;

		DirCoord curr = bestSeam[i];
		DirCoord prev = bestSeam[j];
		DirCoord next = bestSeam[k];

		if ( curr.direction != prev.direction )
			output_.fillToCorner(bestSeam, i, k, costsHelper);
		if ( curr.direction != next.direction )
			output_.fillToCorner(bestSeam, i, j, costsHelper);

		if ( output_.check(curr) )
		{
			output_(curr) = 0.5f;
			if ( curr.direction != prev.direction || curr.moreForwardThan(prev) || curr.moreBackwardThan(prev) )
				for ( curr.advOut(); output_.check(curr); curr.advOut() ) // fill out to edge
					output_(curr) = 0.f;
		}
	}

	// ensure that all <float>::inf pixels are 1.f in mask
	for ( uint j = 0; j < output.height(); ++j )
		for ( uint i = 0; i < output.width(); ++i )
			if ( costs(i + 1, j + 1) == std::numeric_limits<float>::infinity() )
				output(i, j) = 1.f;

	// DEBUG
	if ( output(0,0) == 1 )
	{
		printf("BAD MASK\n");
		costsHelper.display(bestSeam, 0, 16, 0, 16);
		printf(">>>\n");
		output.display(0, 16, 0, 16);
		printf("###\n");
	}

	return output * 255.f;
}

} // end namespace CircSeams

} // namespace TexSynth
