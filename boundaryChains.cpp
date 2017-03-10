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

// TODO: finish DEBUG leveling
#define DEBUGLEVEL 1
#define DEBUG1(output) {if(0) { output }}
#define DEBUG2(output) {if(0) { output }}
#define DEBUG3(output) {if(0) { output }}

#include "boundaryChains.h"

using namespace Diminer;

/*
 * Groups coords together by 8-neighbour connectivity
 *
 * Invariant:= new coords will either: 
 * 	connect to two adjacent chain coords -> insert between
 *      connect to both ends (when ends distinct) -> close loop
 *      connect to one end only -> become end
 *      connect to one non-end only -> ignore (too complex a boundary)
 *      not connect -> ignore
 *
 * NB: This code can ONLY cope with simple boundaries around image regions
 * that are fully 4-neighbour connected and have no holes or thin protrusions.
 * Complications must be simplified from the image before processing.
 */

uint ChainLink::COUNT = 0;

NeighbourType ChainLink::neighbourTypeOf(ChainLinkPtr cc)
{
	int dx = abs(cc->x() - m_c->x);
	int dy = abs(cc->y() - m_c->y);
	NeighbourType r = NeighbourType::NONE;
	if ( (dx == 1 && dy == 0) || (dx == 0 && dy == 1) ) r = NeighbourType::FOUR;
	else if ( dx == 1 && dy == 1 ) r = NeighbourType::EIGHT;
	return r;
}

ChainLinkPtr ChainLink::link(ChainLinkPtr cc)
{
	if ( ! hither ) return hither = cc;
	else if ( ! thither ) return thither = cc;
	return ChainLinkPtr();
}

ChainLinkPtr ChainLink::next(ChainLinkPtr prev)
{
	DEBUG3(std::cout << "Looking for next in sequence for ID" << m_id << " away from ";
	if ( prev ) std::cout << "ID" << prev->m_id << std::endl;
	else std::cout << "end" << std::endl;)

	ChainLinkPtr ptr = ChainLinkPtr();

	if ( hither == prev ) ptr = thither;
	else if ( thither == prev ) ptr = hither;

	DEBUG3(if (ptr) std::cout << "Got: ID" << ptr->m_id << std::endl;
	else std::cout << "End of the line" << std::endl;)

	return ptr;
}

bool ChainLink::replaceLink(ChainLinkPtr old, ChainLinkPtr cc)
{
	if ( hither == old ) { hither = cc; return true; }
	if ( thither == old ) { thither = cc; return true; }
	return false;
}

ChainGrouping::ChainGrouping(ChainLinkPtr cc)
{
	m_isClosedLoop = false;

	m_xmax = cc->x() + 1;
	m_xmin = cc->x() - 1;
	m_ymax = cc->y() + 1;
	m_ymin = cc->y() - 1;

	m_chainLinks.push_back(cc);
	m_chainStart = cc;
	m_chainEnd = cc;
}

bool ChainGrouping::couldInclude(ChainLinkPtr cc)
{
	return !m_isClosedLoop && cc->x() >= m_xmin && cc->x() <= m_xmax && cc->y() >= m_ymin && cc->y() <= m_ymax;
}

bool ChainGrouping::matchedToStart(ChainLinkPtr cc)
{
	return !m_isClosedLoop && matchedToExtremity(cc, m_chainStart);
}

bool ChainGrouping::matchedToEnd(ChainLinkPtr cc)
{
	return !m_isClosedLoop && matchedToExtremity(cc, m_chainEnd);
}

bool ChainGrouping::canBeALoop()
{
	return m_chainLinks.size() > 2;
}

void ChainGrouping::addToExtremity(ChainLinkPtr cc, ChainLinkPtr e)
{
	addLink(cc);
	cc->link(e);

	if ( e == m_chainStart ) m_chainStart = cc;
	else if ( e == m_chainEnd ) m_chainEnd = cc;

	e->link(cc);
}

void ChainGrouping::closeLoop(ChainLinkPtr cc)
{
	addLink(cc);
	cc->link(m_chainStart);
	cc->link(m_chainEnd);
	
	m_chainStart->link(cc);
	m_chainEnd->link(cc);

	m_chainStart = cc;
	m_chainEnd = cc;

	m_isClosedLoop = true;
}

void ChainGrouping::mergeWith(ChainGrouping & other, ChainLinkPtr thisEnd, ChainLinkPtr otherEnd)
{
	/* concatenate other chain onto end of this one
	 * joining together 'thisEnd' of this chain to the 'otherEnd' of the other chain
	 */

	// expand bounding box as necessary
	if ( other.xmax() > m_xmax ) m_xmax = other.xmax();
	if ( other.xmin() < m_xmin ) m_xmin = other.xmin();
	if ( other.xmax() > m_ymax ) m_ymax = other.ymax();
	if ( other.xmin() < m_ymin ) m_ymin = other.ymin();

	// copy across coords
	m_chainLinks.insert(m_chainLinks.end(), other.chainLinks()->begin(), other.chainLinks()->end());

	// link ends
	otherEnd->link(thisEnd);
	thisEnd->link(otherEnd);

	// update end ptrs
	if ( m_chainStart == thisEnd )
		m_chainStart = other.chainStart() == otherEnd ? other.chainEnd() : other.chainStart();
	else
		m_chainEnd = other.chainStart() == otherEnd ? other.chainEnd() : other.chainStart();
}

bool ChainGrouping::addCoordToMiddle(ChainLinkPtr cc)
{
	// no complex checks here - assume that the boundary candidate identification process is sound

	DEBUG2(std::cout << "Starting 'corner' search for insertion..." << std::endl;)

	bool success = false;
	if ( couldInclude(cc) )
	{
		DEBUG2(std::cout << "Bounding box inclusion check passed; continuing..." << std::endl;)

		ChainLinkPtr prev = ChainLinkPtr();
		DEBUG2(std::cout << "Ready 1" << std::endl;)
		ChainLinkPtr curr = m_chainStart;
		DEBUG2(std::cout << "Ready 2" << std::endl;)

		DEBUG2(if ( curr == m_chainStart ) std::cout << "We are at the start" << std::endl;)
		DEBUG2(if ( curr != m_chainEnd ) std::cout << "We are not at the end" << std::endl;)
		DEBUG2(std::cout << "Current ID: " << curr->m_id << std::endl;)
		DEBUG2(if ( prev ) std::cout << "Previous ID: " << prev->m_id << std::endl;)
		DEBUG2(std::cout << "Ready 3" << std::endl;)

		while ( curr && curr != m_chainEnd )
		{
			DEBUG2(std::cout << "Comparing to ID" << curr->m_id << " (prev: ID";)
			DEBUG2(if ( prev ) std::cout << prev->m_id; else std::cout << "---";)

			ChainLinkPtr next = curr->next(prev);

			DEBUG2(std::cout << "; next: ID" << next->m_id << ")" << std::endl;)

			// TODO: flag near misses
			if ( 
				cc->neighbourTypeOf(curr) != NeighbourType::NONE &&
				cc->neighbourTypeOf(next) != NeighbourType::NONE
			) {
				addLink(cc);
				curr->replaceLink(next, cc);
				next->replaceLink(curr, cc);
				cc->link(curr);
				cc->link(next);
				success = true;
				break;
			}

			prev = curr;
			curr = next;
		}
	}

	DEBUG2(if ( success ) std::cout << "Added to middle of chain" << std::endl;)

	return success;
}

void ChainGrouping::addLink(ChainLinkPtr cc)
{
	m_chainLinks.push_back(cc);
	if ( cc->x() <= m_xmin ) m_xmin = cc->x() - 1;
	if ( cc->x() >= m_xmax ) m_xmax = cc->x() + 1;
	if ( cc->y() <= m_ymin ) m_ymin = cc->y() - 1;
	if ( cc->y() >= m_ymax ) m_ymax = cc->y() + 1;
}

bool ChainGrouping::matchedToExtremity(ChainLinkPtr cc, ChainLinkPtr e)
{
	// neighbour of ONLY the end coord
	bool r = e->neighbourTypeOf(cc) != NeighbourType::NONE;
	// if neighbour of next coord in chain, still OK if just diagonally adjacent to that coord
	/*
	 *  #    #   O#O  # O  #O
	 *  #   O#O   #    #O  O#
	 * OOO            OOO
	 *
         * yes  yes  no!  yes  no!
	 */
	if ( r && (e->hither || e->thither) )
	{
		ChainLinkPtr o = e->hither ? e->hither : e->thither;
		r = o->neighbourTypeOf(cc) != NeighbourType::FOUR;
	}
	return r;
}

ChainManager::ChainManager()
{
	DEBUG1(std::cout << "[" << std::endl;)
}

ChainManager::~ChainManager()
{
	DEBUG1(std::cout << "[] ]" << std::endl;)
}

void ChainManager::add(CoordPtr c)
{
	ChainLinkPtr cc = std::make_shared<ChainLink>(c);
	m_chainLinks.push_back(cc);

	/* each coord will:
	 *	not connect -> new chain
	 *	connect to one end each of two chains -> concatenate chains
	 *	connect to one chain -> add to chain
	 *
	 * NB: complex boundaries can result in incorrect boundary shapes
	 */
	uint endCount = 0;
	ChainGrouping * other;

	for ( uint j = 0; j < m_chainGroupings.size(); ++j )
	{
		DEBUG2(std::cout << "Checking ends of group " << j << "..." << std::endl;)
		bool s = m_chainGroupings[j].matchedToStart(cc);
		bool e = m_chainGroupings[j].matchedToEnd(cc);

		DEBUG2(std::cout << "Checking loop..." << std::endl;)
		if ( s && e && m_chainGroupings[j].canBeALoop() )
		{
			// close loop
			DEBUG2(std::cout << "Closing loop" << std::endl;)
			m_chainGroupings[j].closeLoop(cc);
			endCount = 2;
			break;
		}
		else if ( s || e )
		{
			DEBUG2(std::cout << "Matched to end..." << std::endl;)
			ChainLinkPtr ee = e ? m_chainGroupings[j].chainEnd() : m_chainGroupings[j].chainStart();
			if ( endCount == 0 )
			{
				DEBUG2(std::cout << "Adding to end" << std::endl;)
				m_chainGroupings[j].addToExtremity(cc, ee);
				other = &m_chainGroupings[j];
				endCount = 1;
			}
			else
			{
				DEBUG2(std::cout << "Merging!" << std::endl;)
				other->mergeWith(m_chainGroupings[j], cc, ee);
				m_chainGroupings.erase(m_chainGroupings.begin() + j);
				break;
			}
		}
	}

	uint addCount = 0;

	if ( endCount == 0 )
	{
		DEBUG2(std::cout << "Attempting insertion..." << std::endl;)
		// try to add to middle of chains
		for ( uint j = 0; j < m_chainGroupings.size(); ++j )
		{
			DEBUG2(std::cout << "Perhaps into group " << j << "..." << std::endl;)
			if ( m_chainGroupings[j].addCoordToMiddle(cc) )
			{
				DEBUG2(std::cout << "Inserted!" << std::endl;)
				++addCount;
				break;
			}
		}
	}

	DEBUG2(std::cout << "Final checks..." << std::endl;)
	if ( endCount == 0 && addCount == 0 )
	{
		DEBUG2(std::cout << "Creating new group" << std::endl;)
		m_chainGroupings.push_back(ChainGrouping(cc));
	}
	DEBUG2(std::cout << "Finished" << std::endl;)

	// DEBUG
DEBUG1(
	std::cout << "  [" << std::endl;
	for ( uint i = 0; i < m_chainGroupings.size(); ++i )
	{
		std::cout << "    [" << std::endl;
		ChainLinkRefs cs = *m_chainGroupings[i].chainLinks();
		for ( uint j = 0; j < cs.size(); ++j )
		{
			ChainLink c = *cs[j];
			std::cout << "      {ptr: " << cs[j] << ", id: " << c.m_id << ", x: " << c.x() << ", y: " << c.y() << ", hither: ";
			if ( c.hither ) std::cout << c.hither->m_id;
			else std::cout << "null";
			std::cout << ", thither: "; 
			if ( c.thither ) std::cout << c.thither->m_id;
			else std::cout << "null";
			std::cout << "}" << (j < cs.size() - 1 ? "," : "") << std::endl;
		}
		std::cout << "    ]" << (i < m_chainGroupings.size() - 1 ? "," : "") << std::endl;
	}
	std::cout << "  ]," << std::endl;
)
	// DEBUG END
}

// TODO: deal with masks with multiple boundaries
Coords ChainManager::orderedChains()
{
	Coords cs;

	DEBUG2(std::cout << "DEBUG: following the chain" << std::endl;)

	if ( m_chainGroupings.size() > 0 )
	{
		ChainGrouping cg = m_chainGroupings[0];
		ChainLinkPtr curr = cg.chainStart();
		ChainLinkPtr prev = cg.isClosedLoop() ? curr->thither : ChainLinkPtr();
		while ( curr )
		{
			cs.push_back(curr->m_c);	
			ChainLinkPtr next = curr->next(prev);
			prev = curr;
			curr = next;
		}
	}

	return cs;
}

bool ChainManager::isGood(int widthRef, int heightRef)
{
	if ( m_chainGroupings.size() == 0 ) return false;

	ChainGrouping cg = m_chainGroupings[0];

// TODO: allow boundary coords to include edge pixels?
#define ON_BOUNDARY(C, W, H) (C->x() <= 1 || C->x() >= W - 1 || C->y() <= 1 || C->y() >= H - 1)
	return  m_chainGroupings[0].isClosedLoop() || (ON_BOUNDARY(cg.chainStart(), widthRef, heightRef) && ON_BOUNDARY(cg.chainEnd(), widthRef, heightRef));
#undef ON_BOUNDARY
}



