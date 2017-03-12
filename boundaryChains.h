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

#include <stdlib.h>
#include <memory>

#include "diminer.h"

namespace Diminer
{

enum class NeighbourType {NONE, FOUR, EIGHT};

class ChainLink;
typedef std::shared_ptr<ChainLink> ChainLinkPtr;

class ChainLink
{
public:
	ChainLink(CoordPtr c) : m_c(c), hither(ChainLinkPtr()), thither(ChainLinkPtr()) { m_id = COUNT++; }

	int x() { return m_c->x; }
	int y() { return m_c->y; }

	NeighbourType neighbourTypeOf(ChainLinkPtr);
	ChainLinkPtr link(ChainLinkPtr);
	ChainLinkPtr next(ChainLinkPtr prev);
	bool replaceLink(ChainLinkPtr old, ChainLinkPtr cc);

	CoordPtr m_c;
	ChainLinkPtr hither;
	ChainLinkPtr thither;

	int m_id;
	static uint COUNT;
};

typedef std::vector<ChainLink> ChainLinks;
typedef std::vector<ChainLinkPtr> ChainLinkRefs;

class ChainGrouping
{
public:
	ChainGrouping(ChainLinkPtr);
	bool couldInclude(ChainLinkPtr);
	int addToExtremity(ChainLinkPtr, bool, bool);
	bool addCoordToMiddle(ChainLinkPtr);
	void mergeWith(ChainGrouping & other);

	ChainLinkRefs * chainLinks() { return &m_chainLinks; }
	ChainLinkPtr chainStart() { return m_chainStart; }
	ChainLinkPtr chainEnd() { return m_chainEnd; }
	bool isClosedLoop() { return m_isClosedLoop; }
	int xmax() { return m_xmax; }
	int xmin() { return m_xmin; }
	int ymax() { return m_ymax; }
	int ymin() { return m_ymin; }

private:
	int m_xmax, m_ymax, m_xmin, m_ymin;
	ChainLinkPtr m_chainStart;
	ChainLinkPtr m_chainEnd;
	ChainLinkRefs m_chainLinks;
	bool m_isClosedLoop;

	void addLink(ChainLinkPtr);
	bool canBeALoop();
	void closeLoop(ChainLinkPtr);
	bool matchExtremity(ChainLinkPtr, ChainLinkPtr);
	bool matchExtremityButOne(ChainLinkPtr, ChainLinkPtr);
};

class ChainManager
{
public:
	ChainManager();
	~ChainManager();

	void add(CoordPtr c);
	Coords orderedChains();
	bool isGood(int widthRef, int heightRef);
	uint coordCount() { return m_coordCount; }

private:
	std::vector<ChainGrouping> m_chainGroupings;
	ChainLinkRefs m_chainLinks;

	uint m_coordCount;
};


} // end namespace Diminer



