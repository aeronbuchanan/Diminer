/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of Diminer, an digital inpainting resource.
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

#include <iostream>

#include "common.h"
#include "inpainters.h"


int main(int argc, char** argv)
{
	std::cout << "Diminer: Digital Image Inpainting Resources by Aeron Buchanan" << std::endl;

	cimg_usage("Usage: Diminer [options] -i input\n");

	char const * imageFilename = cimg_option("-i", (char*)0, "image file to be inpainted");
	char const * maskFilename = cimg_option("-m", (char*)0, "mask image denoting region to be inpainted (values > 127)");
	char const * outputFilename = cimg_option("-o", (char*)0, "output filename");
	uint inpaintingFunc = cimg_option("-f", 0, "0 = bleed; 1 = weighted; 2 = gradient-weighted");

	if ( !imageFilename )
	{
		std::exit(0);
	}

	uint defaultInpaintingFunc = 0;
	if ( inpaintingFunc > 2 )
	{
		std::cout << "Inpainting function index (" << inpaintingFunc << ") unavailable - switching to default (" << defaultInpaintingFunc << ")" << std::endl;
		inpaintingFunc = defaultInpaintingFunc;
	}

	CImg<uchar> image(imageFilename);
	CImg<uchar> mask;
	if ( maskFilename)
	{
		mask.load(maskFilename);
	}
	else
	{
		mask.assign(image.width(),image.height(),1,1);
		cimg_forXY(image,x,y)
		{
			Color c;
			c.r = image(x,y,0,0);
			c.g = image(x,y,0,1);
			c.b = image(x,y,0,2);
			mask(x,y) = imgMaskTest(c) ? 255 : 0;
		}
	}

	CImg<uchar> regions(image.width(), image.height(), 1, 1, 0); // TODO: float matrix

	CImg<uchar> * maskOrig;

	if ( inpaintingFunc == 2 )
	{
		// TODO: this is a hack to overcome gradient being taken without respect to boundary
		maskOrig = new CImg<uchar>(regions);
		cimg_forXY(mask,x,y)
			(*maskOrig)(x,y) = maskTest(mask(x,y));

		uint r = 2;
		mask.dilate(r * 2 + 1);
	}

	//CImgDisplay main_disp(mask, "Original Mask");

	//*** detect connected regions ***//

	// fill regions
	FillHelper filler(&mask, &regions);

	int count = 0;
	cimg_forXY(mask,x,y)
	{
		if ( mask(x,y) == 255 && regions(x,y) == 0 )
		{
			++count;
			uchar col = count; // 255 / count;
			filler.fill(Coord(x,y), col);
		}
	}

	std::cout << "Found " << count << " regions." << std::endl;

	//CImgDisplay regions_disp(regions,"Regions");

	// Find region 4-boundaries
	std::vector<BoundaryColors> boundaries;
	int W = regions.width() - 1;
	int H = regions.height() - 1;
	for ( int i = 0; i < count; ++i )
	{
		BoundaryColors cs;
		int count = i + 1;
		uchar col = count; // 255 / count;
		cimg_forXY(regions,x,y)
		{
			if ( regions(x,y) == 0 )
			{
				if (
						( y > 0 && regions(x, y - 1) == col ) ||
						( y < H && regions(x, y + 1) == col ) ||
						( x > 0 && regions(x - 1, y) == col ) ||
						( x < W && regions(x + 1, y) == col )
					)
				{
					Color color;
					color.r = image(x,y,0,0);
					color.g = image(x,y,0,1);
					color.b = image(x,y,0,2);
					cs.push_back(CoordCol(Coord(x, y), color));
				}
			}
		}
		boundaries.push_back(cs);
	}

	std::cout << "Calculated boundaries." << std::endl;

#if 1 // DEBUG
		CImgDisplay debug_disp(image, "Debug");
		CImgDisplay maskorig_disp(mask,"Mask");
		AnimDisp debugDisp(&debug_disp);
#endif

	// Inpainting
	std::vector<Inpainter*> inpainters;
	for ( uint i = 0; i < boundaries.size(); ++i )
	{
		switch ( inpaintingFunc )
		{
		case 0:
			inpainters.push_back(new BleedInpainter(&boundaries[i]));
			break;
		case 1:
			inpainters.push_back(new WeightedInpainter(&boundaries[i]));
			break;
		case 2:
			inpainters.push_back(new GradientWeightedInpainter(&boundaries[i], &image, maskOrig));
			break;
		default:
			inpainters.push_back(new BleedInpainter(&boundaries[i]));
		}
	}

	std::cout << "Created inpainters." << std::endl;

	float i = 0;
	float i_total = regions.width() * regions.height();
	std::cout << "Inpainting: ";

	// TODO: allow 'blend mode' where mask value is a blend coefficient

	cimg_forXY(regions,x,y)
	{
		if ( regions(x,y) > 0 )
		{
			// process associated boundary
			uint r = regions(x,y) - 1;
			// TODO: check r is within bounds
			Color c = inpainters[r]->pixelColor(Coord(x, y));
			image(x,y,0,0) = c.r;
			image(x,y,0,1) = c.g;
			image(x,y,0,2) = c.b;
		}

		printf("% 6.0f%% ", 100 * i++ / (i_total - 1));
		for ( uint j = 0; j < 8; ++j ) printf("%c", 8);
	}

	std::cout << "complete." << std::endl;

	CImgDisplay output_disp(image,"Inpainted Original");

	if ( outputFilename )
	{
		image.save(outputFilename);
		mask.save("mask_with_boundary.png");
	}

	while ( !output_disp.is_closed() )
	{
		output_disp.wait(1000);
	}

	return 0;
}

