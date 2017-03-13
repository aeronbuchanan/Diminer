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

#include <iostream>
#include <stdlib.h>

#include "diminer.h"
#include "patch.h"
#include "inpainters.h"
#include "boundaryChains.h"

using namespace Diminer;

int main(int argc, char** argv)
{
	std::cout << "Diminer: Digital Image Inpainting Resources by Aeron Buchanan" << std::endl;

	cimg_usage("Usage: Diminer [options] -i input\n");

	char const * imageFilename = cimg_option("-i", "image.jpg", "image file to be inpainted");
	char const * maskFilename = cimg_option("-m", "", "mask image denoting region to be inpainted (values > 127)");
	char const * outputFilename = cimg_option("-o", "inpainted.jpg", "output filename");
	int inpaintingFunc = cimg_option("-f", 2, "0 = bleed; 1 = weighted; 2 = gradient-weighted");
	float jitter = cimg_option("-j", 0.35, "jitter between 0.f and 1.f");
	float smoothness = cimg_option("-s", 2.f, "smoothness");
	int dilation = cimg_option("-r", 0, "mask dilation radius");
	bool display = cimg_option("-D", false, "display");
        int limitInpainting = cimg_option("-L", 0, "limit inpainting to max extent of image area: 0 = full image; 1 = limited; 2 = limit and crop");

	if ( !imageFilename )
	{
		std::exit(0);
	}

	int defaultInpaintingFunc = 0;
	if ( inpaintingFunc > 2 )
	{
		std::cout << "Inpainting function index (" << inpaintingFunc << ") unavailable - switching to default (" << defaultInpaintingFunc << ")" << std::endl;
		inpaintingFunc = defaultInpaintingFunc;
	}

	CImg<uchar> image(imageFilename);
	CImg<uchar> mask(image.width(), image.height(), 1, 1, 0);
	if (maskFilename)
	{
		CImg<uchar> maskFile(maskFilename);

		cimg_forXY(mask,x,y)
		{
                        bool markedAsMask = maskFile(x,y,0,0) > 0 || maskFile(x,y,0,1) > 0 || maskFile(x,y,0,2) > 0;
			mask(x,y) = markedAsMask ? IS_MASKED : IS_NOT_MASKED; // TODO: use greyscale information if available
                        //std::cout << "(" << x << ", " << y << "): [" << int(mask(x,y,0,0)) << "; " << int(mask(x,y,0,1)) << "; " << int(mask(x,y,0,2)) << "] => " << int(mask(x,y)) << std::endl;
		}
	}
	else
	{
		cimg_forXY(image,x,y)
		{
			Color c;
			c.r = image(x,y,0,0);
			c.g = image(x,y,0,1);
			c.b = image(x,y,0,2);
			mask(x,y) = imgMaskTest(c) ? IS_MASKED : IS_NOT_MASKED; // TODO: use alpha channel if available
		}
	}

	int count = 0;
	cimg_forXY(mask,x,y)
	{
		if ( mask(x,y) == IS_MASKED )
		{
			++count; 
			// TODO: deal with any number of mask
                        if ( count == IS_MASKED ) { std::cerr << "ERROR: too many mask mask" << std::endl; exit(EXIT_FAILURE); }
			uchar regionID = count;
                        mask.draw_fill(x, y, &regionID);

                        /*
                        char * name = (char*)malloc(128);
                        sprintf(name, "debug_mask_%03d.png", count);
                        mask.save(name);
                        free(name);
                        */
		}
	}

	std::cout << "Found " << count << " mask region." << std::endl;

	CImg<uchar> * maskOrig;

	if ( !dilation && inpaintingFunc == 2 )
	{
		dilation = 2;
	}
	if ( dilation )
	{
		// TODO: this is a hack to overcome gradient being taken without respect to boundary
		maskOrig = new CImg<uchar>(mask, false); // force copy
		mask.dilate(dilation * 2 + 1);
	}

        // DEBUG
	/*
        image.save("debug_image.png");
	mask.save("debug_mask.png");
        */

	// Find region 4-boundaries
	std::vector<Coords> boundaries;
	int W = mask.width() - 1;
	int H = mask.height() - 1;

	int x_min = limitInpainting ? W : 0;
	int x_max = limitInpainting ? 0 : W;
	int y_min = limitInpainting ? H : 0;
	int y_max = limitInpainting ? 0 : H;

	for ( int i = 0; i < count; ++i )
	{
		ChainManager cm;
		uchar regionID = i + 1;
		cimg_forXY(mask, x, y)
		{
			// TODO: efficiencies
			if ( mask(x, y) != regionID )
			{
				// TODO: functionize
				if (  (
					( y > 0 && mask(x, y - 1) == regionID ) ||
					( y < H && mask(x, y + 1) == regionID ) ||
					( x > 0 && mask(x - 1, y) == regionID ) ||
					( x < W && mask(x + 1, y) == regionID )
				      )
				)
				{
					// hack to cope with the not-smart-enough boundary chain manager
					int eightNeighbours = 4;
					if (
						( y > 0 && mask(x, y - 1) == regionID ) &&
						( y < H && mask(x, y + 1) == regionID ) &&
						( x > 0 && mask(x - 1, y) == regionID ) &&
						( x < W && mask(x + 1, y) == regionID )
					)
					{
						if ( y > 0 && x > 0 && mask(x - 1, y - 1) == regionID ) eightNeighbours--;
						if ( y > 0 && x < W && mask(x + 1, y - 1) == regionID ) eightNeighbours--;
						if ( y < H && x < W && mask(x + 1, y + 1) == regionID ) eightNeighbours--;
						if ( y < H && x > 0 && mask(x - 1, y + 1) == regionID ) eightNeighbours--;
					}

					if ( eightNeighbours >= 2)
					{
						Color color(
							image(x, y, 0, 0),
							image(x, y, 0, 1),
							image(x, y, 0, 2)
						);
						cm.add(std::make_shared<Coord>(x, y, color));

						if ( limitInpainting )
						{
							if ( x < x_min ) x_min = x;
							if ( x > x_max ) x_max = x;
							if ( y < y_min ) y_min = y;
							if ( y > y_max ) y_max = y;
						}
					}
				}
			}
		}
		if ( ! cm.isGood(image.width(), image.height()) )
			std::cout << "ERROR: boundary ordering failure (fragmented or neither a loop nor spanning the image) - simplify mask and try again?" << std::endl;
		boundaries.push_back(cm.orderedChains()); // TODO: cope with multiple boundaries for a region
	}

	std::cout << "Calculated boundaries." << std::endl;

	if ( display )
	{
		cimg_library::CImgDisplay mask_disp(mask,"Mask");
		cimg_library::CImgDisplay debug_disp(image, "Debug");
	}

	// Inpainting
	int regionID = 0;
	std::vector<Inpainter*> inpainters;
	for ( auto bi = boundaries.begin(); bi != boundaries.end(); bi++ )
	{
		regionID++;
		std::cout << "DEBUG: length of boundary = " << (*bi).size() << std::endl;

		switch ( inpaintingFunc )
		{
		case 0:
			inpainters.push_back(new BleedInpainter(&*bi));
			break;
		case 1:
			inpainters.push_back(new WeightedInpainter(&*bi));
			break;
		case 2:
			inpainters.push_back(new GradientWeightedInpainter(&*bi, &image, maskOrig, &mask, regionID, smoothness, jitter, dilation));
			break;
		default:
			inpainters.push_back(new BleedInpainter(&*bi));
		}
	}

	std::cout << "Created inpainters." << std::endl;

	float i = 0;
	float i_total = (x_max - x_min + 1) * (y_max - y_min + 1);
	printf("Inpainting:      0%% ");

	// TODO: allow 'blend mode' where mask value is a blend coefficient

	for ( int y = y_min; y <= y_max; ++y)
	{
		for ( int x = x_min; x <= x_max; ++x)
		{
			if ( mask(x,y) > 0 )
			{
				// process associated boundary
				int r = mask(x,y) - 1;
				// TODO: check r is within bounds
				Color c = inpainters[r]->pixelColor(std::make_shared<Coord>(x, y));
				image(x,y,0,0) = c.r;
				image(x,y,0,1) = c.g;
				image(x,y,0,2) = c.b;
			}
			++i;
		}
		for ( int j = 0; j < 8; ++j ) printf("%c", 8);
		printf("% 6.0f%% ", 100 * i / (i_total - 1));
	}

	std::cout << "complete." << std::endl;

	if ( limitInpainting == 2 ) image.crop(x_min, y_min, x_max, y_max);

	if ( outputFilename )
	{
		image.save(outputFilename);
		mask.save("mask_with_boundary.png");
	}

	if ( display )
	{
		cimg_library::CImgDisplay output_disp(image,"Inpainted Original");
		while ( !output_disp.is_closed() )
		{
			output_disp.wait(1000);
		}
	}

	return 0;
}

