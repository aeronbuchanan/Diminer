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

#include "diminer.h"
#include "patch.h"
#include "inpainters.h"

using namespace Diminer;

int main(int argc, char** argv)
{
	std::cout << "Diminer: Digital Image Inpainting Resources by Aeron Buchanan" << std::endl;

	cimg_usage("Usage: Diminer [options] -i input\n");

	char const * imageFilename = cimg_option("-i", "image.jpg", "image file to be inpainted");
	char const * maskFilename = cimg_option("-m", "", "mask image denoting region to be inpainted (values > 127)");
	char const * outputFilename = cimg_option("-o", "inpainted.jpg", "output filename");
	uint inpaintingFunc = cimg_option("-f", 2, "0 = bleed; 1 = weighted; 2 = gradient-weighted");
	float jitter = cimg_option("-j", 0.35, "jitter between 0.f and 1.f");
	float smoothness = cimg_option("-s", 2.f, "smoothness");
	int dilation = cimg_option("-r", 0, "mask dilation radius");
	bool display = cimg_option("-D", 0, "display");

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

	CImg<uchar> * maskOrig;

	if ( !dilation && inpaintingFunc == 2 )
	{
		dilation = 2;
	}
	if ( dilation )
	{
		// TODO: this is a hack to overcome gradient being taken without respect to boundary
		maskOrig = new CImg<uchar>(mask);
		mask.dilate(dilation * 2 + 1);
	}

	//*** detect connected regions ***//
	CImg<uchar> regions(mask); // TODO: float matrix

        regions.save("debug_regions_000.png");

	int count = 0;
	cimg_forXY(mask,x,y)
	{
		if ( regions(x,y) == IS_MASKED )
		{
			++count; // TODO: deal with any number of regions
                        if ( count == IS_MASKED ) { std::cerr << "ERROR: too many mask regions" << std::endl; return 1; }
			uchar col = count; // 255 / count;
                        regions.draw_fill(x, y, &col);

                        char * name = (char*)malloc(128);
                        sprintf(name, "debug_regions_%03d.png", count);
                        regions.save(name);
                        free(name);
		}
	}

	std::cout << "Found " << count << " regions." << std::endl;

        // DEBUG
	image.save("debug_image.png");
	mask.save("debug_mask.png");
        regions.save("debug_regions.png");

	// Find region 4-boundaries
	std::vector<BoundaryColors> boundaries;
	int W = regions.width() - 1;
	int H = regions.height() - 1;

	int x_min = W;
	int x_max = 0;
	int y_min = H;
	int y_max = 0;

	for ( int i = 0; i < count; ++i )
	{
		BoundaryColors cs;
		uchar col = i + 1; // 255 / count;
		cimg_forXY(regions, x, y)
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

					if ( x < x_min ) x_min = x;
					if ( x > x_max ) x_max = x;
					if ( y < y_min ) y_min = y;
					if ( y > y_max ) y_max = y;
				}
			}
		}
		boundaries.push_back(cs);
	}

	std::cout << "Calculated boundaries." << std::endl;

	if ( display )
	{
		cimg_library::CImgDisplay mask_disp(mask,"Mask");
		cimg_library::CImgDisplay debug_disp(image, "Debug");
	}

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
			inpainters.push_back(new GradientWeightedInpainter(&boundaries[i], &image, maskOrig, smoothness, jitter, dilation));
			break;
		default:
			inpainters.push_back(new BleedInpainter(&boundaries[i]));
		}
	}

	std::cout << "Created inpainters." << std::endl;

	float i = 0;
	float i_total = (x_max - x_min + 1) * (y_max - y_min + 1);
	std::cout << "Inpainting: ";

	// TODO: allow 'blend mode' where mask value is a blend coefficient

	for ( int y = y_min; y <= y_max; ++y)
	{
		for ( int x = x_min; x <= x_max; ++x)
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
	}

	std::cout << "complete." << std::endl;

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

