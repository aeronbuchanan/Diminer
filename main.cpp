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

	if ( !dilation && inpaintingFunc == 2 )
	{
		dilation = 2;
	}
	if ( dilation )
	{
		// TODO: this is a hack to overcome gradient being taken without respect to boundary
		maskOrig = new CImg<uchar>(regions);
		cimg_forXY(mask,x,y)
			(*maskOrig)(x,y) = maskTest(mask(x,y));

		mask.dilate(dilation * 2 + 1);
	}

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
		CImgDisplay mask_disp(mask,"Mask");
		CImgDisplay debug_disp(image, "Debug");
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

	for ( uint y = y_min; y <= y_max; ++y)
	{
		for ( uint x = x_min; x <= x_max; ++x)
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
		CImgDisplay output_disp(image,"Inpainted Original");
		while ( !output_disp.is_closed() )
		{
			output_disp.wait(1000);
		}
	}

	return 0;
}

