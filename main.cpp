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

#include "smoothGradient.hpp"

using namespace Diminer;

int main(int argc, char** argv)
{
	std::cout << "Diminer: Digital Image Inpainting Resources by Aeron Buchanan" << std::endl;

	cimg_usage("Usage: Diminer [options] -i input\n");

	// TODO: use an inpainter factory
	char * ipFuncs = (char *)malloc(998);
	sprintf(ipFuncs, "%d = bleed; %d = blend; %d = gradient-weighted", Inpainters::BLEED, Inpainters::BLEND, Inpainters::GRADS);

	char const * imageFilename = cimg_option("-i", "image.jpg", "image file to be inpainted");
	char const * maskFilename = cimg_option("-m", "", "mask image denoting region to be inpainted (values > 127)");
	char const * outputFilename = cimg_option("-o", "inpainted.jpg", "output filename");
	int inpaintingFunc = cimg_option("-f", 2, ipFuncs); free(ipFuncs);
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

	SourceImage image(imageFilename);
	MaskImage mask(image.width(), image.height(), 1, 1, 0);
	if (maskFilename)
	{
		MaskImage maskFile(maskFilename);

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

	int regionCount = 0;
	cimg_forXY(mask,x,y)
	{
		if ( mask(x,y) == IS_MASKED )
		{
			++regionCount; 
			// TODO: deal with any number of mask
                        if ( regionCount == IS_MASKED ) { std::cerr << "ERROR: too many mask mask" << std::endl; exit(EXIT_FAILURE); }
			uchar regionID = regionCount;
                        mask.draw_fill(x, y, &regionID);

                        /*
                        char * name = (char*)malloc(128);
                        sprintf(name, "debug_mask_%03d.png", count);
                        mask.save(name);
                        free(name);
                        */
		}
	}

	std::cout << "Found " << regionCount << " mask region." << std::endl;

	MaskImage * maskOrig;

	if ( !dilation && inpaintingFunc == 2 )
	{
		dilation = 2;
	}
	if ( dilation )
	{
		// TODO: this is a hack to overcome gradient being taken without respect to boundary
		maskOrig = new MaskImage(mask, false); // force copy
		mask.dilate(dilation * 2 + 1);
	}

        // DEBUG
	/*
        image.save("debug_image.png");
	mask.save("debug_mask.png");
        */

	BoundaryManager master(&image, &mask, regionCount);
	Boundaries * boundaries = master.boundaries();

	int x_min = limitInpainting ? master.x_min() : 0;
	int x_max = limitInpainting ? master.x_max() : image.width() - 1;
	int y_min = limitInpainting ? master.y_min() : 0;
	int y_max = limitInpainting ? master.y_max() : image.height() - 1;

	std::cout << "Calculated boundaries." << std::endl;

	if ( display )
	{
		cimg_library::CImgDisplay mask_disp(mask,"Mask");
		cimg_library::CImgDisplay debug_disp(image, "Debug");
	}

	// Inpainting
	int regionID = 0;
	std::vector<Inpainter*> inpainters;

	GradImage * grads;
	if ( inpaintingFunc == Inpainters::GRADS )
	{
		grads = new GradImage();
		smoothGradient(&image, maskOrig, grads);
	}

	for ( auto bi = boundaries->begin(); bi != boundaries->end(); bi++ )
	{
		regionID++;
		std::cout << "DEBUG: length of boundary = " << (*bi).size() << std::endl;

		switch ( inpaintingFunc )
		{
		case Inpainters::BLEED:
			inpainters.push_back(new BleedInpainter(&*bi));
			break;
		case Inpainters::BLEND:
			inpainters.push_back(new WeightedInpainter(&*bi));
			break;
		case Inpainters::GRADS:
			inpainters.push_back(new GradientWeightedInpainter(&*bi, &image, maskOrig, &mask, grads, regionID, smoothness, jitter, dilation));
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

