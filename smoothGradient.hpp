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

#include <math.h>
#include <limits>
#include <tuple>

#include "diminer.h"
#include "CImg.h" 

namespace Diminer
{

void smoothGradient(SourceImage const * const _img, MaskImage const * const _mask, GradImage * const _output)
{
	std::cout << "Creating gradient image" << std::endl;

	_output->assign(_img->width(), _img->height(), 1, 3, 0.f);

	// greyscale
	CImg<float> g(_img->width(), _img->height(), 1, 1, 0.f);
	cimg_forXYC((*_img), x, y, k)
		g(x, y) += _img->operator()(x, y, k);

/*
	@misc{NRIGO,
	Author = {Pavel Holoborodko},
	Title = {Noise Robust Gradient Operators.},
	howpublished = {\url{http://www.holoborodko.com/pavel/image-processing/edge-detection/}}
	year = {2009}
	}

	dI/dx = 1/32 [ -1 -2 0 2 1; -2 -4 0 4 2; -1 -2 0 2 1 ]

	(like an extended Sobel kernel)
*/

	std::cout << "Looping..." << std::endl;

	// TODO: only calculate near boundary
	// TODO: optimize recalculations
	for ( int y = 0; y < g.height(); ++y )
	{
		int y_pp = std::max( 0, y - 2);
		int y_p = std::max( 0, y - 1);
		int y_n = std::min( g.height() - 1, y + 1);
		int y_nn = std::min( g.height() - 1, y + 2);

		for ( int x = 0; x < g.width(); ++x )
		{
			int x_pp = std::max( 0, x - 2);
			int x_p = std::max( 0, x - 1);
			int x_n = std::min( g.width() - 1, x + 1);
			int x_nn = std::min( g.width() - 1, x + 2);

			double gxy = g(x, y);

			double gxpyp = (*_mask)(x_p, y_p) ? gxy : g(x_p, y_p);
			double gxppyp = (*_mask)(x_pp, y_p) ? gxpyp : g(x_pp, y_p);
			double gxnyp = (*_mask)(x_n, y_p) ? gxy : g(x_n, y_p);
			double gxnnyp = (*_mask)(x_nn, y_p) ? gxnyp : g(x_nn, y_p);

			double gxpy = (*_mask)(x_p, y) ? gxy : g(x_p, y);
			double gxppy = (*_mask)(x_pp, y) ? gxpy : g(x_pp, y);
			double gxny = (*_mask)(x_n, y) ? gxy : g(x_n, y);
			double gxnny = (*_mask)(x_nn, y) ? gxny : g(x_nn, y);

			double gxpyn = (*_mask)(x_p, y_n) ? gxy : g(x_p, y_n);
			double gxppyn = (*_mask)(x_pp, y_n) ? gxpyn : g(x_pp, y_n);
			double gxnyn = (*_mask)(x_n, y_n) ? gxy : g(x_n, y_n);
			double gxnnyn = (*_mask)(x_nn, y_n) ? gxnyn : g(x_nn, y_n);

			double gx = ( gxnnyn + gxnnyp + 2 * (gxnyn + gxnny + gxnyp) + 4 * gxny ) - ( gxppyp + gxppyn + 2 * (gxppy + gxpyp + gxpyn) + 4 * gxpy );
			gx /= 32;

			//double gxpyp;
			double gxpypp = (*_mask)(x_p, y_pp) ? gxpyp : g(x_p, y_pp, 0);
			//double gxpyn;
			double gxpynn = (*_mask)(x_p, y_nn) ? gxpyn : g(x_p, y_nn, 0);

			double gxyp = (*_mask)(x, y_p) ? gxy : g(x, y_p, 0);
			double gxypp = (*_mask)(x, y_pp) ? gxyp : g(x, y_pp, 0);
			double gxyn = (*_mask)(x, y_n) ? gxy : g(x, y_n, 0);
			double gxynn = (*_mask)(x, y_nn) ? gxyn : g(x, y_nn, 0);

			//double gxnyp;
			double gxnypp = (*_mask)(x_n, y_pp) ? gxnyp : g(x_n, y_pp, 0);
			//double gxnyn;
			double gxnynn = (*_mask)(x_n, y_nn) ? gxnyn : g(x_n, y_nn, 0);

			double gy = ( gxnynn + gxpynn + 2 * (gxnyn + gxynn + gxpyn) + 4 * gxyn ) - ( gxpypp + gxnypp + 2 * (gxypp + gxpyp + gxnyp) + 4 * gxyp );
			gy /= 32;

			double mag = sqrt( gx * gx + gy * gy );
			(*_output)(x, y, 0, 0) = mag;
			(*_output)(x, y, 0, 1) = gx / mag;
			(*_output)(x, y, 0, 2) = gy / mag;
		}
	}

	std::cout << "Done." << std::endl;

	_output->save("debug_gradient.png");
}


} // end namespace Diminer
