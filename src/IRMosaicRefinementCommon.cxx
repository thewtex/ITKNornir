// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: t -*-
// NOTE: the first line of this file sets up source code indentation rules
// for Emacs; it is also a hint to anyone modifying this file.

/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


// File         : mosaic_refinement_common.cxx
// Author       : Pavel A. Koshevoy
// Created      : Mon Nov  3 20:26:25 MST 2008
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Helper functions for automatic mosaic refinement.

// local includes:
#include "IRMosaicRefinementCommon.h"


//----------------------------------------------------------------
// regularize_displacements
//
void
regularize_displacements( // computed displacement vectors of the moving image
                          // grid transform control points, in mosaic space:
  std::vector<vec2d_t> & xy_shift,
  std::vector<double> &  mass,

  image_t::Pointer & dx,
  image_t::Pointer & dy,
  image_t::Pointer & db,

  // median filter radius:
  const unsigned int & median_radius)
{
  // shortcuts:
  image_t::RegionType::SizeType sz = dx->GetLargestPossibleRegion().GetSize();
  unsigned int                  mesh_cols = sz[0];
  unsigned int                  mesh_rows = sz[1];

  // denoise
  if (median_radius > 0)
  {
    dx = median<image_t>(dx, median_radius);
    dy = median<image_t>(dy, median_radius);
    // db = median<image_t>(db, median_radius);
  }

  // extend (fill in gaps):
  typedef itk::ImageRegionConstIteratorWithIndex<image_t> iter_t;
  iter_t                                                  iter(dx, dx->GetLargestPossibleRegion());
  image_t::Pointer                                        dx_blurred = cast<image_t, image_t>(dx);
  image_t::Pointer                                        dy_blurred = cast<image_t, image_t>(dy);
  image_t::Pointer                                        db_blurred = cast<image_t, image_t>(db);

  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    image_t::IndexType index = iter.GetIndex();
    if (!db->GetPixel(index))
    {
      static const double max_w = 3.0;
      double              w = 0.0;
      double              px = 0.0;
      double              py = 0.0;

      // keep expanding the neighborhood until at least one
      // successful sample is found:

      int max_x = std::max(int(index[0]), int(mesh_cols - 1 - index[0]));
      int max_y = std::max(int(index[1]), int(mesh_rows - 1 - index[1]));
      int max_r = std::min(1, std::max(max_x, max_y));
      for (int r = 1; r <= max_r && w < max_w; r++)
      {
        image_t::IndexType ix;
        int                x0 = index[0] - r;
        int                x1 = index[0] + r;
        int                y0 = index[1] - r;
        int                y1 = index[1] + r;

        int d = 2 * r + 1;
        for (int o = 0; o < d; o++)
        {
          ix[0] = x0;
          ix[1] = y0 + o + 1;
          if (ix[0] >= 0 && ix[0] < int(mesh_cols) && ix[1] >= 0 && ix[1] < int(mesh_rows) && db->GetPixel(ix))
          {
            px += dx->GetPixel(ix);
            py += dy->GetPixel(ix);
            w += 1.0;
          }

          ix[0] = x1;
          ix[1] = y0 + o;
          if (ix[0] >= 0 && ix[0] < int(mesh_cols) && ix[1] >= 0 && ix[1] < int(mesh_rows) && db->GetPixel(ix))
          {
            px += dx->GetPixel(ix);
            py += dy->GetPixel(ix);
            w += 1.0;
          }

          ix[0] = x0 + o;
          ix[1] = y0;
          if (ix[0] >= 0 && ix[0] < int(mesh_cols) && ix[1] >= 0 && ix[1] < int(mesh_rows) && db->GetPixel(ix))
          {
            px += dx->GetPixel(ix);
            py += dy->GetPixel(ix);
            w += 1.0;
          }

          ix[0] = x0 + o + 1;
          ix[1] = y1;
          if (ix[0] >= 0 && ix[0] < int(mesh_cols) && ix[1] >= 0 && ix[1] < int(mesh_rows) && db->GetPixel(ix))
          {
            px += dx->GetPixel(ix);
            py += dy->GetPixel(ix);
            w += 1.0;
          }
        }
      }

      if (w != 0.0)
      {
        dx_blurred->SetPixel(index, px / w);
        dy_blurred->SetPixel(index, py / w);
        db_blurred->SetPixel(index, 1);
      }
    }
  }

  // blur:
  dx_blurred = smooth<image_t>(dx_blurred, 1.0);
  dy_blurred = smooth<image_t>(dy_blurred, 1.0);
  // db_blurred = smooth<image_t>(db_blurred, 1.0);

  dx = dx_blurred;
  dy = dy_blurred;
  db = db_blurred;

  // update the mesh displacement field:
  iter = iter_t(dx, dx->GetLargestPossibleRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    image_t::IndexType index = iter.GetIndex();
    unsigned int       i = index[0] + index[1] * mesh_cols;

    xy_shift[i][0] = dx->GetPixel(index);
    xy_shift[i][1] = dy->GetPixel(index);
    mass[i] += db->GetPixel(index);
  }
}
