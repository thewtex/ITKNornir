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


// File         : the_aa_bbox.hxx
// Author       : Pavel A. Koshevoy
// Created      : Mon Jun  7 22:14:00 MDT 2004
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Axis aligned bounding box.

// system includes:
#include <assert.h>
#include <algorithm>

// local includes:
#include "IRAABBox.h"


///----------------------------------------------------------------
// the_aa_bbox_t::operator <<
// 
the_aa_bbox_t &
the_aa_bbox_t::operator << (const p3x1_t & pt)
{
  if (min_.x() > pt.x()) min_.x() = pt.x();
  if (max_.x() < pt.x()) max_.x() = pt.x();
  if (min_.y() > pt.y()) min_.y() = pt.y();
  if (max_.y() < pt.y()) max_.y() = pt.y();
  if (min_.z() > pt.z()) min_.z() = pt.z();
  if (max_.z() < pt.z()) max_.z() = pt.z();
  
  return *this;
}

//----------------------------------------------------------------
// the_aa_bbox_t::operator +=
// 
the_aa_bbox_t &
the_aa_bbox_t::operator += (const the_aa_bbox_t & bbox)
{
  if (bbox.is_empty()) return *this;
  return (*this) << bbox.min_ << bbox.max_;
}

//----------------------------------------------------------------
// the_aa_bbox_t::operator *=
// 
the_aa_bbox_t &
the_aa_bbox_t::operator *= (const float & s)
{
  if (is_empty()) return *this;
  
  p3x1_t o = min_ + 0.5 * (max_ - min_);
  min_ = o + s * (min_ - o);
  max_ = o + s * (max_ - o);
  return *this;
}

//----------------------------------------------------------------
// the_aa_bbox_t::operator +=
// 
the_aa_bbox_t &
the_aa_bbox_t::operator += (const float & r)
{
  if (is_empty())
  {
    min_.assign(-r, -r, -r);
    max_.assign( r,  r,  r);
  }
  else
  {
    v3x1_t shift(r, r, r);
    min_ -= shift;
    max_ += shift;
  }
  
  assert((min_.x() < max_.x()) &&
	 (min_.y() < max_.y()) &&
	 (min_.z() < max_.z()));
  
  return *this;
}

//----------------------------------------------------------------
// the_aa_bbox_t::corners
// 
void
the_aa_bbox_t::corners(p3x1_t * corner) const
{
  corner[0] = p3x1_t(max_.x(), min_.y(), max_.z());
  corner[1] = p3x1_t(min_.x(), min_.y(), max_.z());
  corner[2] = p3x1_t(min_.x(), min_.y(), min_.z());
  corner[3] = p3x1_t(max_.x(), min_.y(), min_.z());
  corner[4] = p3x1_t(max_.x(), max_.y(), max_.z());
  corner[5] = p3x1_t(min_.x(), max_.y(), max_.z());
  corner[6] = p3x1_t(min_.x(), max_.y(), min_.z());
  corner[7] = p3x1_t(max_.x(), max_.y(), min_.z());
}

//----------------------------------------------------------------
// the_aa_bbox_t::is_empty
// 
bool
the_aa_bbox_t::is_empty() const
{
  return ((min_.x() == FLT_MAX) &&
	  (min_.y() == FLT_MAX) &&
	  (min_.z() == FLT_MAX) &&
	  (max_.x() == -FLT_MAX) &&
	  (max_.y() == -FLT_MAX) &&
	  (max_.z() == -FLT_MAX));
}

//----------------------------------------------------------------
// the_aa_bbox_t::is_linear
// 
bool
the_aa_bbox_t::is_linear() const
{
  if (is_empty()) return false;
  
  float dx = max_.x() - min_.x();
  float dy = max_.y() - min_.y();
  float dz = max_.z() - min_.z();
  return (((dx != 0.0) && (dy == 0.0) && (dz == 0.0)) ||
	  ((dx == 0.0) && (dy != 0.0) && (dz == 0.0)) ||
	  ((dx == 0.0) && (dy == 0.0) && (dz != 0.0)));
}

//----------------------------------------------------------------
// the_aa_bbox_t::is_planar
// 
bool
the_aa_bbox_t::is_planar() const
{
  if (is_empty()) return false;
  
  float dx = max_.x() - min_.x();
  float dy = max_.y() - min_.y();
  float dz = max_.z() - min_.z();
  return (((dx == 0.0) && (dy != 0.0) && (dz != 0.0)) ||
	  ((dx != 0.0) && (dy == 0.0) && (dz != 0.0)) ||
	  ((dx != 0.0) && (dy != 0.0) && (dz == 0.0)));
}

//----------------------------------------------------------------
// the_aa_bbox_t::is_spacial
// 
bool
the_aa_bbox_t::is_spacial() const
{
  if (is_empty()) return false;
  
  float dx = max_.x() - min_.x();
  float dy = max_.y() - min_.y();
  float dz = max_.z() - min_.z();
  return ((dx != 0.0) && (dy != 0.0) && (dz != 0.0));
}

//----------------------------------------------------------------
// the_aa_bbox_t::radius
// 
float
the_aa_bbox_t::radius(const p3x1_t & center) const
{
  if (is_empty()) return 0.0;
  
  p3x1_t corner[] =
  {
    p3x1_t(min_.x(), min_.y(), min_.z()),
    p3x1_t(min_.x(), min_.y(), max_.z()),
    p3x1_t(min_.x(), max_.y(), min_.z()),
    p3x1_t(min_.x(), max_.y(), max_.z()),
    p3x1_t(max_.x(), min_.y(), min_.z()),
    p3x1_t(max_.x(), min_.y(), max_.z()),
    p3x1_t(max_.x(), max_.y(), min_.z()),
    p3x1_t(max_.x(), max_.y(), max_.z())
  };
  
  float max_dist = -FLT_MAX;
  for (unsigned int i = 0; i < 8; i++)
  {
    float dist = ~(center - corner[i]);
    if (dist > max_dist) max_dist = dist;
  }
  
  return max_dist;
}

//----------------------------------------------------------------
// the_aa_bbox_t::radius
//
float
the_aa_bbox_t::radius(const p3x1_t & center,
		      const unsigned int & axis_w_id) const
{
  assert(axis_w_id < 3);
  unsigned int axis_u_id = (axis_w_id + 1) % 3;
  unsigned int axis_v_id = (axis_u_id + 1) % 3;
  
  p2x1_t uv_center(center[axis_u_id],
		   center[axis_v_id]);
  
  p2x1_t uv_corner[] = {
    p2x1_t(min_[axis_u_id], min_[axis_v_id]),
    p2x1_t(min_[axis_u_id], max_[axis_v_id]),
    p2x1_t(max_[axis_u_id], min_[axis_v_id]),
    p2x1_t(max_[axis_u_id], max_[axis_v_id])
  };
  
  float max_dist = -FLT_MAX;
  for (unsigned int i = 0; i < 4; i++)
  {
    v2x1_t uv_vec = uv_center - uv_corner[i];
    float dist = ~uv_vec;
    if (dist > max_dist) max_dist = dist;
  }
  
  return max_dist;
}

//----------------------------------------------------------------
// the_aa_bbox_t::contains
// 
void
the_aa_bbox_t::contains(const p3x1_t & pt,
			bool & contained_in_x,
			bool & contained_in_y,
			bool & contained_in_z) const
{
  // check whether the point is within the bounding box limits:
  contained_in_x = ((pt.x() <= max_.x()) &&
		    (pt.x() >= min_.x()));
  contained_in_y = ((pt.y() <= max_.y()) &&
		    (pt.y() >= min_.y()));
  contained_in_z = ((pt.z() <= max_.z()) &&
		    (pt.z() >= min_.z()));
}

//----------------------------------------------------------------
// the_aa_bbox_t::intersects
// 
bool
the_aa_bbox_t::intersects(const the_aa_bbox_t & b) const
{
  bool disjoint =
    (min_[0] > b.max_[0] || b.min_[0] > max_[0]) ||
    (min_[1] > b.max_[1] || b.min_[1] > max_[1]) ||
    (min_[2] > b.max_[2] || b.min_[2] > max_[2]);
  
  return !disjoint;
}

//----------------------------------------------------------------
// the_aa_bbox_t::clamp
// 
void
the_aa_bbox_t::clamp(const the_aa_bbox_t & confines)
{
  if (confines.is_empty()) return;
  
  min_[0] = std::min(confines.max_[0], std::max(confines.min_[0], min_[0]));
  min_[1] = std::min(confines.max_[1], std::max(confines.min_[1], min_[1]));
  min_[2] = std::min(confines.max_[2], std::max(confines.min_[2], min_[2]));
  
  max_[0] = std::min(confines.max_[0], std::max(confines.min_[0], max_[0]));
  max_[1] = std::min(confines.max_[1], std::max(confines.min_[1], max_[1]));
  max_[2] = std::min(confines.max_[2], std::max(confines.min_[2], max_[2]));
}

//----------------------------------------------------------------
// the_aa_bbox_t::intersects_ray
// 
bool
the_aa_bbox_t::intersects_ray(const p3x1_t & o,
			      const v3x1_t & d,
			      float & t_min,
			      float & t_max) const
{
  float t[3][2];
  
  for (unsigned int i = 0; i < 3; i++)
  {
    float inv_d = 1.0f / d[i];
    unsigned int e = inv_d < 0.0f;
    
    t[i][e]           = (min_[i] - o[i]) * inv_d;
    t[i][(e + 1) % 2] = (max_[i] - o[i]) * inv_d;
  }
  
  t_min = std::max(t[0][0], std::max(t[1][0], t[2][0]));
  t_max = std::min(t[0][1], std::min(t[1][1], t[2][1]));
  
  return t_min < t_max;
}

//----------------------------------------------------------------
// the_aa_bbox_t::largest_dimension
// 
unsigned int
the_aa_bbox_t::largest_dimension() const
{
  float d[] =
  {
    max_[0] - min_[0],
    max_[1] - min_[1],
    max_[2] - min_[2]
  };
  
  unsigned int axis = 0;
  if (d[1] > d[0]) axis = 1;
  if (d[2] > d[axis]) axis = 2;
  
  return axis;
}

//----------------------------------------------------------------
// the_aa_bbox_t::smallest_dimension
// 
unsigned int
the_aa_bbox_t::smallest_dimension() const
{
  float d[] =
  {
    max_[0] - min_[0],
    max_[1] - min_[1],
    max_[2] - min_[2]
  };
  
  unsigned int axis = 0;
  if (d[1] < d[0]) axis = 1;
  if (d[2] < d[axis]) axis = 2;
  
  return axis;
}

//----------------------------------------------------------------
// the_aa_bbox_t::dump
// 
void
the_aa_bbox_t::dump(ostream & strm) const
{
  strm << "the_aa_bbox_t(" << (void *)this << ")" << endl
       << "{" << endl
       << "  min_ = " << min_ << endl
       << "  max_ = " << max_ << endl
       << "}" << endl;
}
