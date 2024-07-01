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


// File         : itkGridTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2006/11/20 22:48
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A discontinuous transform -- a uniform grid of vertices is
//                mapped to an image. At each vertex, in addition to image
//                space coordinates, a second set of coordinates is stored.
//                This is similar to texture mapped OpenGL triangle meshes,
//                where the texture coordinates correspond to the image space
//                vertex coordinates.

#ifndef __itkGridTransform_h
#define __itkGridTransform_h

// system includes:
#include <math.h>
#include <iostream>
#include <assert.h>

// ITK includes:
#include <itkTransform.h>
#include <itkIdentityTransform.h>
#include <itkMacro.h>
#include <itkImage.h>

// local includes:
#include "IRGridTransform.h"


//----------------------------------------------------------------
// itk::GridTransform
//
namespace itk
{
class GridTransform : public Transform<double, 2, 2>
{
public:
  // standard typedefs:
  typedef GridTransform            Self;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef Transform<double, 2, 2> Superclass;

  // Base inverse transform type:
  typedef Superclass                         InverseTransformType;
  typedef SmartPointer<InverseTransformType> InverseTransformPointer;

  // RTTI:
  itkTypeMacro(GridTransform, Transform);

  // macro for instantiation through the object factory:
  itkNewMacro(Self);

  /** Standard scalar type for this class. */
  typedef double ScalarType;

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 2);

  // shortcuts:
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::JacobianType   JacobianType;

  typedef Superclass::InputPointType  InputPointType;
  typedef Superclass::OutputPointType OutputPointType;

  // virtual:
  OutputPointType
  TransformPoint(const InputPointType & x) const
  {
    OutputPointType y;
    if (is_inverse())
    {
      pnt2d_t uv;
      uv[0] = (x[0] - transform_.tile_min_[0]) / transform_.tile_ext_[0];
      uv[1] = (x[1] - transform_.tile_min_[1]) / transform_.tile_ext_[1];
      transform_.transform_inv(uv, y);
    }
    else
    {
      transform_.transform(x, y);
      y[0] *= transform_.tile_ext_[0];
      y[1] *= transform_.tile_ext_[1];
      y[0] += transform_.tile_min_[0];
      y[1] += transform_.tile_min_[1];
    }

    // ITK does not handle NaN numbers:
    if (y[0] != y[0])
    {
      y[0] = std::numeric_limits<double>::max();
      y[1] = y[0];
    }

    return y;
  }

  // Inverse transformations:
  // If y = Transform(x), then x = BackTransform(y);
  // if no mapping from y to x exists, then an exception is thrown.
  InputPointType
  BackTransformPoint(const OutputPointType & y) const
  {
    InputPointType x;
    if (is_inverse())
    {
      transform_.transform(y, x);
      x[0] *= transform_.tile_ext_[0];
      x[1] *= transform_.tile_ext_[1];
      x[0] += transform_.tile_min_[0];
      x[1] += transform_.tile_min_[1];
    }
    else
    {
      pnt2d_t uv;
      uv[0] = (y[0] - transform_.tile_min_[0]) / transform_.tile_ext_[0];
      uv[1] = (y[1] - transform_.tile_min_[1]) / transform_.tile_ext_[1];
      transform_.transform_inv(uv, x);
    }

    // ITK does not handle NaN numbers:
    if (x[0] != x[0])
    {
      x[0] = std::numeric_limits<double>::max();
      x[1] = x[0];
    }

    return x;
  }

  using InputDiffusionTensor3DType = typename Superclass::InputDiffusionTensor3DType;
  using OutputDiffusionTensor3DType = typename Superclass::OutputDiffusionTensor3DType;
  OutputDiffusionTensor3DType
  TransformDiffusionTensor3D(const InputDiffusionTensor3DType & inputTensor, const InputPointType & point) const override
  {
    itkExceptionMacro("TransformDiffusionTensor3D( const InputDiffusionTensor3DType & ) is "
                      "unimplemented for "
                      << this->GetNameOfClass());
  }

  using InputVectorPixelType = typename Superclass::InputVectorPixelType;
  using OutputVectorPixelType = typename Superclass::OutputVectorPixelType;
  OutputVectorPixelType
  TransformDiffusionTensor3D(const InputVectorPixelType & inputTensor, const InputPointType & point) const override
  {
    itkExceptionMacro("TransformDiffusionTensor3D( const InputVectorPixelType & ) is "
                      "unimplemented for "
                      << this->GetNameOfClass());
  }


  // virtual:
  void
  SetFixedParameters(const ParametersType & params)
  {
    this->m_FixedParameters = params;
  }

  // virtual:
  const ParametersType &
  GetFixedParameters() const
  {
    ParametersType params = this->m_FixedParameters;

    params[1] = transform_.rows_;
    params[2] = transform_.cols_;
    params[3] = transform_.tile_min_[0];
    params[4] = transform_.tile_min_[1];
    params[5] = transform_.tile_ext_[0];
    params[6] = transform_.tile_ext_[1];

    // update the parameters vector:
    GridTransform * fake = const_cast<GridTransform *>(this);
    fake->m_FixedParameters = params;
    return this->m_FixedParameters;
  }

  // virtual:
  void
  SetParameters(const ParametersType & params)
  {
    this->m_Parameters = params;

    size_t rows = int(this->m_FixedParameters[1]);
    size_t cols = int(this->m_FixedParameters[2]);

    pnt2d_t tile_min;
    tile_min[0] = this->m_FixedParameters[3];
    tile_min[1] = this->m_FixedParameters[4];

    pnt2d_t tile_max;
    tile_max[0] = tile_min[0] + this->m_FixedParameters[5];
    tile_max[1] = tile_min[1] + this->m_FixedParameters[6];

    std::vector<pnt2d_t> xy((rows + 1) * (cols + 1));

    unsigned int num_points = xy.size();
    assert(2 * num_points == params.size());
    for (unsigned int i = 0; i < num_points; i++)
    {
      xy[i][0] = params[i * 2];
      xy[i][1] = params[i * 2 + 1];
    }

    transform_.setup(rows, cols, tile_min, tile_max, xy);
  }

  // virtual:
  const ParametersType &
  GetParameters() const
  {
    ParametersType params(GetNumberOfParameters());
    unsigned int   num_pts = params.size() / 2;
    unsigned int   num_cols = (transform_.cols_ + 1);
    for (unsigned int i = 0; i < num_pts; i++)
    {
      unsigned int row = i / num_cols;
      unsigned int col = i % num_cols;

      const pnt2d_t & xy = transform_.vertex(row, col).xy_;
      unsigned int    idx = 2 * i;
      params[idx] = xy[0];
      params[idx + 1] = xy[1];
    }

    // update the parameters vector:
    GridTransform * fake = const_cast<GridTransform *>(this);
    fake->m_Parameters = params;
    return this->m_Parameters;
  }

  // virtual:
  Superclass::NumberOfParametersType
  GetNumberOfParameters() const override
  {
    return 2 * transform_.grid_.mesh_.size();
  }

  // virtual: return an inverse of this transform.
  InverseTransformPointer
  GetInverse() const
  {
    GridTransform::Pointer inv = GridTransform::New();
    inv->setup(transform_, !is_inverse());
    return inv.GetPointer();
  }

  // setup the transform:
  void
  setup(const the_grid_transform_t & transform, const bool & is_inverse = false)
  {
    transform_ = transform;
    GetParameters();
    GetFixedParameters();
    this->m_FixedParameters[0] = is_inverse ? 1.0 : 0.0;
  }

  // inverse transform flag check:
  inline bool
  is_inverse() const
  {
    return this->m_FixedParameters[0] != 0.0;
  }

  void
  ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & jacobian) const override
  {
    // FIXME: 20061227 -- this function was written and not tested:

    // these scales are necessary to account for the fact that
    // the_grid_transform_t expects uv in the [0,1]x[0,1] range,
    // where as we remap it into the image tile physical coordinates
    // according to tile_min_ and tile_ext_:
    double Su = transform_.tile_ext_[0];
    double Sv = transform_.tile_ext_[1];

    unsigned int idx[3];
    double       jac[12];
    jacobian.SetSize(2, GetNumberOfParameters());
    jacobian.Fill(0.0);
    if (transform_.jacobian(point, idx, jac))
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        unsigned int addr = idx[i] * 2;
        jacobian(0, addr) = Su * jac[i * 2];
        jacobian(0, addr + 1) = Su * jac[i * 2 + 1];
        jacobian(1, addr) = Sv * jac[i * 2 + 6];
        jacobian(1, addr + 1) = Sv * jac[i * 2 + 7];
      }
    }
  }

protected:
  GridTransform()
  {
    this->m_FixedParameters.SetSize(7);

    // initialize the inverse flag:
    this->m_FixedParameters[0] = 0.0;

    // grid size:
    this->m_FixedParameters[1] = 0.0;
    this->m_FixedParameters[2] = 0.0;

    // tile bbox:
    this->m_FixedParameters[3] = std::numeric_limits<double>::max();
    this->m_FixedParameters[4] = this->m_FixedParameters[3];
    this->m_FixedParameters[5] = -(this->m_FixedParameters[3]);
    this->m_FixedParameters[6] = -(this->m_FixedParameters[3]);
  }

private:
  // disable default copy constructor and assignment operator:
  GridTransform(const Self & other);
  const Self &
  operator=(const Self & t);

public:
  // the actual transform:
  the_grid_transform_t transform_;

}; // class GridTransform

} // namespace itk


#endif //  __itkGridTransform_h
