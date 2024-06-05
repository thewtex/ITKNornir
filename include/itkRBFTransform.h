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


// File         : itkRBFTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2007/01/23 10:15
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A Thin Plate Spline transform.

#ifndef __itkRBFTransform_h
#define __itkRBFTransform_h

// system includes:
#include <iostream>
#include <assert.h>

// ITK includes:
#include <itkTransform.h>
#include <itkMacro.h>

// local includes:
#include "itkInverseTransform.h"


//----------------------------------------------------------------
// itk::RBFTransform
//
// Radial Basis Function transform:
//
// Let
//   A = (u - uc) / Xmax
//   B = (v - vc) / Ymax
//   Ci = (u - ui) / Xmax
//   Di = (v - vi) / Ymax
//
// where Xmax, Ymax are normalization parameters (typically the width and
// height of the image), uc, vc corresponds to the center or rotation
// of the image expressed in the coordinate system of the mosaic,
// and ui, vi are the mosaic space coordinates of the control points.
//
// The transform is defined as
//   x(u, v) = Xmax * (a0 + a1 * A + a2 * B + F)
//   y(u, v) = Ymax * (b0 + b1 * A + b2 * B + G)
//
// where
//   F = sum(i in [0, k-1], fi * Q(Ci, Di));
//   G = sum(i in [0, k-1], gi * Q(Ci, Di));
//   Q  = r2 * ln(r2)
//   r2 = Ci^2 + Di^2
//
namespace itk
{
class RBFTransform : public Transform<double, 2, 2>
{
public:
  // standard typedefs:
  typedef RBFTransform             Self;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef Transform<double, 2, 2> Superclass;

  // Base inverse transform type:
  typedef Superclass                         InverseTransformType;
  typedef SmartPointer<InverseTransformType> InverseTransformPointer;

  // RTTI:
  itkTypeMacro(RBFTransform, Transform);

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
  TransformPoint(const InputPointType & uv) const;

  // Inverse transformations:
  // If y = Transform(x), then x = BackTransform(y);
  // if no mapping from y to x exists, then an exception is thrown.
  InputPointType
  BackTransformPoint(const OutputPointType & y) const;

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
    return this->m_FixedParameters;
  }

  // virtual:
  void
  SetParameters(const ParametersType & params)
  {
    this->m_Parameters = params;
  }

  // virtual:
  const ParametersType &
  GetParameters() const
  {
    return this->m_Parameters;
  }

  IdentifierType
  GetNumberOfParameters() const override
  {
    return this->m_Parameters.size();
  }

  // virtual: return an inverse of this transform.
  InverseTransformPointer
  GetInverse() const;

  // setup the transform parameters:
  void
  setup(                              // image bounding box expressed in the image space,
                                      // defines transform normalization parameters:
    const OutputPointType & tile_min, // tile space
    const OutputPointType & tile_max, // tile space

    // landmark correspondences:
    const unsigned int      num_pts, // number of pairs
    const InputPointType *  uv,      // mosaic space
    const OutputPointType * xy);     // tile space

  void
  ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & jacobian) const override;

#if 0
    // helper required for numeric inverse transform calculation;
    // evaluate F = T(x), J = dT/dx (another Jacobian):
    void eval(const std::vector<double> & x,
	      std::vector<double> & F,
	      std::vector<std::vector<double> > & J) const;
#endif

  // number of landmark correspondences:
  inline unsigned int
  num_points() const
  {
    return (this->m_FixedParameters.size() - 4) / 2;
  }

  // calculate parameter vector indeces for various transform parameters:
  inline static unsigned int
  index_a(const unsigned int & i)
  {
    return i * 2;
  }

  inline static unsigned int
  index_b(const unsigned int & i)
  {
    return i * 2 + 1;
  }

  inline static unsigned int
  index_f(const unsigned int & i)
  {
    return i * 2 + 6;
  }

  inline static unsigned int
  index_g(const unsigned int & i)
  {
    return i * 2 + 7;
  }

  inline static unsigned int
  index_uv(const unsigned int & i)
  {
    return i * 2 + 4;
  }

  // accessors to the normalization parameters Xmax, Ymax:
  inline const double &
  GetXmax() const
  {
    return this->m_FixedParameters[0];
  }

  inline const double &
  GetYmax() const
  {
    return this->m_FixedParameters[1];
  }

  // accessors to the warp origin expressed in the mosaic coordinate system:
  inline const double &
  GetUc() const
  {
    return this->m_FixedParameters[2];
  }

  inline const double &
  GetVc() const
  {
    return this->m_FixedParameters[3];
  }

  // mosaic space landmark accessor:
  inline const double *
  uv(const unsigned int & i) const
  {
    return &(this->m_FixedParameters[index_uv(i)]);
  }

  // polynomial coefficient accessors:
  inline const double &
  a(const unsigned int & i) const
  {
    return this->m_Parameters[index_a(i)];
  }

  inline const double &
  b(const unsigned int & i) const
  {
    return this->m_Parameters[index_b(i)];
  }

  // radial basis function accessors:
  inline const double &
  f(const unsigned int & i) const
  {
    return this->m_Parameters[index_f(i)];
  }

  inline const double &
  g(const unsigned int & i) const
  {
    return this->m_Parameters[index_g(i)];
  }

  // the Radial Basis Function kernel:
  inline static double
  kernel(const double * uv, const double * uv_i, const double & Xmax, const double & Ymax)
  {
    double U = (uv[0] - uv_i[0]) / Xmax;
    double V = (uv[1] - uv_i[1]) / Ymax;
    double R = U * U + V * V;
    return R == 0 ? 0 : R * log(R);
  }

protected:
  RBFTransform();

  // virtual:
  void
  PrintSelf(std::ostream & s, Indent indent) const;

private:
  // disable default copy constructor and assignment operator:
  RBFTransform(const Self & other);
  const Self &
  operator=(const Self & t);

}; // class RBFTransform

} // namespace itk


#endif // __itkRBFTransform_h
