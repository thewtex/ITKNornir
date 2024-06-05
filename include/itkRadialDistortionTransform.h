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


// File         : itkRadialDistortionTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2005/06/03 10:16
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A radial distortion transform.

#ifndef __itkRadialDistortionTransform_h
#define __itkRadialDistortionTransform_h

// system includes:
#include <iostream>
#include <cassert>

// ITK includes:
#include <itkTransform.h>
#include <itkMacro.h>

// local includes:
#include "itkInverseTransform.h"


//----------------------------------------------------------------
// itk::RadialDistortionTransform
//
// Let
//    A = (a + ta * Rmax - ac)
//    B = (b + tb * Rmax - bc)
// where ta, tb specify percentage of translation along the
// a and b axis respectively, and
//    ac = 0.5 * (a_min + a_max)
//    bc = 0.5 * (b_min + b_max)
// define the center of distortion (specified by a_min, a_max,
// b_min, b_max -- the bounding box of some image).
//
// The transform is defined as
//    x(a, b) = ac + A * S
//    y(a, b) = bc + B * S
//
// where
//    S = sum(n in [0, N - 1], k[n] * pow(R2 / Rmax2, n));
//    R2 = A^2 + B^2
//    Rmax2 = Rmax * Rmax
//
namespace itk
{
template <class TScalar = double, unsigned int N = 2>
class RadialDistortionTransform : public Transform<TScalar, 2, 2>
{
public:
  // standard typedefs:
  typedef RadialDistortionTransform Self;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef Transform<TScalar, 2, 2> Superclass;

  // Base inverse transform type:
  typedef Superclass                         InverseTransformType;
  typedef SmartPointer<InverseTransformType> InverseTransformPointer;

  // static constant for the number of polynomial coefficients:
  itkStaticConstMacro(Nk, unsigned int, N);

  // RTTI:
  itkTypeMacro(RadialDistortionTransform, Transform);

  // macro for instantiation through the object factory:
  itkNewMacro(Self);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  // shortcuts:
  typedef typename Superclass::ParametersType ParametersType;
  typedef typename Superclass::JacobianType   JacobianType;

  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  // virtual:
  OutputPointType
  TransformPoint(const InputPointType & p) const;

  // virtual: Inverse transformations:
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
    return N + 2;
  }

  void
  ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & jacobian) const override;

  // virtual: return an inverse of this transform.
  InverseTransformPointer
  GetInverse() const
  {
    typedef InverseTransform<Self>     InvTransformType;
    typename InvTransformType::Pointer inv = InvTransformType::New();
    inv->SetForwardTransform(this);
    return inv.GetPointer();
  }

  // setup the fixed transform parameters:
  void
  setup( // image bounding box expressed in the physical space:
    const double a_min,
    const double a_max,
    const double b_min,
    const double b_max,

    // normalization parameter (image radius in physical space):
    const double Rmax = 0.0)
  {
    double & ac_ = this->m_FixedParameters[0];
    double & bc_ = this->m_FixedParameters[1];
    double & rmax_ = this->m_FixedParameters[2];

    // center of distortion:
    ac_ = 0.5 * (a_min + a_max);
    bc_ = 0.5 * (b_min + b_max);

    // setup the normalization parameter:
    if (Rmax != 0.0)
    {
      rmax_ = Rmax;
    }
    else
    {
      const double w = a_max - a_min;
      const double h = b_max - b_min;
      rmax_ = sqrt(w * w + h * h) / 2.0;
    }
  }

  // setup the translation parameters:
  void
  setup_translation( // translation is expressed in the physical space:
    const double ta_Rmax = 0.0,
    const double tb_Rmax = 0.0)
  {
    const double & Rmax = this->m_FixedParameters[2];
    assert(Rmax != 0.0);

    // store ta,tb as translation relative to Rmax:
    double & ta = this->m_Parameters[N];
    double & tb = this->m_Parameters[N + 1];
    ta = ta_Rmax / Rmax;
    tb = tb_Rmax / Rmax;
  }

  // helper required by BackTransform:
  // evaluate F = T(x), J = dT/dx (another Jacobian):
  void
  eval(const std::vector<ScalarType> & x, std::vector<ScalarType> & F, std::vector<std::vector<ScalarType>> & J) const;

  // accessors to the fixed normalization parameter:
  inline const double &
  GetRmax() const
  {
    return this->m_FixedParameters[2];
  }

  // generate a mask of shared parameters:
  static void
  setup_shared_params_mask(bool shared, std::vector<bool> & mask)
  {
    mask.assign(N + 2, false);
    for (unsigned int i = 0; i < N; i++)
    {
      mask[i] = shared;
    }
  }

protected:
  RadialDistortionTransform();

  // virtual:
  void
  PrintSelf(std::ostream & s, Indent indent) const;

private:
  // disable default copy constructor and assignment operator:
  RadialDistortionTransform(const Self & other);
  const Self &
  operator=(const Self & t);

}; // class RadialDistortionTransform

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkRadialDistortionTransform.hxx"
#endif

#endif // __itkRadialDistortionTransform_h
