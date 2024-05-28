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


// File         : itkNormalizeImageFilterWithMask.h
// Author       : Pavel A. Koshevoy
// Created      : 2006/06/15 17:09
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : An enhanced version of the itk::NormalizeImageFilter
//                adding support for a mask image.

#ifndef __itkNormalizeImageFilterWithMask_h
#define __itkNormalizeImageFilterWithMask_h

// ITK includes:
#include <itkImageToImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkSpatialObject.h>
#include <itkEventObject.h>

// local includes:
#include "itkStatisticsImageFilterWithMask.h"

namespace itk {

/** \class NormalizeImageFilterWithMask 
 * \brief Normalize an image by setting its mean to zero and variance to one.
 *
 * NormalizeImageFilterWithMask shifts and scales an image so that the pixels
 * in the image have a zero mean and unit variance. This filter uses
 * StatisticsImageFilterWithMask to compute the mean and variance of the input
 * and then applies ShiftScaleImageFilter to shift and scale the pixels.
 *
 * NB: since this filter normalizes the data to lie within -1 to 1,
 * integral types will produce an image that DOES NOT HAVE a unit variance.
 * \ingroup MathematicalImageFilters
 */
template<class TInputImage,class TOutputImage>
class ITK_EXPORT NormalizeImageFilterWithMask : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard Self typedef */
  typedef NormalizeImageFilterWithMask Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(NormalizeImageFilterWithMask, ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef typename TInputImage::Pointer InputImagePointer;
  typedef typename TOutputImage::Pointer OutputImagePointer;
  
  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject<TInputImage::ImageDimension> ImageMaskType;
  typedef typename ImageMaskType::Pointer ImageMaskPointer;
  
  /** Set/Get the image mask. */
  itkSetObjectMacro( ImageMask, ImageMaskType );
  itkGetConstObjectMacro( ImageMask, ImageMaskType );

protected:
  NormalizeImageFilterWithMask();

  /** GenerateData. */
  void  GenerateData ();

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

private:
  NormalizeImageFilterWithMask(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename StatisticsImageFilterWithMask<TInputImage>::Pointer m_StatisticsFilter;
  typename ShiftScaleImageFilter<TInputImage,TOutputImage>::Pointer m_ShiftScaleFilter;

  mutable ImageMaskPointer m_ImageMask;
} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizeImageFilterWithMask.hxx"
#endif

#endif
