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


// File         : itkStatisticsImageFilterWithMask.h
// Author       : Pavel A. Koshevoy
// Created      : 2006/06/15 17:09
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : An enhanced version of the itk::StatisticsImageFilter
//                adding support for a mask image.

#ifndef __itkStatisticsImageFilterWithMask_h
#define __itkStatisticsImageFilterWithMask_h

// ITK includes:
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>
#include <itkArray.h>
#include <itkSimpleDataObjectDecorator.h>
#include <itkSpatialObject.h>


namespace itk {

/** \class StatisticsImageFilterWithMask 
 * \brief Compute min. max, variance and mean of an Image.
 *
 * StatisticsImageFilterWithMask computes the minimum, maximum, sum, mean, variance
 * sigma of an image.  The filter needs all of its input image.  It
 * behaves as a filter with an input and output. Thus it can be inserted
 * in a pipline with other filters and the statistics will only be
 * recomputed if a downstream filter changes.
 *
 * The filter passes its input through unmodified.  The filter is
 * threaded. It computes statistics in each thread then combines them in
 * its AfterThreadedGenerate method.
 *
 * \ingroup MathematicalStatisticsImageFilters
 */
template<class TInputImage>
class ITK_EXPORT StatisticsImageFilterWithMask : 
    public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard Self typedef */
  typedef StatisticsImageFilterWithMask Self;
  typedef ImageToImageFilter<TInputImage,TInputImage>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(StatisticsImageFilterWithMask, ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef typename TInputImage::Pointer InputImagePointer;

  typedef typename TInputImage::RegionType RegionType ;
  typedef typename TInputImage::SizeType SizeType ;
  typedef typename TInputImage::IndexType IndexType ;
  typedef typename TInputImage::PixelType PixelType ;
  typedef typename TInputImage::PointType PointType ;
  
  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension ) ;

  /** Type to use for computations. */
  typedef typename NumericTraits<PixelType>::RealType RealType;

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  /** Type of DataObjects used for scalar outputs */
  typedef SimpleDataObjectDecorator<RealType>  RealObjectType;
  typedef SimpleDataObjectDecorator<PixelType> PixelObjectType;
  
  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject<TInputImage::ImageDimension> ImageMaskType;
  typedef typename ImageMaskType::Pointer ImageMaskPointer;
  
  /** Set/Get the image mask. */
  itkSetObjectMacro( ImageMask, ImageMaskType );
  itkGetConstObjectMacro( ImageMask, ImageMaskType );
  
  /** Return the computed Minimum. */
  PixelType GetMinimum() const
    { return this->GetMinimumOutput()->Get(); }
  PixelObjectType* GetMinimumOutput();
  const PixelObjectType* GetMinimumOutput() const;
  
  /** Return the computed Maximum. */
  PixelType GetMaximum() const
    { return this->GetMaximumOutput()->Get(); }
  PixelObjectType* GetMaximumOutput();
  const PixelObjectType* GetMaximumOutput() const;

  /** Return the computed Mean. */
  RealType GetMean() const
    { return this->GetMeanOutput()->Get(); }
  RealObjectType* GetMeanOutput();
  const RealObjectType* GetMeanOutput() const;

  /** Return the computed Standard Deviation. */
  RealType GetSigma() const
    { return this->GetSigmaOutput()->Get(); }
  RealObjectType* GetSigmaOutput();
  const RealObjectType* GetSigmaOutput() const;

  /** Return the computed Variance. */
  RealType GetVariance() const
    { return this->GetVarianceOutput()->Get(); }
  RealObjectType* GetVarianceOutput();
  const RealObjectType* GetVarianceOutput() const;

  /** Return the compute Sum. */
  RealType GetSum() const
    { return this->GetSumOutput()->Get(); }
  RealObjectType* GetSumOutput();
  const RealObjectType* GetSumOutput() const;

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  virtual DataObjectPointer MakeOutput(unsigned int idx);

protected:
  StatisticsImageFilterWithMask();
  ~StatisticsImageFilterWithMask(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Pass the input through unmodified. Do this by Grafting in the AllocateOutputs method. */
  void AllocateOutputs();      

  /** Initialize some accumulators before the threads run. */
  void BeforeThreadedGenerateData ();
  
  /** Do final mean and variance computation from data accumulated in threads. */
  void AfterThreadedGenerateData ();
  
  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData (const RegionType& 
                              outputRegionForThread,
                              int threadId) ;

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion(DataObject *data);

private:
  StatisticsImageFilterWithMask(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  mutable ImageMaskPointer m_ImageMask;

  Array<RealType>  m_ThreadSum;
  Array<RealType>  m_SumOfSquares;
  Array<long>      m_Count;
  Array<PixelType> m_ThreadMin;
  Array<PixelType> m_ThreadMax;

} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticsImageFilterWithMask.hxx"
#endif

#endif
