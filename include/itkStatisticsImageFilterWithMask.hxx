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


// File         : itkStatisticsImageFilterWithMask.txx
// Author       : Pavel A. Koshevoy
// Created      : 2006/06/15 17:09
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : An enhanced version of the itk::StatisticsImageFilter
//                adding support for a mask image.

#ifndef _itkStatisticsImageFilterWithMask_txx
#define _itkStatisticsImageFilterWithMask_txx

// ITK includes:
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNumericTraits.h>
#include <itkProgressReporter.h>

namespace itk {

template<class TInputImage>
StatisticsImageFilterWithMask<TInputImage>
::StatisticsImageFilterWithMask(): m_ThreadSum(1), m_SumOfSquares(1), m_Count(1), m_ThreadMin(1), m_ThreadMax(1)
{
  // first output is a copy of the image, DataObject created by
  // superclass
  //
  // allocate the data objects for the outputs which are
  // just decorators around pixel types
  for (int i=1; i < 3; ++i)
    {
    typename PixelObjectType::Pointer output
      = static_cast<PixelObjectType*>(this->MakeOutput(i).GetPointer());
    this->ProcessObject::SetNthOutput(i, output.GetPointer());
    }
  // allocate the data objects for the outputs which are
  // just decorators around real types
  for (int i=3; i < 7; ++i)
    {
    typename RealObjectType::Pointer output
      = static_cast<RealObjectType*>(this->MakeOutput(i).GetPointer());
    this->ProcessObject::SetNthOutput(i, output.GetPointer());
    }

  this->GetMinimumOutput()->Set( NumericTraits<PixelType>::max() );
  this->GetMaximumOutput()->Set( NumericTraits<PixelType>::NonpositiveMin() );
  this->GetMeanOutput()->Set( NumericTraits<RealType>::max() );
  this->GetSigmaOutput()->Set( NumericTraits<RealType>::max() );
  this->GetVarianceOutput()->Set( NumericTraits<RealType>::max() );
  this->GetSumOutput()->Set( NumericTraits<RealType>::Zero );
}


template<class TInputImage>
DataObject::Pointer
StatisticsImageFilterWithMask<TInputImage>
::MakeOutput(unsigned int output)
{
  switch (output)
    {
   case 0:
      return static_cast<DataObject*>(TInputImage::New().GetPointer());
      break;
    case 1:
      return static_cast<DataObject*>(PixelObjectType::New().GetPointer());
      break;
    case 2:
      return static_cast<DataObject*>(PixelObjectType::New().GetPointer());
      break;
    case 3:
    case 4:
    case 5:
    case 6:
      return static_cast<DataObject*>(RealObjectType::New().GetPointer());
      break;
    default:
      // might as well make an image
      return static_cast<DataObject*>(TInputImage::New().GetPointer());
      break;
    }
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::PixelObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMinimumOutput()
{
  return static_cast<PixelObjectType*>(this->ProcessObject::GetOutput(1));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::PixelObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMinimumOutput() const
{
  return static_cast<const PixelObjectType*>(this->ProcessObject::GetOutput(1));
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::PixelObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMaximumOutput()
{
  return static_cast<PixelObjectType*>(this->ProcessObject::GetOutput(2));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::PixelObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMaximumOutput() const
{
  return static_cast<const PixelObjectType*>(this->ProcessObject::GetOutput(2));
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMeanOutput()
{
  return static_cast<RealObjectType*>(this->ProcessObject::GetOutput(3));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetMeanOutput() const
{
  return static_cast<const RealObjectType*>(this->ProcessObject::GetOutput(3));
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetSigmaOutput()
{
  return static_cast<RealObjectType*>(this->ProcessObject::GetOutput(4));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetSigmaOutput() const
{
  return static_cast<const RealObjectType*>(this->ProcessObject::GetOutput(4));
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetVarianceOutput()
{
  return static_cast<RealObjectType*>(this->ProcessObject::GetOutput(5));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetVarianceOutput() const
{
  return static_cast<const RealObjectType*>(this->ProcessObject::GetOutput(5));
}


template<class TInputImage>
typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetSumOutput()
{
  return static_cast<RealObjectType*>(this->ProcessObject::GetOutput(6));
}

template<class TInputImage>
const typename StatisticsImageFilterWithMask<TInputImage>::RealObjectType*
StatisticsImageFilterWithMask<TInputImage>
::GetSumOutput() const
{
  return static_cast<const RealObjectType*>(this->ProcessObject::GetOutput(6));
}

template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    InputImagePointer image =
      const_cast< typename Superclass::InputImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::AllocateOutputs()
{
  // Pass the input through as the output
  InputImagePointer image =
    const_cast< TInputImage * >( this->GetInput() );
  this->GraftOutput( image );

  // Nothing that needs to be allocated for the remaining outputs
}

template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::BeforeThreadedGenerateData()
{
  int numberOfThreads = this->GetNumberOfThreads();

  // Resize the thread temporaries
  m_Count.SetSize(numberOfThreads);
  m_SumOfSquares.SetSize(numberOfThreads);
  m_ThreadSum.SetSize(numberOfThreads);
  m_ThreadMin.SetSize(numberOfThreads);
  m_ThreadMax.SetSize(numberOfThreads);
  
  // Initialize the temporaries
  m_Count.Fill(NumericTraits<long>::Zero);
  m_ThreadSum.Fill(NumericTraits<RealType>::Zero);
  m_SumOfSquares.Fill(NumericTraits<RealType>::Zero);
  m_ThreadMin.Fill(NumericTraits<PixelType>::max());
  m_ThreadMax.Fill(NumericTraits<PixelType>::NonpositiveMin());
}

template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::AfterThreadedGenerateData()
{
  int i;
  long count;
  RealType sumOfSquares;
    
  int numberOfThreads = this->GetNumberOfThreads();

  PixelType minimum;
  PixelType maximum;
  RealType  mean;
  RealType  sigma;
  RealType  variance;
  RealType  sum;
  
  sum = sumOfSquares = NumericTraits<RealType>::Zero;
  count = 0;

  // Find the min/max over all threads and accumulate count, sum and
  // sum of squares
  minimum = NumericTraits<PixelType>::max();
  maximum = NumericTraits<PixelType>::NonpositiveMin();
  for( i = 0; i < numberOfThreads; i++)
    {
    count += m_Count[i];
    sum += m_ThreadSum[i];
    sumOfSquares += m_SumOfSquares[i];

    if (m_ThreadMin[i] < minimum)
      {
      minimum = m_ThreadMin[i];
      }
    if (m_ThreadMax[i] > maximum)
      {
      maximum = m_ThreadMax[i];
      }
    }
  // compute statistics
  mean = sum / static_cast<RealType>(count);

  // unbiased estimate
  variance = (sumOfSquares - (sum*sum / static_cast<RealType>(count)))
    / (static_cast<RealType>(count) - 1);
  sigma = sqrt(variance);

  // Set the outputs
  this->GetMinimumOutput()->Set( minimum );
  this->GetMaximumOutput()->Set( maximum );
  this->GetMeanOutput()->Set( mean );
  this->GetSigmaOutput()->Set( sigma );
  this->GetVarianceOutput()->Set( variance );
  this->GetSumOutput()->Set( sum );
}

template<class TInputImage>
void
StatisticsImageFilterWithMask<TInputImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread,
                       int threadId) 
{
  RealType realValue;
  PixelType value;
  ImageRegionConstIteratorWithIndex<TInputImage> it (this->GetInput(), outputRegionForThread);
  
  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  
  const TInputImage * image = this->GetInput();
  
  // do the work
  while (!it.IsAtEnd())
    {
    bool skip = false;
    if (m_ImageMask)
      {
	// test the mask:
	PointType point;
	image->TransformIndexToPhysicalPoint( it.GetIndex(), point );
	skip = !m_ImageMask->IsInside( point );
      }
    
    if (!skip)
      {
      value = it.Get();
      realValue = static_cast<RealType>( value );
      if (value < m_ThreadMin[threadId])
        {
	m_ThreadMin[threadId] = value;
        }
      if (value > m_ThreadMax[threadId])
        {
	m_ThreadMax[threadId] = value;
        }
      
      m_ThreadSum[threadId] += realValue;
      m_SumOfSquares[threadId] += (realValue * realValue);
      m_Count[threadId]++;
      }
    
    ++it;
    progress.CompletedPixel();
    }
}

template <class TImage>
void 
StatisticsImageFilterWithMask<TImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Minimum: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(this->GetMinimum()) << std::endl;
  os << indent << "Maximum: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(this->GetMaximum()) << std::endl;
  os << indent << "Sum: "      << this->GetSum() << std::endl;
  os << indent << "Mean: "     << this->GetMean() << std::endl;
  os << indent << "Sigma: "    << this->GetSigma() << std::endl;
  os << indent << "Variance: " << this->GetVariance() << std::endl;
}


}// end namespace itk
#endif
