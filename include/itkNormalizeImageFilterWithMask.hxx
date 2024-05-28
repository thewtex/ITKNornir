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


// File         : itkNormalizeImageFilterWithMask.txx
// Author       : Pavel A. Koshevoy
// Created      : 2006/06/15 17:09
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : An enhanced version of the itk::NormalizeImageFilter
//                adding support for a mask image.

#ifndef _itkNormalizeImageFilterWithMask_txx
#define _itkNormalizeImageFilterWithMask_txx

// ITK includes:
#include <itkImageRegionIterator.h>
#include <itkShiftScaleImageFilter.h>
#include <itkProgressAccumulator.h>
#include <itkCommand.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
NormalizeImageFilterWithMask<TInputImage, TOutputImage>
::NormalizeImageFilterWithMask()
{
  m_StatisticsFilter = 0;
  m_StatisticsFilter = StatisticsImageFilterWithMask<TInputImage>::New();
  m_ShiftScaleFilter = ShiftScaleImageFilter<TInputImage,TOutputImage>::New();
}

template <class TInputImage, class TOutputImage>
void
NormalizeImageFilterWithMask<TInputImage, TOutputImage>
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

template <class TInputImage, class TOutputImage>
void 
NormalizeImageFilterWithMask<TInputImage, TOutputImage>
::GenerateData()
{
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  progress->RegisterInternalFilter(m_StatisticsFilter,.5f);
  progress->RegisterInternalFilter(m_ShiftScaleFilter,.5f);

  // Gather statistics
  
  m_StatisticsFilter->SetInput(this->GetInput());
  m_StatisticsFilter->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  m_StatisticsFilter->SetImageMask(this->m_ImageMask);
  m_StatisticsFilter->Update();

  // Set the parameters for Shift
  m_ShiftScaleFilter->SetShift(-m_StatisticsFilter->GetMean());
  m_ShiftScaleFilter->SetScale(NumericTraits<typename StatisticsImageFilterWithMask<TInputImage>::RealType>::One
                               / m_StatisticsFilter->GetSigma());
  m_ShiftScaleFilter->SetInput(this->GetInput());
  
  m_ShiftScaleFilter->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  m_ShiftScaleFilter->Update();

  // Graft the mini pipeline output to this filters output
  this->GraftOutput(m_ShiftScaleFilter->GetOutput());
}

} // end namespace itk

#endif
