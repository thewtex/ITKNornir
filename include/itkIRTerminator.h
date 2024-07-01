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


// File         : itk_terminator.hxx
// Author       : Pavel A. Koshevoy
// Created      : Sun Sep 24 18:16:00 MDT 2006
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : a ITK filter terminator convenience class

#ifndef ITK_TERMINATOR_HXX_
#define ITK_TERMINATOR_HXX_

// local includes:
#include "IRTerminator.h"

#ifdef USE_ITK_TERMINATORS
#define itk_terminator_t the_terminator_t


//----------------------------------------------------------------
// terminator_t
//
template <typename process_t>
class terminator_t : public the_terminator_t
{
public:
  terminator_t(process_t * proc):
    itk_terminator_t(proc->GetNameOfClass()),
    process_(proc)
  {}
  
  // virtual:
  void terminate()
  {
    process_->AbortGenerateDataOn();
    itk_terminator_t::terminate();
  }
  
private:
  typename process_t::Pointer process_;
};


#endif // USE_ITK_TERMINATORS
#endif // ITK_TERMINATOR_HXX_
