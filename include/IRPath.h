// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: nil -*-
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


// File         : IRPath.h
// Author       : Bradley C. Grimm
// Created      : 2009/09/15 11:24 AM
// Copyright    : (C) 2009 University of Utah
// License      : GPLv2
// Description  : Helpful path manipulation tools.

#ifndef __IR_PATH_HXX__
#define __IR_PATH_HXX__

#if defined(WIN32)
  #pragma warning ( disable : 4996 )
#endif
#include "dirent.h"
#if defined(WIN32)
  #pragma warning ( default : 4996 )
#endif

#include <itkIRText.h>

class IRPath
{
public:
  static the_text_t CleanPath(const the_text_t &path)
  {
    the_text_t directory = path;
    CleanSlashes(directory);
    if ( !directory.is_empty() && !directory.match_tail("/") )
      directory += the_text_t("/");
    return directory;
  }

  static the_text_t DirectoryFromPath(const the_text_t &path)
  {
    the_text_t directory = path;
    CleanSlashes(directory);

    directory = directory.splitAt('/', std::numeric_limits<unsigned int>::max())[0];
    if ( !directory.is_empty() && !directory.match_tail("/") )
      return directory + "/";
    else
      return directory;
  }

  static the_text_t FilenameFromPath(the_text_t &path)
  {
    CleanSlashes(path);
    
    the_text_t directory;
    if ( path.contains('/') )
      directory = path.splitAt('/', std::numeric_limits<unsigned int>::max())[1];
    else
      directory = path;
    return directory;
  }

  static void CleanSlashes(the_text_t &path)
  {
    if ( path.contains('\\') )
      path.replace('\\','/');
  }

  static std::list<the_text_t> DirectoryContents(const the_text_t &path)
  {
    std::list<the_text_t> files;
  
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(path.text())) == NULL) {
      cout << "Error(" << errno << ") opening " << path << endl;
      return files;
    }
  
    while ((dirp = readdir(dp)) != NULL) {
      files.push_back( dirp->d_name );
    }
    closedir(dp);

    // Remove .
    files.pop_front();

    // Remove ..
    files.pop_front();

    return files;
  }

};

#endif


