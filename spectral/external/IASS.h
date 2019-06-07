/***************************************************************************
 *   Copyright (C) 2007-2008 by Fraunhofer ITWM Kaiserslautern             *
 *                                                                         *
 *   @author Björn Wagner, bjoern.wagner@itwm.fraunhofer.de                *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Lesser General Public            *
 *   License as published by the Free Software Foundation; either          *
 *   version 3.0 of the License, or (at your option) any later version.    *
 *                                                                         *
 *   This library is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston,                 *
 *   MA 02111-1307 USA                                                     *
 *                                                                         *
 *   Simple methodes for handling 3d iass datasets.                        *
 *                                                                         *
 *   Note: Color images are not supported.                                 *
 *   Note: The label tag is not supported.                                 *
 *                                                                         *
 *   ToDo: More detailed error handling, e.g. can't open, invalid type ... *
 ***************************************************************************/
#ifndef __IASS_h__
#define __IASS_h__

#ifndef _MSC_VER // use this header only with non microsoft compilers
#	include <stdint.h>
#else
#	include "msc_stdint.h"
#endif

// std headers
#include <string>
#include <iostream>

namespace IASS
{
  /// supported pixel types
  enum pixel_type
  { MONO = 0,                   /*< monochrome, one bit per pixle runlenght
                                    encoded - 0=background, false=foreground */
    GREY_8 = 1,                 /*< greyscale 8 bit per pixel, unsigned
                                    integer */
    GREY_16 = 2,                /*< greyscale 16 bit per pixel, unsigned
                                   integer */
    GREY_32 = 3,                /*< greyscale 32 bit per pixel, unsigned
                                    integer */
    GREY_F = 4,                 /*< greyscale 32 bit per pixel, floating
                                    point */
    COMPLEX_F = 5,              /*< greyscale 64 bit per pixel, complex,
                                    ordered pair of floats */
    INVALID
  };

  /// complex pixel
  struct StComplex
  {
    float real;
    float imaginary;
  };

  /// iass header structure
  struct StHeader
  {
		// constructor
		StHeader();		

		// reset header
		void Reset();

	  /// value triple
		template < typename T > struct StTriple
		{
		  T x;                      /*< x coord */
	  	T y;                      /*< y coord */
  	  T z;                      /*< z coord */
		};

	  typedef StTriple < int16_t >  size_type;     /*< type for size triple */
 		typedef StTriple < intmax_t > stride_type;   /*< type for stride triple */
	 	typedef StTriple < double >   spacing_type;  /*< type for spacing triple */

    pixel_type   type;           /*< pixel/data type */
		size_t       bpp;			 			 /*< bytes per pixel */
    size_type    size;           /*< sample size */
    stride_type  stride;         /*< data arrangement, specifies the distance 
                                     between two adjacent pixels for each
                                     direction (x/y/z) in bytes */
    spacing_type spacing;        /*< dimensions of one pixel in meters */
    std::string  creator;        /*< creator of the sample/file */
    std::string  history;        /*< history of the sample */
  };

	// ----------------------------------------------------------------------------
	// simple io functions
	// ----------------------------------------------------------------------------

  /** Load an iass file.
 	 * @param [out] header iass header
 	 * @param [out] pData data buffer, allocated by function
 	 * @param [in] szFile filename
 	 * @returns true if file was loaded sucessfully, false otherwise
 	 */
  bool Load(StHeader & header, uint8_t * &pData, const char *szFile);

  /** Save 3d data as an iass file.
 	 * @param [in] header iass header
 	 * @param [in] pData data buffer
 	 * @param [in] szFile filename, if the file name ends with .iass.gz the file
                        will be compressed with gzip.
 	 * @returns true if file was written sucessfully, false otherwise
 	 */
  bool Save(const StHeader & header, const uint8_t * pData, const char *szFile);

	// ----------------------------------------------------------------------------
	// advanced io functions
	// ----------------------------------------------------------------------------

  /** Read header of an iass dateset from an input stream.
 	 * @param [out] header iass header
 	 * @param [inout] inputStream input stream.
 	 * @returns true if the header was read sucessfully, false otherwise
 	 */
	bool ReadHeader(StHeader& header, std::istream& inputStream);

  /** Write the header of an iass dateset to an ouput stream.
 	 * @param [in] header iass header
 	 * @param [inout] outputStream ouput stream.
 	 * @returns true if the header was written sucessfully, false otherwise
 	 */
	bool WriteHeader(const StHeader& header, std::ostream& outputStream);

  /** Read the pixel data of an iass dateset from an input stream.	
 	 * @param [in] pData data buffer
 	 * @param [out] header iass header
 	 * @param [inout] inputStream input stream.
 	 * @returns true if the data was read sucessfully, false otherwise
 	 */
	bool ReadData(uint8_t* pData, const StHeader& header, std::istream& inputStream);

  /** Write the pixel data of an iass dateset to an output stream.	
 	 * @param [out] pData data buffer
 	 * @param [out] header iass header
 	 * @param [inout] outputStream input stream.
 	 * @returns true if the data was written sucessfully, false otherwise
 	 */
	bool WriteData(const uint8_t* pData, const StHeader& header, std::ostream& outputStream);
}

#endif
