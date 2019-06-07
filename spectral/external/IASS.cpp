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
 ***************************************************************************/
// local headers
#include "IASS.h"

#ifdef _IASS_COMPRESSION_SUPPORT // support compressed iass files
#include "gzstream.h"
#endif

// std headers
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <memory>
#include <algorithm>
#include <cctype>
#include <cstring>

namespace IASS
{
	StHeader::StHeader()
		: type(INVALID), bpp(0), size(), stride(), spacing(), creator(""), history("")
	{
	
	}

	void StHeader::Reset()
	{
    type = INVALID;
		bpp = 0;
    size.x = 0;
    size.y = 0;
    size.z = 0;
    stride.x = 0;
    stride.y = 0;
    stride.z = 0;
    spacing.x = 0.0;
    spacing.y = 0.0;
    spacing.z = 0.0;
    creator = "";
    history = "";
	}

  // simple buffer class for autodeletion of allocated memory
  class CBuffer
  {
  public:
    CBuffer():m_pData(0)
    {
    }

    CBuffer(size_t size):m_pData(0)
    {
      m_pData = new uint8_t[size];
    }

    ~CBuffer()
    {
      delete[]m_pData;
    }

    void Allocate(size_t size)
    {
      delete[]m_pData;
      m_pData = new uint8_t[size];
    }

    uint8_t *Ptr()
    {
      return m_pData;
    }

    uint8_t *Release()
    {
      uint8_t *pTmp = m_pData;

      m_pData = 0;
      return pTmp;
    }

  private:
    uint8_t * m_pData;

    // no copying
  CBuffer(const CBuffer & b):m_pData(0)
    {
    };
    void operator=(const CBuffer & b)
    {
    };
  };

	bool ReadHeader(StHeader& header, std::istream& inputStream)
	{
		header.Reset();

		if(!inputStream.good())
			return false;

      std::string sLine;

	    // read magic
      getline(inputStream, sLine);
      if (!(sLine.substr(0, 9) == "SVstatmat" || sLine.substr(0, 4) == "a4iL")
          || inputStream.fail())
        return false;

      // read header
      do
      {
        getline(inputStream, sLine);

        if (inputStream.fail())
          return false;

        if (sLine.substr(0, 11) == "# SPACING: ")
        {
          std::istringstream in(sLine.substr(11, sLine.length()));
          in >> header.spacing.x >> header.spacing.y >> header.spacing.z;
        }
        else if (sLine.substr(0, 11) == "# CREATOR: ")
          header.creator = sLine.substr(11, sLine.length());

        else if (sLine.substr(0, 11) == "# HISTORY: ")
          header.history = sLine.substr(11, sLine.length());
        else if (sLine.substr(0, 8) == "# TYPE: ")
        {
          if (isdigit(sLine[8]))
          {
            std::string sTmp = sLine.substr(8, sLine.length());
            std::istringstream input_helper(sTmp);
            int type;

            input_helper >> type;
            header.type = static_cast < pixel_type > (type);

            if (header.type < 0 || header.type >= INVALID)
              return false;
          }
          else
          {
            if (sLine.substr(8, sLine.length()) == "MONO")
              header.type = MONO;
            else if (sLine.substr(8, sLine.length()) == "GREY_8")
              header.type = GREY_8;
            else if (sLine.substr(8, sLine.length()) == "GREY_16")
              header.type = GREY_16;
            else if (sLine.substr(8, sLine.length()) == "GREY_32")
              header.type = GREY_32;
            else if (sLine.substr(8, sLine.length()) == "GREY_F")
              header.type = GREY_F;
            else if (sLine.substr(8, sLine.length()) == "COMPLEX_F")
              header.type = COMPLEX_F;
            else
              return false;
          }
        }
      }
      while (sLine[0] == '#');

      if (inputStream.fail() || inputStream.eof())
        return false;
      else
      {
        std::istringstream input_helper(sLine);
        input_helper >> header.size.x >> header.size.y >> header.size.z;
        if (inputStream.fail())
          return false;
      }

      if (header.type == MONO)
      	header.bpp = 1;
     	else if (header.type == GREY_8)
        header.bpp = 1;
      else if (header.type == GREY_16)
        header.bpp = 2;
      else if (header.type == GREY_32)
        header.bpp = 4;
      else if (header.type == GREY_F)
        header.bpp = 4;
      else if (header.type == COMPLEX_F)
	      header.bpp = 8;

		return true;
	}

	bool WriteHeader(const StHeader& header, std::ostream& outputStream)
	{
    // write magic
    outputStream << "SVstatmat publicIASS v1" << std::endl;
    if (outputStream.fail())
      return false;

    // write type
    if (header.type == MONO)
      outputStream << "# TYPE: MONO" << std::endl;
    else if (header.type == GREY_8)
      outputStream << "# TYPE: GREY_8" << std::endl;
    else if (header.type == GREY_16)
      outputStream << "# TYPE: GREY_16" << std::endl;
    else if (header.type == GREY_32)
      outputStream << "# TYPE: GREY_32" << std::endl;
    else if (header.type == GREY_F)
      outputStream << "# TYPE: GREY_F" << std::endl;
    else if (header.type == COMPLEX_F)
      outputStream << "# TYPE: COMPLEX_F" << std::endl;
    else
      return false;

    if (outputStream.fail())
      return false;

    // write spacing
    StHeader::spacing_type spacing;
    spacing.x = (header.spacing.x == 0.0) ? 1.0 : header.spacing.x;
    spacing.y = (header.spacing.y == 0.0) ? 1.0 : header.spacing.y;
    spacing.z = (header.spacing.z == 0.0) ? 1.0 : header.spacing.z;

    outputStream << "# SPACING: " << spacing.x << " " << spacing.
      y << " " << spacing.z << std::endl;

    if (outputStream.fail())
      return false;

    // write creator
    if (header.creator.length() != 0)
    {
      outputStream << "# CREATOR: " << header.creator << std::endl;
      if (outputStream.fail())
        return false;
      }

      // write history
      if (header.creator.length() != 0)
      {
        outputStream << "# HISTORY: " << header.history << std::endl;
        if (outputStream.fail())
          return false;
      }

      // write size
      outputStream << header.size.x << " " << header.size.y << " " << header.size.
        z << std::endl;
      if (outputStream.fail())
        return false;

			return true;
	}

	bool ReadData(uint8_t* pData, const StHeader& header, std::istream& inputStream)
	{
      if (header.type == MONO)
      {
        // read length of packed data
        std::streamsize nLengthRleStream;

        inputStream >> nLengthRleStream;

        char cLinefeed;

        inputStream.get(cLinefeed);

        if (inputStream.fail())
          return false;

        if (cLinefeed != '\n')
          return false;

        CBuffer rleBuffer(nLengthRleStream);
        size_t nSliceSize = header.size.y * header.size.z;
        CBuffer sliceBuffer(nSliceSize);

        // position in the packed data buffer
        std::streamsize nPosRleStream = 0;

        // slice index
        size_t nIndex = 0;

        // position in output (slice) stream
        size_t nPosOutStream = 0;

        // Current pixel value
        uint8_t nCurrValue = 1;

        // read packed data
		inputStream.read(reinterpret_cast < char *>(rleBuffer.Ptr()), nLengthRleStream);

        if (inputStream.fail())
          return false;

        // now start decoding
        size_t nCurrLength = *rleBuffer.Ptr();

        do
        {
          // This loop is entered if another pixel run would exceed the size
          // of the slice buffer
          while (nPosOutStream + nCurrLength >= nSliceSize)
          {
            // How many pixels are left to fill up buffer
            size_t nRestLength = nSliceSize - nPosOutStream;

            // Set these remaining pixels
            memset(sliceBuffer.Ptr() + nPosOutStream,
                   static_cast < uint8_t > (nCurrValue % 2), nRestLength);

            for (uint16_t y = 0; y < header.size.y; ++y)
            {
              uint8_t *p = pData + header.stride.x * nIndex + header.stride.y * y;
              uint8_t *q = sliceBuffer.Ptr() + y * header.size.z;

              for (uint16_t z = 0; z < header.size.z;
                   ++z, p += header.stride.z, ++q)
                *p = *q;
            }

            nCurrLength -= nRestLength;
            nPosOutStream = 0;
            ++nIndex;
          }

          // Set pixels of required length
          if (nCurrLength > 0)
          {
            memset(sliceBuffer.Ptr() + nPosOutStream, static_cast < uint8_t > (nCurrValue % 2), nCurrLength);
            nPosOutStream += nCurrLength;
          }

          nCurrValue++;
          nPosRleStream++;
          nCurrLength = *(rleBuffer.Ptr() + nPosRleStream);
        }
        while (nPosRleStream < nLengthRleStream);
      }
      else
      {
        CBuffer lineBuffer(header.size.z * header.bpp);

        for (uint16_t x = 0; x < header.size.x; ++x)
          for (uint16_t y = 0; y < header.size.y; ++y)
          {
            uint8_t *p = pData + header.stride.x * x + header.stride.y * y;
            uint8_t *q = lineBuffer.Ptr();
            inputStream.read(reinterpret_cast < char *>(lineBuffer.Ptr()),  static_cast<std::streamsize>(header.size.z * header.bpp));

            if (inputStream.fail())
              return false;

            for (uint16_t z = 0; z < header.size.z; ++z, p += header.stride.z, q += header.bpp)
              memcpy(p, q, header.bpp);
          }
      }
		return true;
	}

	bool WriteData(const uint8_t* pData, const StHeader& header, std::ostream& outputStream)
	{
    if (header.type == MONO)
    {
      size_t nDataSize = header.size.x * header.size.y * header.size.z;
		  size_t nMax8bit = std::numeric_limits < uint8_t >::max();
      size_t nNumOf8bitNumbers = nMax8bit + 1;

      // container storing run length encoded data
      CBuffer rleBuffer(nDataSize + 1);

      // encode
      uint8_t nLastValue = 0;
      uint8_t nCurrValue = 0;
      std::streamsize nCounter = 0;

      if (*pData == 0)
      {
        nLastValue = 0;
        *rleBuffer.Ptr() = 0;
        nCounter++;
      }
      else
        nLastValue = 1;

      size_t nChordLength = 0;

      for (uint16_t x = 0; x < header.size.x; ++x)
        for (uint16_t y = 0; y < header.size.y; ++y)
          for (uint16_t z = 0; z < header.size.z; ++z)
          {
            if (pData[x + header.stride.y * y + header.stride.z * z] > 0)
              nCurrValue = 1;
            else
              nCurrValue = 0;

            if (nLastValue == nCurrValue)
              nChordLength++;
            else
            {
              while (nChordLength >= nNumOf8bitNumbers)
              {
                *(rleBuffer.Ptr() + nCounter) = (uint8_t) nMax8bit;
                *(rleBuffer.Ptr() + nCounter + 1) = 0;
                nCounter += 2;
                nChordLength -= nMax8bit;
              }
              *(rleBuffer.Ptr() + nCounter) = (uint8_t) nChordLength;
              nCounter++;
              nChordLength = 1;
              nLastValue = nCurrValue;
            }
          }

      while (nChordLength >= nNumOf8bitNumbers)
      {
        *(rleBuffer.Ptr() + nCounter) = (uint8_t) nMax8bit;
        *(rleBuffer.Ptr() + nCounter + 1) = 0;
        nCounter += 2;
        nChordLength -= nMax8bit;
      }
      *(rleBuffer.Ptr() + nCounter) = (uint8_t) nChordLength;
      nCounter++;

      // write to file
      outputStream << nCounter;
      if (outputStream.fail())
        return false;

      outputStream.put('\n');
      if (outputStream.fail())
        return false;

      // write data
      outputStream.write(reinterpret_cast < char *>(rleBuffer.Ptr()), nCounter);

      if (outputStream.fail())
        return false;
    }
    else
    {
      CBuffer lineBuffer(header.size.z * header.bpp);

      for (uint16_t x = 0; x < header.size.x; ++x)
        for (uint16_t y = 0; y < header.size.y; ++y)
        {
          const uint8_t *q =
            pData + header.stride.x * x + header.stride.y * y;
          uint8_t *p = lineBuffer.Ptr();

          for (uint16_t z = 0; z < header.size.z;
               ++z, q += header.stride.z, p += header.bpp)
            memcpy(p, q, header.bpp);
          outputStream.write(reinterpret_cast < char *>(lineBuffer.Ptr()), static_cast<std::streamsize>(header.size.z * header.bpp));
          if (outputStream.fail())
            return false;
        }
    }
		return true;
	}

  bool Load(StHeader & header, uint8_t * &pData, std::istream& inputStream)
	{
		try
		{
			header.Reset();

 	  	if (!inputStream.good())
				return false;	

	    CBuffer pixelBuffer;

			// read header
			if(!ReadHeader(header, inputStream))
				return false;

	    pixelBuffer.Allocate(header.size.x * header.size.y * header.size.z * header.bpp);
  	  header.stride.x = header.bpp;
    	header.stride.y = header.stride.x * header.size.x;
	    header.stride.z = header.stride.y * header.size.y;

			if(!ReadData(pixelBuffer.Ptr(), header, inputStream))
				return false;

	  	pData = pixelBuffer.Release();
  	 	return true;
		}
    catch(...)
    {
      return false;
    }
	return true;
	}

  bool Load(StHeader & header, uint8_t * &pData, const char *szFile)
  {
    try
    {
			header.Reset();
	
			if(szFile==NULL)
				return false;

#ifdef _IASS_COMPRESSION_SUPPORT // support compressed iass files
			if(igzstream::IsZipped(szFile))
			{
				// open compressed file stream, will be closed when stream 
        // object is destructed
				igzstream inputStream(szFile);
	  	  // check file stream
  	  	if (inputStream.rdbuf()->is_open())
			    return Load(header, pData, inputStream);
			}
			else
			{	
#endif
				// open default file stream, will be closed when stream 
        // object is destructed
				std::ifstream inputStream(szFile, std::ifstream::binary);
	    	if (inputStream.rdbuf()->is_open())
					return Load(header, pData, inputStream);
#ifdef _IASS_COMPRESSION_SUPPORT // support compressed iass files
			}
#endif
    }
    catch(...)
    {
      return false;
    }
	return true;
  }

	bool Save(const StHeader & header, const uint8_t * pData, std::ostream& outputStream)
	{
    try
    {
      if (pData == NULL)
        return false;

      if (!outputStream.good())
        return false;

			if(!WriteHeader(header, outputStream))
				return false;

			return WriteData(pData, header, outputStream);
    }
    catch(...)
    {
      return false;
    }
	}

	bool Save(const StHeader & header, const uint8_t * pData, const char *szFile)
	{
    try
    {
      if (pData == NULL)
        return false;

      if (szFile == NULL)
        return false;

#ifdef _IASS_COMPRESSION_SUPPORT // support compressed iass files
			std::string sFileName(szFile);
			std::transform(sFileName.begin(), sFileName.end(), sFileName.begin(), (int(*)(int)) std::tolower);
			std::string::size_type nPos = sFileName.rfind(".iass.gz");
			if(nPos!=std::string::npos && sFileName.length()-nPos==8)
			{
	      // open compressed file stream
  	   	ogzstream outputStream(szFile);

	      // check file stream
  	    if (!outputStream.rdbuf()->is_open())
    	    return false;

				return Save(header, pData, outputStream);
			}
			else
			{
#endif
	      // open file stream
  	    std::ofstream outputStream(szFile, std::ifstream::binary);

	      // check file stream
  	    if (!outputStream.is_open())
    	    return false;

				return Save(header, pData, outputStream);
#ifdef _IASS_COMPRESSION_SUPPORT // support compressed iass files
			}
#endif
    }
    catch(...)
    {
      return false;
    }
  }
}
