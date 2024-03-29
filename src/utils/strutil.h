﻿/**********************************************************************************
*  MIT License                                                                    *
*                                                                                 *
*  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
*                                                                                 *
*  Geodesy tools for conversions between (historical) datums                      *
*  (c) Chris Veness 2005-2019                                                     *
*  www.movable-type.co.uk/scripts/latlong-convert-coords.html                     *
*  www.movable-type.co.uk/scripts/geodesy-library.html#latlon-ellipsoidal-datum   *
*                                                                                 *
*  Permission is hereby granted, free of charge, to any person obtaining a copy   *
*  of this software and associated documentation files (the "Software"), to deal  *
*  in the Software without restriction, including without limitation the rights   *
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      *
*  copies of the Software, and to permit persons to whom the Software is          *
*  furnished to do so, subject to the following conditions:                       *
*                                                                                 *
*  The above copyright notice and this permission notice shall be included in all *
*  copies or substantial portions of the Software.                                *
*                                                                                 *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       *
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  *
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  *
*  SOFTWARE.                                                                      *
***********************************************************************************/
#ifndef STRUTIL_H
#define STRUTIL_H

#include <string>
#include <vector>

namespace geodesy
{
   class strutil
   {
   public:
      static std::string strip(const std::string& str);
      static bool start_with(const std::string& str, const std::string& prefix);
      static bool ends_with(const std::string& str, const std::string& suffix);

      static std::vector<std::string> split(const std::string& str, wchar_t sep);
      static std::vector<std::string> split_filter_empty(const std::string& str, wchar_t sep);
      static std::vector<std::string> split_regex(const std::string& str, const std::string& sep);

      static std::string padLeft(const std::string& data, const size_t& totalWidth, const char& padding);
      static std::string padLeft(const std::string& data, const size_t& totalWidth, const std::string& padding);
      static std::string padRight(const std::string& data, const size_t& totalWidth, const char& padding);
      static std::string padRight(const std::string& data, const size_t& totalWidth, const std::string& padding);
   };
}


#endif // STRUTIL_H
