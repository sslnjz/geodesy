/**********************************************************************************
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
#include "strutil.h"

#include <regex>

using namespace geodesy;

std::vector<std::wstring> strutil::split(const std::wstring& str, wchar_t sep)
{
	std::vector<std::wstring> res;
   auto start = str.begin();
   const auto end = str.end();
   auto next = find(start, end, sep);

	while (next != end)
	{
		if (start < next)
			res.emplace_back(start, next);

		start = next + 1;
		next = find(start, end, sep);
	}

	if (start < next)
		res.emplace_back(start, next);

	return res;
}

std::vector<std::wstring> strutil::split_filter_empty(const std::wstring& str, wchar_t sep)
{
	std::vector<std::wstring> res;
	std::wstring::const_iterator start = str.begin();
	std::wstring::const_iterator end = str.end();
	std::wstring::const_iterator next = find(start, end, sep);

	while (next != end)
	{
		if (start < next)
			res.emplace_back(start, next);

		start = next + 1;
		next = find(start, end, sep);
	}

	if (start < next)
		res.emplace_back(start, next);

	return res;
}


std::vector<std::wstring> strutil::split_regex(const std::wstring& str, const std::wstring& sep)
{
	std::vector<std::wstring> res;
   const std::wregex rgx(sep);

	std::wsregex_token_iterator iter(str.begin(), str.end(), rgx, -1);
   const std::wsregex_token_iterator end;

	while (iter != end) 
	{
		res.push_back(*iter);
		++iter;
	}

	return res;
}

std::wstring strutil::strip(const std::wstring& str)
{
	std::wstring res;

	if (!str.empty())
	{
		const wchar_t* st = str.c_str();
		const wchar_t* ed = st + str.size() - 1;

		while (st <= ed)
		{
			if (!isspace(*st))
				break;
			st++;
		}

		while (ed >= st)
		{
			if (!isspace(*ed))
				break;
			ed--;
		}

		if (ed >= st)
			res.assign(st, ed - st + 1);
	}

	return res;
}

bool strutil::start_with(const std::wstring& str, const std::wstring& prefix)
{
   const size_t prefix_len = prefix.size();

	if (str.size() < prefix_len)
		return false;

	for (size_t i = 0; i < prefix_len; i++)
	{
		if (str[i] != prefix[i])
			return false;
	}

	return true;
}

