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

