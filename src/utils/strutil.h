
#ifndef STRUTIL_H
#define STRUTIL_H

#include <string>
#include <vector>

namespace geodesy
{
   class strutil
   {
   public:
      static std::wstring strip(const std::wstring& str);
      static bool start_with(const std::wstring& str, const std::wstring& prefix);

      static std::vector<std::wstring> split(const std::wstring& str, wchar_t sep);
      static std::vector<std::wstring> split_filter_empty(const std::wstring& str, wchar_t sep);
      static std::vector<std::wstring> split_regex(const std::wstring& str, const std::wstring& sep);
   };
}


#endif // STRUTIL_H
