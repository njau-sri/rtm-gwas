#ifndef TYPES_H
#define TYPES_H

#include <cctype>
#include <cstddef>
#include <string>
#include <vector>

// %d    int
// %u    unsigned int
// %lld  long long
// %llu  unsigned long long
// %td   ptrdiff_t
// %zu   size_t

using std::int8_t;
using std::int16_t;
using std::int32_t;
using std::int64_t;

using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;

using isize_t = std::ptrdiff_t;
using usize_t = std::size_t;

inline isize_t length(const std::string &s)
{
    return (isize_t) s.length();
}

template<typename T>
inline isize_t length(const std::vector<T> &v)
{
    return (isize_t) v.size();
}

#endif // TYPES_H
