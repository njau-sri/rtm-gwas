#include "stringutil.h"

#include <cctype>
#include <cstring>
#include <algorithm>

int Token::compare(const Token &t) const
{
    return len_ == t.len_ ? std::memcmp(str_, t.str_, len_) : (len_ < t.len_ ? -1 : 1);
}

bool contains(const std::string &a, const std::string &b)
{
    return a.find(b) != std::string::npos;
}

isize_t index(const std::string &a, const std::string &b)
{
    auto pos = a.find(b);
    return pos != std::string::npos ? (isize_t) pos : -1;
}

bool starts_with(const std::string &s, const std::string &prefix)
{
    return length(s) >= length(prefix) && s.compare(0, length(prefix), prefix) == 0;
}

bool ends_with(const std::string &s, const std::string &suffix)
{
    return length(s) >= length(suffix) && s.compare(length(s) - length(suffix), length(suffix), suffix) == 0;
}

std::vector<std::string> split(const std::string &s, const std::string &sep)
{
    std::vector<std::string> v;

    split(s, sep, v);

    return  v;
}

void split(const std::string &s, const std::string &sep, std::vector<std::string> &v)
{
    auto ptr = s.data();
    auto i = s.find_first_not_of(sep);
    auto j = s.find_first_of(sep, i);

    while (j != std::string::npos) {
        v.emplace_back(ptr + i, j - i);
        i = s.find_first_not_of(sep, j);
        j = s.find_first_of(sep, i);
    }

    if (i != std::string::npos)
        v.emplace_back(ptr + i, s.size() - i);
}

void split(const std::string &s, const std::string &sep, std::vector<Token> &v)
{
    auto ptr = s.data();
    auto i = s.find_first_not_of(sep);
    auto j = s.find_first_of(sep, i);

    while (j != std::string::npos) {
        v.emplace_back(ptr + i, j - i);
        i = s.find_first_not_of(sep, j);
        j = s.find_first_of(sep, i);
    }

    if (i != std::string::npos)
        v.emplace_back(ptr + i, s.size() - i);
}

std::string join(const std::vector<std::string> &v, const std::string &sep)
{
    std::string s;

    auto itr = v.begin();
    if (itr != v.end()) {
        s += *itr;
        ++itr;
    }

    for (; itr != v.end(); ++itr) {
        s += sep;
        s += *itr;
    }

    return s;
}

std::string to_upper(const std::string &s)
{
    auto t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return std::toupper(c); });
    return t;
}

std::string to_lower(const std::string &s)
{
    auto t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return std::tolower(c); });
    return t;
}
