#include "util.h"


bool starts_with(const std::string &s1, const std::string &s2)
{
    return s1.size() >= s2.size() && s1.compare(0, s2.size(), s2) == 0;
}

bool ends_with(const std::string &s1, const std::string &s2)
{
    return s1.size() >= s2.size() && s1.compare(s1.size() - s2.size(), s2.size(), s2) == 0;
}

std::vector<std::string> split(const std::string &str, const std::string &sep)
{
    std::vector<std::string> v;

    auto beg = str.data();
    auto i = str.find_first_not_of(sep);
    auto j = str.find_first_of(sep, i);

    while (j != std::string::npos) {
        v.emplace_back(beg + i, j - i);
        i = str.find_first_not_of(sep, j);
        j = str.find_first_of(sep, i);
    }

    if (i != std::string::npos)
        v.emplace_back(beg + i, str.size() - i);

    return  v;
}

std::string join(const std::vector<std::string> &vs, const std::string &sep)
{
    std::string s;

    auto itr = vs.begin();
    if (itr != vs.end()) {
        s += *itr;
        ++itr;
    }

    for (; itr != vs.end(); ++itr) {
        s += sep;
        s += *itr;
    }

    return s;
}
