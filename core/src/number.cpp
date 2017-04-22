#include <cstring>
#include <limits>
#include <iostream>
#include "number.h"

template<>
int number(const char *str, bool *ok)
{
    errno = 0;
    if ( ok ) *ok = true;

    char *ptr;
    auto ret = std::strtol(str, &ptr, 10);

    if (ptr == str || ptr != str + std::strlen(str) || errno == ERANGE ||
            ret < std::numeric_limits<int>::min()  ||
            ret > std::numeric_limits<int>::max()) {
        if ( ok ) *ok = false;
        ret = std::numeric_limits<int>::min();
        std::cerr << "ERROR: failed to convert string to int: " << str << "\n";
    }

    return ret;
}

template<>
long number(const char *str, bool *ok)
{
    errno = 0;
    if ( ok ) *ok = true;

    char *ptr;
    auto ret = std::strtol(str, &ptr, 10);

    if (ptr == str || ptr != str + std::strlen(str) || errno == ERANGE) {
        if ( ok ) *ok = false;
        ret = std::numeric_limits<long>::min();
        std::cerr << "ERROR: failed to convert string to long: " << str << "\n";
    }

    return ret;
}

template<>
long long number(const char *str, bool *ok)
{
    errno = 0;
    if ( ok ) *ok = true;

    char *ptr;
    auto ret = std::strtoll(str, &ptr, 10);

    if (ptr == str || ptr != str + std::strlen(str) || errno == ERANGE) {
        if ( ok ) *ok = false;
        ret = std::numeric_limits<long long>::min();
        std::cerr << "ERROR: failed to convert string to long long: " << str << "\n";
    }

    return ret;
}

template<>
double number(const char *str, bool *ok)
{
    errno = 0;
    if ( ok ) *ok = true;

    char *ptr;
    auto ret = std::strtod(str, &ptr);

    if (ptr == str || ptr != str + std::strlen(str) || errno == ERANGE) {
        if ( ok ) *ok = false;
        ret = std::numeric_limits<double>::quiet_NaN();
        std::cerr << "ERROR: failed to convert string to double: " << str << "\n";
    }

    return ret;
}

template<>
int number(const std::string &str, bool *ok)
{
    return number<int>(str.c_str(), ok);
}

template<>
long number(const std::string &str, bool *ok)
{
    return number<long>(str.c_str(), ok);
}

template<>
long long number(const std::string &str, bool *ok)
{
    return number<long long>(str.c_str(), ok);
}

template<>
double number(const std::string &str, bool *ok)
{
    return number<double>(str.c_str(), ok);
}
