#include <limits>
#include <iostream>
#include "util.h"

template<>
int number(const string &s, bool *ok)
{
    int a = 0;
    size_t p = 0;

    if ( ok ) *ok = true;

    try {
        a = std::stoi(s, &p);
    }
    catch (const std::exception &e) {
        p = 0;
        (void) e;
    }

    if ( s.empty() || p != s.size() ) {
        if ( ok ) *ok = false;
        a = std::numeric_limits<int>::min();
        std::cerr << "ERROR: failed to convert string to INT: " << s << "\n";
    }

    return a;
}

template<>
long number(const string &s, bool *ok)
{
    long a = 0;
    size_t p = 0;

    if ( ok ) *ok = true;

    try {
        a = std::stol(s, &p);
    }
    catch (const std::exception &e) {
        p = 0;
        (void) e;
    }

    if ( s.empty() || p != s.size() ) {
        if ( ok ) *ok = false;
        a = std::numeric_limits<long>::min();
        std::cerr << "ERROR: failed to convert string to LONG: " << s << "\n";
    }

    return a;
}

template<>
long long number(const string &s, bool *ok)
{
    long long a = 0;
    size_t p = 0;

    if ( ok ) *ok = true;

    try {
        a = std::stoll(s, &p);
    }
    catch (const std::exception &e) {
        p = 0;
        (void) e;
    }

    if ( s.empty() || p != s.size() ) {
        if ( ok ) *ok = false;
        a = std::numeric_limits<long long>::min();
        std::cerr << "ERROR: failed to convert string to LONG LONG: " << s << "\n";
    }

    return a;
}

template<>
double number(const string &s, bool *ok)
{
    double a = 0.0;
    size_t p = 0;

    if ( ok ) *ok = true;

    try {
        a = std::stod(s, &p);
    }
    catch (const std::exception &e) {
        p = 0;
        (void) e;
    }

    if ( s.empty() || p != s.size() ) {
        if ( ok ) *ok = false;
        a = std::numeric_limits<double>::quiet_NaN();
        std::cerr << "ERROR: failed to convert string to DOUBLE: " << s << "\n";
    }

    return a;
}
