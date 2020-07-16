#ifndef PRINT_H
#define PRINT_H

#include <cstdio>
#include <string>

namespace internal {

    template<typename T>
    inline T argument(T a) { return a; }

    template<typename T>
    inline const T* argument(const std::basic_string<T> &a) { return a.c_str(); }

} // internal

// stdout

template<typename ... Ts>
inline void print(const char *format, const Ts& ... args)
{
    std::fprintf(stdout, format, internal::argument(args) ...);
}

inline void print(const char *str)
{
    std::fputs(str, stdout);
}

// stderr

template<typename ... Ts>
inline void eprint(const char *format, const Ts& ... args)
{
    std::fprintf(stderr, format, internal::argument(args) ...);
}

inline void eprint(const char *str)
{
    std::fputs(str, stderr);
}

// file

template<typename ... Ts>
inline void fprint(std::FILE *stream, const char *format, const Ts& ... args)
{
    std::fprintf(stream, format, internal::argument(args) ...);
}

inline void fprint(std::FILE *stream, const char *str)
{
    std::fputs(str, stream);
}

// string

template<typename ... Ts>
inline std::string sprint(const char *format, const Ts& ... args)
{
    std::string s;

    int n = std::snprintf(nullptr, 0, format, internal::argument(args) ...);

    if (n > 0) {
        s.resize(n);
        std::snprintf(&s[0], n+1, format, internal::argument(args) ...);
    }

    return s;
}

inline std::string sprint(const char *str)
{
    return std::string(str);
}

#endif // PRINT_H
