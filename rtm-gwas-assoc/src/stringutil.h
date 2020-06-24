#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include "types.h"

class Token
{
public:
    using size_type = isize_t;

    Token() = default;

    Token(const char *str, isize_t len) : str_(str), len_(len) {}

    isize_t length() const { return len_; }

    bool empty() const { return len_ == 0; }

    char operator[](isize_t pos) const { return str_[pos]; }

    const char* data() const { return str_; }

    std::string str() const { return { str_, (std::string::size_type) len_ }; }

    bool operator==(const Token &t) const { return compare(t) == 0; }

    bool operator!=(const Token &t) const { return compare(t) != 0; }

    bool operator<(const Token &t) const { return compare(t) < 0; }

    bool operator>(const Token &t) const { return compare(t) > 0; }

    int compare(const Token &t) const;

private:
    const char *str_ = nullptr;
    isize_t len_ = 0;
};

inline isize_t length(const Token &t)
{
    return t.length();
}

bool contains(const std::string &a, const std::string &b);

isize_t index(const std::string &a, const std::string &b);

bool starts_with(const std::string &s, const std::string &prefix);

bool ends_with(const std::string &s, const std::string &suffix);

std::vector<std::string> split(const std::string &s, const std::string &sep);

void split(const std::string &s, const std::string &sep, std::vector<std::string> &v);

void split(const std::string &s, const std::string &sep, std::vector<Token> &v);

std::string join(const std::vector<std::string> &v, const std::string &sep);

std::string to_upper(const std::string &s);

std::string to_lower(const std::string &s);

#endif // STRINGUTIL_H
