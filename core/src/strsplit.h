#ifndef STRSPLIT_H
#define STRSPLIT_H

#include <cstring>
#include <string>

class Token
{
public:
    using size_type = std::string::size_type;

    Token() : s_(nullptr), l_(0) {}

    explicit Token(const std::string &s) : s_(s.data()), l_(s.size()) {}

    Token(const char *s, size_type l) : s_(s), l_(l) {}

    template<typename Iterator>
    Token(Iterator begin, Iterator end) : s_(&*begin), l_(end-begin) {}

    size_type size() const { return l_; }

    bool empty() const { return l_ == 0; }

    char operator[](size_type i) const { return s_[i]; }

    const char* data() const { return s_; }

    explicit operator std::string() const { return std::string(s_, l_); }

    bool operator==(const Token &t) const { return compare(t) == 0; }

    bool operator==(const std::string &s) const { return operator==(Token(s)); }

    bool operator<(const Token &t) const { return compare(t) < 0; }

    bool operator<(const std::string &s) const { return operator<(Token(s)); }

private:
    int compare(const Token &t) const
    {
        return l_ == t.l_ ? std::memcmp(s_, t.s_, l_) : (l_ < t.l_ ? -1 : 1);
    }

private:
    const char *s_;
    size_type l_;
};

template<typename Delimiter, typename Iterator, typename Sequence>
void strsplit(Delimiter delim, Iterator first, Iterator last, Sequence &result)
{
    while (first != last && delim(*first))
        ++first;

    auto second = first;

    while (second != last) {
        if (delim(*second)) {
            result.emplace_back(first, second);
            while (++second != last && delim(*second))
                continue;
            first = second;
        }
        else {
            ++second;
        }
    }

    if (first != second)
        result.emplace_back(first, second);
}

#endif // STRSPLIT_H
