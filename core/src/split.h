#ifndef SPLIT_H
#define SPLIT_H

#include <cstring>
#include <string>

class Token
{
public:
    Token() : s_(nullptr), n_(0) {}

    explicit Token(const std::string &s) : s_(s.data()), n_(s.size()) {}

    Token(const char *s, size_t n) : s_(s), n_(n) {}

    template<typename Iterator>
    Token(Iterator begin, Iterator end) : s_(&*begin), n_(end-begin) {}

    size_t size() const { return n_; }

    bool empty() const { return n_ == 0; }

    char operator[](size_t i) const { return s_[i]; }

    const char* data() const { return s_; }

    explicit operator std::string() const { return std::string(s_, n_); }

    bool operator==(const Token &t) const { return compare(t) == 0; }

    bool operator==(const std::string &s) const { return operator==(Token(s)); }

    bool operator<(const Token &t) const { return compare(t) < 0; }

    bool operator<(const std::string &s) const { return operator<(Token(s)); }

private:
    int compare(const Token &t) const
    {
        return n_ == t.n_ ? std::memcmp(s_, t.s_, n_) : (n_ < t.n_ ? -1 : 1);
    }

private:
    const char *s_;
    size_t n_;
};

template<typename Delimiter, typename Iterator, typename Sequence>
void split(Delimiter delim, Iterator first, Iterator last, Sequence &result)
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

#endif // SPLIT_H
