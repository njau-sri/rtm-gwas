#ifndef STRSPLIT_H
#define STRSPLIT_H


#include <cstring>
#include <string>


class Token
{
public:
    using size_type = std::size_t;

    Token() : str_(nullptr), len_(0) {}

    Token(const char *str, size_type len) : str_(str), len_(len) {}

    size_type size() const { return len_; }

    bool empty() const { return len_ == 0; }

    char operator[](size_type pos) const { return str_[pos]; }

    const char* data() const { return str_; }

    std::string to_string() const { return { str_, len_ }; }

    bool operator==(const Token &t) const { return compare(t) == 0; }

    bool operator!=(const Token &t) const { return compare(t) != 0; }

    bool operator<(const Token &t) const { return compare(t) < 0; }

    bool operator>(const Token &t) const { return compare(t) > 0; }

    int compare(const Token &t) const { return len_ == t.len_ ? std::memcmp(str_, t.str_, len_) : (len_ < t.len_ ? -1 : 1); }

private:
    const char *str_;
    size_type len_;
};

template<typename ContainerT>
void split(const std::string &str, const std::string &sep, ContainerT &vec)
{
    auto dat = str.data();
    auto i = str.find_first_not_of(sep);
    auto j = str.find_first_of(sep, i);

    while (j != std::string::npos) {
        vec.emplace_back(dat + i, j - i);
        i = str.find_first_not_of(sep, j);
        j = str.find_first_of(sep, i);
    }

    if (i != std::string::npos)
        vec.emplace_back(dat + i, str.size() - i);
}

#endif // STRSPLIT_H
