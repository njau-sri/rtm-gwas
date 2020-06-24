#ifndef CFILE_H
#define CFILE_H

#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

// ascii 0, utf 1-5, binary 6, error -1
int check_ascii_file(const std::string &filename);

int read_all_text(const std::string &filename, std::string &content);

int write_all_text(const std::string &content, const std::string &filename);

class CFile
{
public:
    CFile() = default;
    CFile(const char *filename, const char *mode) : fp_(std::fopen(filename, mode)) {}
    CFile(const std::string &filename, const char *mode) : CFile(filename.c_str(), mode) {}

    void open(const char *filename, const char *mode)
    {
        if (fp_) std::fclose(fp_);
        fp_ = std::fopen(filename, mode);
    }

    void open(const std::string &filename, const char *mode)
    {
        open(filename.c_str(), mode);
    }

    CFile(CFile&&) = delete;
    CFile(const CFile&) = delete;
    CFile& operator=(CFile&&) = delete;
    CFile& operator=(const CFile&) = delete;

    ~CFile() { if (fp_) std::fclose(fp_); }

    inline operator std::FILE*() const { return fp_; }

private:
    std::FILE *fp_ = nullptr;
};

class CFileLineReader
{
public:
    CFileLineReader(const char *filename)
        : file_(filename, "rb"), buffer_(new char[buffer_size_])
    {
        if (file_)
            std::setbuf(file_, NULL);
    }

    CFileLineReader(const std::string &filename) : CFileLineReader(filename.c_str()) {}

    bool operator!() const { return !file_; }

    bool read(std::string &str)
    {
        str.clear();

        bool eol = false;
        auto buffer = buffer_.get();

        do {
            if (end_ >= len_) {
                len_ = (int) std::fread(buffer, 1, buffer_size_, file_);
                beg_ = end_ = 0;
            }
            if (len_ > end_) {
                auto ptr = std::memchr(buffer + end_, '\n', (std::size_t) len_ - end_);
                eol = ptr != NULL;
                end_ = eol ? (int) ((char*) ptr - buffer) : len_;
                if (beg_ < end_)
                    str.append(buffer + beg_, buffer + end_);
                beg_ = end_ = end_ + 1;
            }
            else
                eol = true;
        } while (!eol);

        if (!str.empty() && str.back() == '\r')
            str.pop_back();

        return !str.empty() || len_ > 0;
    }

private:
    CFile file_;
    std::unique_ptr<char[]> buffer_;
    static constexpr int buffer_size_ = 1024 * 1024;
    int len_ = 0;
    int beg_ = 0;
    int end_ = 0;
};

#endif // CFILE_H
