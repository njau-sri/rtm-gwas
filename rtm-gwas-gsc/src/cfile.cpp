#include "cfile.h"

#include <memory>

#include "print.h"

int check_ascii_file(const std::string &filename)
{
    CFile file(filename, "rb");

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return -1;
    }

    // Unicode Transformation Format (UTF)
    // https://www.unicode.org/faq/utf_bom.html
    // Byte Order Mark (BOM)
    // 1    00 00 FE FF    UTF-32 BE
    // 2    FF FE 00 00    UTF-32 LE
    // 3    FE FF          UTF-16 BE
    // 4    FF FE          UTF-16 LE
    // 5    EF BB BF       UTF-8 with BOM

    int a = std::fgetc(file);

    // empty
    if (a == EOF)
        return 0;

    int b = std::fgetc(file);

    // 1 byte
    if (b == EOF)
        return a == 0x00 ? 6 : 0;

    int c = std::fgetc(file);

    // 2 bytes
    if (c == EOF) {
        if (a == 0xFE && b == 0xFF)
            return 3;
        if (a == 0xFF && b == 0xFE)
            return 4;
        if (a == 0x00 || b == 0x00)
            return 6;
        return 0;
    }

    int d = std::fgetc(file);

    // 3 bytes
    if (d == EOF) {
        if (a == 0xEF && b == 0xBB && c == 0xBF)
            return 5;
        if (a == 0xFE && b == 0xFF)
            return c == 0x00 ? 6 : 3;
        if (a == 0xFF && b == 0xFE)
            return c == 0x00 ? 6 : 4;
        if (a == 0x00 || b == 0x00 || c == 0x00)
            return 6;
        return 0;
    }

    // >= 4 bytes

    if (a == 0x00 && b == 0x00 && c == 0xFE && d == 0xFF)
        return 1;

    if (a == 0xFF && b == 0xFE && c == 0x00 && d == 0x00)
        return 2;

    if (a == 0xFE && b == 0xFF)
        return 3;

    if (a == 0xFF && b == 0xFE)
        return 4;

    if (a == 0xEF && b == 0xBB && c == 0xBF)
        return 5;

    for (int i = 0; i < 8192; ++i) {
        if ((c = std::fgetc(file)) == EOF)
            return 0;
        if (c == 0x00)
            return 6;
    }

    return 0;
}

int read_all_text(const std::string &filename, std::string &content)
{
    CFile file(filename, "r");

    if (!file) {
        eprint("ERROR: can't open file for reading: %s\n", filename);
        return 1;
    }

    static const std::size_t buffer_size = 1024 * 1024;
    std::unique_ptr<char[]> buffer(new char[buffer_size]);

    std::size_t n = 0;
    while ((n = std::fread(buffer.get(), 1, buffer_size, file)))
        content.append(buffer.get(), n);

    return 0;
}

int write_all_text(const std::string &content, const std::string &filename)
{
    CFile file(filename, "w");

    if (!file) {
        eprint("ERROR: can't open file for writing: %s\n", filename);
        return 1;
    }

    if (std::fwrite(content.c_str(), 1, content.size(), file) != content.size()) {
        eprint("ERROR: not all data was written: %s\n", filename);
        return 1;
    }

    return 0;
}
