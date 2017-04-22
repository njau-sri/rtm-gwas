#ifndef NUMBER_H
#define NUMBER_H

#include <string>

template<typename T>
T number(const char *str, bool *ok = nullptr);

template<typename T>
T number(const std::string &str, bool *ok = nullptr);

#endif // NUMBER_H
