#ifndef CMDLINE_H
#define CMDLINE_H

#include <map>
#include <string>

class CmdLine
{
public:

    void add(const std::string &arg);

    void add(const std::string &arg, const std::string &val);

    int parse(int argc, char *argv[]);

    bool has(const std::string &arg) const;

    std::string get(const std::string &arg) const;

private:
    int take(const std::string &key, const std::string &val);

private:
    std::map<std::string, std::string> arg_;
    std::map<std::string, bool> flag_;
};

#endif // CMDLINE_H
