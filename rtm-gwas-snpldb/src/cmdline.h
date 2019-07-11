#ifndef CMDLINE_H
#define CMDLINE_H

#include <map>
#include <string>

class CmdLine
{
public:

    void add(const std::string &arg, const std::string &msg);

    void add(const std::string &arg, const std::string &msg, const std::string &val);

    void show() const;

    void parse(int argc, char *argv[]);

    bool has(const std::string &arg) const;

    std::string get(const std::string &arg) const;

private:
    void take(const std::string &key, const std::string &val);

private:
    std::map<std::string, std::string> arg_;
    std::map<std::string, std::string> arg_msg_;
    std::map<std::string, bool> flag_;
    std::map<std::string, std::string> flag_msg_;
    std::string usage_;
};

#endif // CMDLINE_H
