#ifndef CMDLINE_H
#define CMDLINE_H

#include <map>
#include "main.h"

class CmdLine
{
public:
    CmdLine(const string &usage);

    void add(const string &name, const string &desc);

    void add(const string &name, const string &desc, const string &value);

    void help() const;

    void parse(int argc, char **argv);

    bool has(const string &name) const;

    string get(const string &name) const;

private:
    std::map<string,string> m_argval;
    std::map<string,string> m_argdesc;
    std::map<string,bool> m_flagval;
    std::map<string,string> m_flagdesc;
    string m_usage;
};

#endif // CMDLINE_H
