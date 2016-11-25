#include <iomanip>
#include <iostream>
#include "cmdline.h"

CmdLine::CmdLine(const string &usage) : m_usage(usage) {}

void CmdLine::add(const string &name, const string &desc)
{
    m_flagval[name] = false;
    m_flagdesc[name] = desc;
}

void CmdLine::add(const string &name, const string &desc, const string &value)
{
    m_argval[name] = value;
    m_argdesc[name] = desc;
}

void CmdLine::help() const
{
    std::cerr << std::left;
    std::cerr << "Usage: " << m_usage << "\n";

    for (auto &e : m_argval)
        std::cerr << "  " << std::setw(10) << e.first << " <" << e.second << ">\n"
                  << "    " << m_argdesc.at(e.first) << "\n";

    for (auto &e : m_flagval)
        std::cerr << "  " << e.first << "\n" << "    " << m_flagdesc.at(e.first) << "\n";
}

void CmdLine::parse(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        if (m_argval.count(argv[i]) != 0) {
            if (++i != argc && m_argval.count(argv[i]) == 0 && m_flagval.count(argv[i]) == 0)
                m_argval[argv[i-1]] = argv[i];
            else
                std::cerr << "ERROR: missing an argument for: " << argv[--i] << "\n";
        }
        else if (m_flagval.count(argv[i]) != 0) {
            m_flagval[argv[i]] = true;
        }
        else {
            std::cerr << "ERROR: unrecognized command line argument: " << argv[i] << "\n";
        }
    }
}

bool CmdLine::has(const string &name) const
{
    return m_flagval.at(name);
}

string CmdLine::get(const string &name) const
{
    return m_argval.at(name);
}
