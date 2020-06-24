#include "cmdline.h"
#include "print.h"
#include "stringutil.h"

void CmdLine::add(const std::string &arg)
{
    flag_[arg] = false;
}

void CmdLine::add(const std::string &arg, const std::string &val)
{
    arg_[arg] = val;
}

bool CmdLine::has(const std::string &arg) const
{
    return flag_.at(arg);
}

std::string CmdLine::get(const std::string &arg) const
{
    return arg_.at(arg);
}

void CmdLine::parse(int argc, char *argv[])
{
    std::string key;
    std::vector<std::string> val;

    for (int i = 1; i < argc; ++i) {
        if ( starts_with(argv[i], "--") ) {
            if ( ! key.empty() ) {
                take(key, join(val,","));
                key.clear();
                val.clear();
            }
            key = argv[i];
            val.clear();
        }
        else {
            if ( ! key.empty() )
                val.push_back(argv[i]);
            else
                eprint("ERROR: unrecognized command line argument: %s\n", argv[i]);
        }
    }

    if ( ! key.empty() )
        take(key, join(val,","));
}

void CmdLine::take(const std::string &key, const std::string &val)
{
    if (arg_.count(key) != 0) {
        if (!val.empty())
            arg_[key] = val;
        else
            eprint("ERROR: missing an argument for: %s\n", key);
    }
    else if (flag_.count(key) != 0) {
        if (val.empty())
            flag_[key] = true;
        else
            eprint("ERROR: argument is not allowed for: %s\n", key);
    }
    else
        eprint("ERROR: unrecognized command line: %s\n", key);
}
