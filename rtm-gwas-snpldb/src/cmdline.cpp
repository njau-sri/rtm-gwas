#include "cmdline.h"

#include <set>

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

int CmdLine::parse(int argc, char *argv[])
{
    std::string key;
    std::vector<std::string> val;
    std::set<std::string> dupchk;

    for (int i = 1; i < argc; ++i) {
        if (starts_with(argv[i], "--")) {
            if (!key.empty()) {
                if (take(key, join(val, ",")) != 0)
                    return 1;
            }

            key = argv[i];
            val.clear();

            if (!dupchk.insert(key).second) {
                eprint("ERROR: duplicate argument is not allowed: %s\n", key);
                return 1;
            }
        }
        else {
            if (key.empty()) {
                eprint("ERROR: unrecognized command line argument: %s\n", argv[i]);
                return 1;
            }
            val.push_back(argv[i]);
        }
    }

    if (!key.empty())
        return take(key, join(val, ","));

    return 0;
}

int CmdLine::take(const std::string &key, const std::string &val)
{
    if (arg_.count(key) != 0) {
        if (val.empty()) {
            eprint("ERROR: missing an argument for: %s\n", key);
            return 1;
        }
        arg_[key] = val;
        return 0;
    }

    if (flag_.count(key) != 0) {
        if (!val.empty()) {
            eprint("ERROR: argument is not allowed for: %s\n", key);
            return 1;
        }
        flag_[key] = true;
        return 0;
    }

    eprint("ERROR: unrecognized command line argument: %s\n", key);

    return 1;
}
