#! /usr/bin/env python

from time import strftime, localtime

hermes_common_path = "../../hermes_common"

def subs(s, d):
    for k in d:
        s = s.replace(k, d[k])
    return s

def convert_cpp():
    lines = open(hermes_common_path + "/_hermes_common_api_new.h").readlines()
    line = lines[0]
    while not line.startswith("extern"):
        del lines[0]
        line = lines[0]
    repl_dict = {}
    while line.startswith("extern"):
        line = line.strip()
        repl_dict[line.replace("extern", "static")] = line
        del lines[0]
        line = lines[0]

    lines = open("_hermes2d.cpp").readlines()
    if len(lines) > 2 and lines[1].startswith("/* Corrected by convert_api.py"):
        return
    f = open("_hermes2d.cpp", "w")
    f.write(lines[0])
    del lines[0]
    f.write("/* Corrected by convert_api.py on %s */\n" % \
            strftime("%a %b %d %H:%M:%S %Y", localtime()))
    orig = "".join(lines)
    orig = subs(orig, repl_dict)
    f.write(orig)

if __name__ == "__main__":
    convert_cpp()
