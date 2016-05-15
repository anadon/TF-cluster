#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;

typedef unsigned char u8;

vector<string> tokenizeString(char *input, 
                                      const char *delimitors = " \t,");

int verifyInput(int argc, char **argv);

#endif