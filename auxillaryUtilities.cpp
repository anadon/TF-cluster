#include <string.h>

#include "auxillaryUtilities.hpp"


vector<string> tokenizeString(char *input, const char *delimitors){
  vector<string> tr;
  char *tokenStart;
  
  
  tokenStart = strtok(input, delimitors);
  
  while(tokenStart){
    tr.push_back(string(tokenStart));
    tokenStart = strtok(NULL, delimitors);
  }
  
  return tr;
}


int verifyInput(int argc, char **argv){
  
  if( argc > 2 ) return EINVAL;
  
  if(argc == 2) freopen(argv[1], "r", stdin);
  
  return 0;
}