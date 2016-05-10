#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <errno.h>
#include <cblas.h>
#include <vector>
#include <queue>
#include <utility>
#include <string.h>

using namespace std;
/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project 
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/



vector<string> tokenizeString(char *input, 
                                      const char *delimitors = " \t,"){
  vector<string> tr;
  char *tokenStart;
  
  
  tokenStart = strtok(input, delimitors);
  
  while(tokenStart){
    tr.push_back(string(tokenStart));
    tokenStart = strtok(NULL, delimitors);
  }
  
  return tr;
}


class pairTable{
  public:
  vector<vector<pair<string, double>>> data;
  vector<string> columnLabels;
  
  void loadFromFile(FILE *input){
    char *line = NULL;
    vector<string> tokens;

    getline(&line, NULL, stdin);
    tokens = tokenizeString(line);
    
    while(feof(stdin)){
      free(line); line = NULL;
      for(size_t i = 0; i < tokens.size(); i+=2){
        data[i].push_back(pair<string, double>(tokens[i], atof(tokens[i+1].c_str())));
      }
    }
    
    free(line); line = NULL;
    
    sortColumns();
  }
  
  
  void sortColumns(){
    for(size_t i = 0; i < data.size(); i++) data[i] = sortColumn(data[i]);
  }
  
  //quickmerge, high to low
  vector<pair<string, double>> sortColumn(vector<pair<string, double>> toSort){
    if(toSort.size() < 2) return toSort;
    size_t i = 0, j;
    while(toSort[i].second < toSort[i+1].second);
    
    vector<pair<string, double>> left, right, merged;
    
    //TODO
    for(j = 0; j < i; j++) left.push_back(toSort[i]);
    for(;j < toSort.size(); j++) right.push_back(toSort[i]);
    
    left = sortColumn(left);
    right = sortColumn(right);
    
    i = j = 0;
    while(i < left.size() && j < right.size()){
      if(left[i] > right[j]){
        merged.push_back(left[i++]);
      }else{
        merged.push_back(right[j++]);
      }
    }
    
    while(i < left.size()) merged.push_back(left[i++]);
    while(j < right.size()) merged.push_back(right[j++]);
    
    return merged;
  }
  
  
};


pairTable prune(const pairTable toPrune, unsigned char keepTopN = 100){
  pairTable pruned;
  
  for(size_t i = 0; i < toPrune.columnLabels.size(); i++){
    pruned.columnLabels.push_back(toPrune.columnLabels[i]);
    pruned.data.push_back(vector<pair<string, double>>());
    for(unsigned char j = 0; j < keepTopN; j++){
      pruned.data[i].push_back(toPrune.data[i][j]);
    }
  }
  
  return pruned;
}


int verifyInput(int argc, char **argv){
  
  if( argc > 2 ) return EINVAL;
  
  if(argc == 2) freopen(argv[1], "r", stdin);
  
  return 0;
}


int main(int argc, char **argv){
  pairTable inputTable;
  
  int error = verifyInput(argc, argv);
  if(error) return error;
  
  pairTable corrData;
  
  corrData.loadFromFile(stdin);
  
  corrData = prune(corrData);
  
  return 0;
}