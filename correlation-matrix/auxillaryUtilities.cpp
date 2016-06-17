
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <queue>
#include <string>
#include <thread>
#include <utility>


#include "auxillaryUtilities.hpp"
#include "statistics.hpp"


using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::iostream;
using std::ofstream;
using std::pair;
using std::queue;
using std::string;
using std::thread;
using std::vector;


vector<string> tokenizeString(const string toTokenize, const char *delimitors){
  vector<string> tr;
  size_t begin, end;
  
  begin = 0;
  
  do{
    end = toTokenize.find_first_of(delimitors, begin);
    string pToken = toTokenize.substr(begin, end - begin);
    if(string("") != pToken)  tr.push_back(pToken);
    begin = end+1;
  }while(end+1);

  return tr;
}


int verifyInput(int argc, char **argv, string geneListFile, string expressionFile){

  if( argc > 1 ) return EINVAL;

  int errorCode = 0;
  vector<string> geneList;
  vector< vector<string> > corrMatrix;
  string tmp, headerGeneName, line;
  
  ifstream geneListReader, exprCoeffReader;
  
  geneListReader.open(geneListFile);
  
  if(!geneListReader.good()){
    cerr << "could not access \"" << geneListFile << "\"" << endl;
    return EIO;
  }
  
  while(geneListReader.good()){
    string geneName;
    geneListReader >> geneName;
    if(string("") == geneName)  continue;
    geneList.push_back(geneName);
  }
  geneListReader.close();
  
  
  exprCoeffReader.open(expressionFile);
  
  if(!exprCoeffReader.good()){
    cerr << "could not access \"" << expressionFile << "\"" << endl;
    return EIO;
  }
  
  while(exprCoeffReader.good()){
    getline(exprCoeffReader, line);
    corrMatrix.push_back(tokenizeString(line, " \t,\n"));
  }
  exprCoeffReader.close();
  
  for(size_t i = 0; i < corrMatrix.size(); i++){
    for(size_t j = 1; j < corrMatrix[i].size(); j++){
      double weight;
      weight = atof(corrMatrix[i][j].c_str());
      if(0 == weight){
        cerr << "Error: Line " << i << " Token " << j << 
                " has invalid value \"" << corrMatrix[i][j] << "\"" 
             << endl;
        errorCode = EINVAL;
      }
    
    }
  }


  return errorCode;
}


struct loadFromFileReturn loadFromFile(string geneListFile, 
                                                string expressionFile){
  
  vector<string> colHeaders, rowHeaders;
  string tmp, headerGeneName, line;
  ifstream geneListReader, exprCoeffReader;
  struct loadFromFileReturn tr = {NULL, 0, 0, vector<string>()};
  
  geneListReader.open(geneListFile);
  
  while(geneListReader.good()){
    string geneName;
    geneListReader >> geneName;
    if(string("") == geneName)  continue;
    colHeaders.push_back(geneName);
  }
  geneListReader.close();
  
  
  exprCoeffReader.open(expressionFile);
  
  while(exprCoeffReader.good()){
    getline(exprCoeffReader, line);
    if(string("") != line){
      vector<string> toParse = tokenizeString(line);
      tr.genes.push_back(toParse[0]);
      tr.numCols = toParse.size()-1;
      tr.matrix = (f64**) realloc(tr.matrix, (tr.numRows + 1) * sizeof(*tr.matrix));
      tr.matrix[tr.numRows] = (f64*) malloc(tr.numCols * sizeof(**tr.matrix));
      for(size_t i = 0; i < tr.numCols; i++)
        tr.matrix[tr.numRows][i] = atof(toParse[i+1].c_str());
      inplaceCenterMean(tr.matrix[tr.numRows], tr.numCols);
      tr.numRows++;
    }
  }
  
  exprCoeffReader.close();
  
  return tr;
}


struct correlationMatrix calculateCorrelationMatrix(const struct loadFromFileReturn &input){

  struct correlationMatrix tr;
  pthread_t *workers;
  struct corrHelpStruct *instructions;
  int *toIgnore;
  
  tr.numElements = (input.numRows * (input.numRows -1)) >> 1;
  
  tr.UDMatrix = (double*) malloc(sizeof(*tr.UDMatrix) * tr.numElements);
  
  //TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  cu32 numCPUs = std::thread::hardware_concurrency() < tr.numElements ? 
                 std::thread::hardware_concurrency() : tr.numElements ;
  
  workers = (pthread_t*) malloc(sizeof(*workers) * numCPUs);
  instructions = (struct corrHelpStruct*) malloc(sizeof(*instructions) * numCPUs);
  
  for(size_t i = 0; i < numCPUs; i++){
    instructions[i].side = input.numRows;
    instructions[i].start = (i * tr.numElements) / numCPUs;
    instructions[i].end = ((i+1) * tr.numElements) / numCPUs;
    instructions[i].matrix = input.matrix;
    instructions[i].vecLeng = input.numCols;
    instructions[i].results = tr.UDMatrix;
  }
  
  for(size_t i = 0; i < numCPUs; i++)
    pthread_create(&workers[i], NULL, correlationHelper, &instructions[i]);
  
  for(size_t i = 0; i < numCPUs; i++)
    pthread_join(workers[i], (void**) &toIgnore);
  
  free(workers);
  free(instructions);
  
  
  return tr;
}


void *correlationHelper(void *protoArgs){
  size_t i, j, itr;
  const struct corrHelpStruct *args = (struct corrHelpStruct*) protoArgs;
  f64 *results = args->results;
  cf64 **matrix = (cf64**) args->matrix;
  csize_t side = args->side;
  csize_t end = args->end;
  csize_t vecLeng = args->vecLeng;
  
  itr = args->start;
  
  //needed because of triangular addressing
  i = args->side - sqrtl((args->side * args->side) - 2*itr);//y
  j = itr - i*(args->side - i/2.0);//x
  
  if(i >= args->side){
    cerr << "i is too large!" << endl; exit(1);
  }
  if(j >= args ->side){
    cerr << "j is too large!" << endl; exit(1);
  }
  
  //primer loop
  for(;j < side && itr < end; j++)
    results[itr++] = getCenteredCorrelation(matrix[i], matrix[j], vecLeng);
  
  
  //main loop
  for(i++;itr < end; i++){
    for(j = i+1; j < side && itr < end; j++){
      results[itr++] = getCenteredCorrelation(matrix[i], matrix[j], vecLeng);
    }
  }
  
  return NULL;
}



void convertCoeffToSigmaValue(f64 *UDMatrix, csize_t size){
  
  cf64 sigma = calculateCenteredStandardDeviation(UDMatrix, size);
  
  for(size_t i = 0; i < size; i++)
    UDMatrix[i] /= sigma;
}



struct config loadConfig(){
  struct config settings;
  ifstream configFile;
  string line;
  bool errorInFile = false;
  
  settings.topPick = 100;
  settings.geneListFile = "";
  settings.expressionFile = "";
  settings.kickSize = 0;
  settings.tripleLink1 = -1;
  settings.tripleLink2 = -1;
  settings.tripleLink3 = -1;
  
  configFile.open("SCCM_pipe.cfg");
  
  if(!configFile.is_open()){
    cerr << "ERROR: cannot open required file \"SCCM_pipe.cfg\"" << endl;
    exit(1);
  }
  
  while(getline(configFile, line)){
    if(0 == line.size()) continue;
    if('#' == line[0]) continue;
    if(string::npos != line.find("geneList")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting geneList to " << readValue << endl;
      settings.geneListFile = readValue;
      continue;
    }
    if(string::npos != line.find("expression")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting expression to " << readValue << endl;
      settings.expressionFile = readValue;
      continue;
    }
    if(string::npos != line.find("topPick")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting topPick to " << readValue << endl;
      settings.topPick = atoi(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("kickSize")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting kickSize to " << readValue << endl;
      settings.kickSize = atoi(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink1")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting tripleLink1 to " << readValue << endl;
      settings.tripleLink1 = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink2")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#") + 1;
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting tripleLink2 to " << readValue << endl;
      settings.tripleLink2 = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink3")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting tripleLink3 to " << readValue << endl;
      settings.tripleLink3 = atof(readValue.c_str());
      continue;
    }
    cerr << "Skipped line \"" << line << "\"" << endl;
  }
  configFile.close();
  
  if(settings.geneListFile != ""){
    cerr << "geneList is in \"" << settings.geneListFile << "\"" << endl;
  }else{
    cerr << "ERROR: geneList file not set" << endl;
    errorInFile = true;
  }
  if(settings.expressionFile != ""){
    cerr << "expression is in \"" << settings.expressionFile << "\"" << endl;
  }else{
    cerr << "ERROR: expression file not set" << endl;
    errorInFile = true;
  }
  
  cerr << "topPick is " << (short) settings.topPick << endl;
  
  cerr << "kickSize is " << settings.kickSize << endl;
  
  if(settings.tripleLink1 < 0){
    cerr << "ERROR: tripleLink1 is less than 0" << endl;
    errorInFile = true;
  }else if(settings.tripleLink1 < settings.tripleLink2){
    cerr << "ERROR: tripleLink1 is smaller than tripleLink2" << endl;
    errorInFile = true;
  }else{
    cerr << "tripleLink1 is " << settings.tripleLink1 << endl;
  }
  if(settings.tripleLink2 < 0){
    cerr << "ERROR: tripleLink2 is less than 0" << endl;
    errorInFile = true;
  }else if(settings.tripleLink2 < settings.tripleLink3){
    cerr << "ERROR: tripleLink2 is smaller than tripleLink3" << endl;
    errorInFile = true;
  }else{
    cerr << "tripleLink2 is " << settings.tripleLink2 << endl;
  }
  if(settings.tripleLink3 < 0){
    cerr << "ERROR: tripleLink3 is smaller than 0" << endl;
    errorInFile = true;
  }else{
    cerr << "tripleLink3 is " << settings.tripleLink3 << endl;
  }
  
  if(errorInFile) exit(1);

  
  return settings;
}


void printEdges(const struct loadFromFileReturn &fileData, 
                              const struct correlationMatrix &corrData){
  size_t itr;
  
  itr = 0;
  for(size_t i = 0; i < fileData.numRows; i++)
    for(size_t j = i+1; j < fileData.numRows; j++)
      fprintf(stdout, "%s\t%s\t%fl\n", fileData.genes[i].c_str(), 
                  fileData.genes[j].c_str(), corrData.UDMatrix[itr++]);
  
  fflush(stdout);
}


