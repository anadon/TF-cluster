#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>


#include "auxillaryUtilities.hpp"
#include "graph.t.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"


using std::string;
using std::ifstream;
using std::ofstream;
using std::iostream;
using std::vector;
using std::ifstream;
using std::endl;
using std::cout;
using std::cerr;


std::vector<string> tokenizeString(char *input, const char *delimitors){
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

  if( argc > 1 ) return EINVAL;

/*
  if( argc > 2 ) return EINVAL;

  if(argc == 2) freopen(argv[1], "r", stdin);
*/

  return 0;
}


void loadFromFile(graph<struct geneData, double> *geneNetwork, 
                            string geneListFile, string expressionFile){
  
  vector<string> geneList;
  string tmp, headerGeneName;
  
  ifstream geneListReader, exprCoeffReader;
  
  geneListReader.open(geneListFile);
  
  if(!geneListReader.good()){
    cerr << "could not access \"" << geneListFile << "\"" << endl;
    exit(1);
  }
  
  while(geneListReader.good()){
    struct geneData newGene;
    geneListReader >> newGene.name;
    geneList.push_back(newGene.name);
    geneNetwork->addVertex(newGene);
  }
  geneListReader.close();
  
  exprCoeffReader.open(expressionFile);
  
  if(!exprCoeffReader.good()){
    cerr << "could not access \"" << expressionFile << "\"" << endl;
    exit(1);
  }
  
  while(exprCoeffReader.good()){
    string rightVertexName;
    struct geneData testValue;
    vertex<struct geneData, double> *leftVertex;
    
    exprCoeffReader >> testValue.name;
    
    const size_t loopSize = geneList.size();
    leftVertex = geneNetwork->getVertexForValue(testValue);
    
    for(size_t i = 0; i < loopSize; i++){
      double weight;
      
      exprCoeffReader >> weight;
      testValue.name = geneList[i];
      
      geneNetwork->addEdge(leftVertex, geneNetwork->getVertexForValue(testValue), weight);
    }
  }
  exprCoeffReader.close();

}


void pruneGraph(graph<struct geneData, double> *geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    quickMergeEdges(geneNetwork->vertexArray[i], geneNetwork->vertexArray[i]->numEdges);
  }

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    for(size_t j = geneNetwork->vertexArray[i]->numEdges-1; j >= keepTopN; j--)
      geneNetwork->removeEdge(geneNetwork->vertexArray[i]->edges[j]);
  }
}


void quickMergeEdges(vertex<struct geneData, double> *toPrune,
                                                    const size_t size){
  size_t numRising;
  size_t i;

  numRising = 0;

  for(i = 0; i < size-1; i++)
    if(toPrune->edges[i]->weight < toPrune->edges[i+1]->weight) numRising++;

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    edge<struct geneData, double> *tmp;
    for(i = 0; i < size/2; i++){
      tmp = toPrune->edges[i];
      toPrune->edges[i] = toPrune->edges[(size-1) - i];
      toPrune->edges[(size-1) - i] = tmp;
    }
  }

  vector<size_t> indiciesOfInterest;
  indiciesOfInterest.push_back(0);

  i=0;
  while(i < size-1)
    if(toPrune->edges[i]->weight < toPrune->edges[i+1]->weight)
      indiciesOfInterest.push_back(++i);
  indiciesOfInterest.push_back(size);

  while(indiciesOfInterest.size() > 2){
    vector<size_t> newIndiciesOfInterest;
    for(i = 0; i < indiciesOfInterest.size()-2; i+=2){
      mergeHelper(toPrune->edges, indiciesOfInterest[i],
                  indiciesOfInterest[i+1],
                  indiciesOfInterest[i+2]);
      newIndiciesOfInterest.push_back(indiciesOfInterest[i]);
    }
    indiciesOfInterest = newIndiciesOfInterest;
    indiciesOfInterest.push_back(size);
  }

}


void mergeHelper(edge<struct geneData, double> **toSort, const size_t leftIndex,
                        const size_t rightIndex, const size_t endIndex){
  size_t leftParser, rightParser, mergedParser;
  edge<struct geneData, double> **sortSpace;

  sortSpace = (edge<struct geneData, double>**) malloc(sizeof(*sortSpace) * (endIndex - leftIndex));

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser] > toSort[rightParser] ?
        toSort[leftParser++] : toSort[rightParser++];

  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));

}


double calculateSigma(const edge<struct geneData, double> ** links,
                                                const size_t numEdges){
  double mean;
  double meanOfSquareDifferences;

  mean = 0;
  for(size_t i = 0; i < numEdges; i++)
    mean += links[i]->weight;
  mean /= numEdges;

  meanOfSquareDifferences = 0;
  for(size_t i = 0; i < numEdges; i++){
    double difference = links[i]->weight - mean;
    meanOfSquareDifferences += (difference * difference);
  }
  meanOfSquareDifferences /= numEdges;

  return sqrt(meanOfSquareDifferences);
}


void convertCoeffToSigmaValue(graph<struct geneData, double> *geneNetwork, double sigma){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    geneNetwork->edgeArray[i]->weight /= sigma;
}


void removeLowEdges(graph<struct geneData, double> *geneNetwork, 
                                                  const double &cutOff){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    if(geneNetwork->edgeArray[i]->weight < cutOff)
      geneNetwork->removeEdge(geneNetwork->edgeArray[i]);
}


void removeWeakVerticies(graph<struct geneData, double> *geneNetwork){
  //first, we need to remove all nodes which do not have 3 available
  //links.  It's technically possible, and technically possible for new
  //ones to occur after a given removal so this is expensive...
  bool disconnectedVerticiesFound;
  do{
    disconnectedVerticiesFound = false;
    for(size_t i = 0; i < geneNetwork->numVertexes; i++){
      if(3 > geneNetwork->vertexArray[i]->numEdges){
        geneNetwork->removeVertex(geneNetwork->vertexArray[i]);
        disconnectedVerticiesFound = true;
      }else{
        double threeSigmaReq, twoSigmaReq, oneSigmaReq;
        threeSigmaReq = geneNetwork->vertexArray[i]->edges[0]->weight;
        twoSigmaReq   = geneNetwork->vertexArray[i]->edges[1]->weight;
        oneSigmaReq   = geneNetwork->vertexArray[i]->edges[2]->weight;
        if(3.0 > threeSigmaReq || 2 > twoSigmaReq || 1 > oneSigmaReq){
          geneNetwork->removeVertex(geneNetwork->vertexArray[i]);
          disconnectedVerticiesFound = true;
        }
      }
    }
  }while(disconnectedVerticiesFound);
}


void simpleError(const char *message){
  cerr << "Error: " << message << endl;
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
    simpleError("cannot open required file \"SCCM_pipe.cfg\"");
    exit(1);
  }
  
  while(getline(configFile, line)){
    if(0 == line.size()) continue;
    if('#' == line[0]) continue;
    if(string::npos != line.find("geneListFile")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.geneListFile = line.substr(startIndex, endIndex - startIndex);
      continue;
    }
    if(string::npos != line.find("expressionFile")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.expressionFile = line.substr(startIndex, endIndex - startIndex);
      continue;
    }
    if(string::npos != line.find("topPick")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.topPick = atoi(line.substr(startIndex, endIndex - startIndex).c_str());
      continue;
    }
    if(string::npos != line.find("kickSize")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.kickSize = atoi(line.substr(startIndex, endIndex - startIndex).c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink1")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.tripleLink1 = atof(line.substr(startIndex, endIndex - startIndex).c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink2")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.tripleLink2 = atof(line.substr(startIndex, endIndex - startIndex).c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink3")){
      size_t startIndex = line.find("=");
      size_t endIndex = line.find("#");
      settings.tripleLink3 = atof(line.substr(startIndex, endIndex - startIndex).c_str());
      continue;
    }
  }
  configFile.close();
  
  if(settings.geneListFile != ""){
    cout << "geneList is in \"" << settings.geneListFile << "\"" << endl;
  }else{
    simpleError("geneList file not set");
    errorInFile = true;
  }
  if(settings.expressionFile != ""){
    cout << "expression is in \"" << settings.expressionFile << "\"" << endl;
  }else{
    simpleError("expression file not set");
    errorInFile = true;
  }
  
  cout << "topPick: " << settings.topPick << endl;
  
  cout << "kickSize: " << settings.kickSize << endl;
  
  if(settings.tripleLink1 < 0){
    simpleError("tripleLink1 is less than 0");
    errorInFile = true;
  }else if(settings.tripleLink1 < settings.tripleLink2){
    simpleError("tripleLink1 is smaller than tripleLink2");
    errorInFile = true;
  }else{
    cout << "tripleLink1 is " << settings.tripleLink1 << endl;
  }
  if(settings.tripleLink2 < 0){
    simpleError("tripleLink2 is less than 0");
    errorInFile = true;
  }else if(settings.tripleLink2 < settings.tripleLink3){
    simpleError("tripleLink2 is smaller than tripleLink3");
    errorInFile = true;
  }else{
    cout << "tripleLink2 is " << settings.tripleLink2 << endl;
  }
  if(settings.tripleLink3 < 0){
    simpleError("tripleLink3 is smaller than 0");
    errorInFile = true;
  }else{
    cout << "tripleLink3 is " << settings.tripleLink3 << endl;
  }
  
  if(errorInFile) exit(1);

  
  return settings;
}
