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
using std::cout;


vector<string> tokenizeString(string toTokenize, const char *delimitors){
  vector<string> tr;
  char *tokenStart, *inputCopy;
  size_t begin, end, diffLen;
  
  begin = 0;
  
  while(string::npos != begin){
    end = toTokenize.find_first_of(delimitors, begin);
    if(string::npos != end) diffLen = end - begin;
    else diffLen = end;
    tr.push_back(toTokenize.substr(begin, diffLen);
    begin = end;
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


void loadFromFile(graph<geneData, double> *geneNetwork, 
                            string geneListFile, string expressionFile){
  
  vector<string> geneList;
  string tmp, headerGeneName;
  
  ifstream geneListReader, exprCoeffReader;
  
  geneListReader.open(geneListFile);
  
  if(!geneListReader.good()){
    cout << "could not access \"" << geneListFile << "\"" << endl;
    exit(1);
  }
  
  while(geneListReader.good()){
    string geneName;
    geneListReader >> geneName;
    geneData newGene(geneName);
    if(newGene.name == ""){
      cout << "WARNING: empty gene name!" << endl;
      continue;
    }
    geneList.push_back(newGene.name);
    cout << "Registered gene \"" << newGene.name << "\"" << endl;
    geneNetwork->addVertex(newGene);
  }
  geneListReader.close();
  
  exprCoeffReader.open(expressionFile);
  
  if(!exprCoeffReader.good()){
    cout << "could not access \"" << expressionFile << "\"" << endl;
    exit(1);
  }
  
  while(exprCoeffReader.good()){
    string rightVertexName;
    vertex<geneData, double> *leftVertex, *rightVertex;
    vector<string> lineTokens;
    string line;
    
    getline(exprCoeffReader, line);
    lineTokens = tokenizeString(line, " \t,");
    
    geneData testValue(lineTokens[0]);
    
    leftVertex = geneNetwork->getVertexForValue(testValue);
    if(NULL == leftVertex){
      cout << "Cannot find gene \"" << testValue.name 
           << "\" in list of genes; Skipping row" << endl;
      continue;
    }
    
    for(size_t i = 1; i < lineTokens.size(); i++){
      double weight;
      string weightString;
      
      weightString = lineTokens[i];
      weight = atof(weightString.c_str());
      if(0 == weight){
        cout << "Likely error converting \"" << weightString 
             << "\" to double" << endl;
      }
      
      testValue.name = geneList[i];
      rightVertex = geneNetwork->getVertexForValue(testValue);
      
      if(NULL == rightVertex){
        cout << "Error looking up known good vertex!" << endl;
        exit(1);
      }
      geneNetwork->addEdge(leftVertex, rightVertex, weight);
    }
  }
  exprCoeffReader.close();

}


void pruneGraph(graph<geneData, double> *geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    quickMergeEdges(geneNetwork->vertexArray[i], geneNetwork->vertexArray[i]->numEdges);
  }

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    for(size_t j = geneNetwork->vertexArray[i]->numEdges-1; j >= keepTopN; j--)
      geneNetwork->removeEdge(geneNetwork->vertexArray[i]->edges[j]);
  }
}


void quickMergeEdges(vertex<geneData, double> *toPrune,
                                                    const size_t size){
  size_t numRising;
  size_t i;

  numRising = 0;

  for(i = 0; i < size-1; i++){
    double prev, next;
    prev = toPrune->edges[i]->weight;
    next = toPrune->edges[i+1]->weight;
    if(prev < next) numRising++;
  }

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    edge<geneData, double> *tmp;
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


void mergeHelper(edge<geneData, double> **toSort, const size_t leftIndex,
                        const size_t rightIndex, const size_t endIndex){
  size_t leftParser, rightParser, mergedParser;
  edge<geneData, double> **sortSpace;

  sortSpace = (edge<geneData, double>**) malloc(sizeof(*sortSpace) * (endIndex - leftIndex));

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser] > toSort[rightParser] ?
        toSort[leftParser++] : toSort[rightParser++];

  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));

}


double calculateSigma(const edge<geneData, double> ** links,
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


void convertCoeffToSigmaValue(graph<geneData, double> *geneNetwork, double sigma){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    geneNetwork->edgeArray[i]->weight /= sigma;
}


void removeLowEdges(graph<geneData, double> *geneNetwork, 
                                                  const double &cutOff){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    if(geneNetwork->edgeArray[i]->weight < cutOff)
      geneNetwork->removeEdge(geneNetwork->edgeArray[i]);
}


void removeWeakVerticies(graph<geneData, double> *geneNetwork){
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
  cout << "Error: " << message << endl;
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
    if(string::npos != line.find("geneList")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting geneList to " << readValue << endl;
      settings.geneListFile = readValue;
      continue;
    }
    if(string::npos != line.find("expression")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting expression to " << readValue << endl;
      settings.expressionFile = readValue;
      continue;
    }
    if(string::npos != line.find("topPick")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting topPick to " << readValue << endl;
      settings.topPick = atoi(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("kickSize")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting kickSize to " << readValue << endl;
      settings.kickSize = atoi(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink1")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting tripleLink1 to " << readValue << endl;
      settings.tripleLink1 = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink2")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#") + 1;
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting tripleLink2 to " << readValue << endl;
      settings.tripleLink2 = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink3")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cout << "setting tripleLink3 to " << readValue << endl;
      settings.tripleLink3 = atof(readValue.c_str());
      continue;
    }
    cout << "Skipped line \"" << line << "\"" << endl;
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
  
  cout << "topPick is " << (short) settings.topPick << endl;
  
  cout << "kickSize is " << settings.kickSize << endl;
  
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


void printClusters(vector< graph<geneData, double>* > clusters){
  for(size_t i = 0; i < clusters.size(); i++){
    cout << "Cluster " << (i+1) << ": " << endl;
    for(size_t j = 0; j < clusters[i]->numVertexes; j++){
      cout << clusters[i]->getVertexes()[j]->value.name << endl;
    }
    cout << endl;
  }
}


void printEdges(graph<geneData, double> *corrData){
  string left, right;
  double weight;
  
  cout << "Here are all the edges in the graph:" << endl;
  
  for(size_t i = 0; i < corrData->numEdges; i++){
    left = corrData->getEdges()[i]->left->value.name;
    right = corrData->getEdges()[i]->right->value.name;
    weight = corrData->getEdges()[i]->weight;
    cout << left << " <--" << weight << "--> " << right << endl;
  }
  
  cout << endl;
}