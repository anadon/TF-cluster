/*TODO liscense/programmer here*/

////////////////////////////////////////////////////////////////////////
//INCLUDES
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <queue>
#include <string>
//#include <thread>
#include <utility>


#include "auxillaryUtilities.hpp"
#include "graph.t.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"


////////////////////////////////////////////////////////////////////////
//USING IN NAMESPACE STD
////////////////////////////////////////////////////////////////////////

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::iostream;
using std::ofstream;
using std::pair;
using std::queue;
using std::string;
//using std::thread;
using std::vector;
using std::make_pair;


////////////////////////////////////////////////////////////////////////
//PRIVATE STRUCTS
////////////////////////////////////////////////////////////////////////

struct quickMergeDoubleSizeTPair_recurse_struct{
  pair<f64, size_t> *toSort;
  pair<f64, size_t> *workSpace;
  size_t size;
};


////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DECLARATIONS
////////////////////////////////////////////////////////////////////////

void *addTopEdgesHelper(void *protoArgs);

void sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort, 
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex, 
                                          pair<f64, size_t> *sortSpace);

void sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort, 
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex, 
                                          pair<f64, size_t> *sortSpace);


////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

//TODO: this can be made more complete
int verifyInput(int argc, char **argv, const string geneListFile){

  if( argc != 2 ) return EINVAL;

  int errorCode = 0;
  vector<string> geneList;
  vector< vector<string> > corrMatrix;
  string tmp, headerGeneName, line;
  
  ifstream geneListReader, exprCoeffReader;
  
  geneListReader.open(geneListFile);
  
  if(!geneListReader.good()){
    cout << "could not access \"" << geneListFile << "\"" << endl;
    return EIO;
  }
  
  while(geneListReader.good()){
    string geneName;
    geneListReader >> geneName;
    if(string("") == geneName)  continue;
    geneList.push_back(geneName);
  }
  geneListReader.close();

  return errorCode;
}


void quickMergeEdges(vertex<geneData, f64> *toPrune, csize_t size){
  size_t numRising;
  size_t i;

  numRising = 0;

  for(i = 0; i < size-1; i++){
    f64 prev, next;
    prev = toPrune->getEdges()[i]->weight;
    next = toPrune->getEdges()[i+1]->weight;
    if(prev < next) numRising++;
  }

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    edge<geneData, f64> *tmp;
    for(i = 0; i < size/2; i++){
      tmp = toPrune->getEdges()[i];
      toPrune->getEdges()[i] = toPrune->getEdges()[(size-1) - i];
      toPrune->getEdges()[(size-1) - i] = tmp;
    }
  }

  vector<size_t> indiciesOfInterest;
  indiciesOfInterest.push_back(0);

  i=0;
  for(i = 0; i < size-1; i++){
    if(toPrune->getEdges()[i]->weight < toPrune->getEdges()[i+1]->weight)
      indiciesOfInterest.push_back(i+1);
  }
  indiciesOfInterest.push_back(size);

  while(indiciesOfInterest.size() > 2){
    vector<size_t> newIndiciesOfInterest;
    for(i = 0; i < indiciesOfInterest.size()-2; i+=2){
      mergeHelper(toPrune->getEdges(), indiciesOfInterest[i],
                  indiciesOfInterest[i+1],
                  indiciesOfInterest[i+2]);
      newIndiciesOfInterest.push_back(indiciesOfInterest[i]);
    }
    indiciesOfInterest = newIndiciesOfInterest;
    indiciesOfInterest.push_back(size);
  }
}


void mergeHelper(edge<geneData, f64> **toSort, csize_t leftIndex,
                        csize_t rightIndex, csize_t endIndex){
  size_t leftParser, rightParser, mergedParser;
  edge<geneData, f64> **sortSpace;

  sortSpace = (edge<geneData, f64>**) malloc(sizeof(*sortSpace) * (endIndex - leftIndex));

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser] > toSort[rightParser] ?
        toSort[leftParser++] : toSort[rightParser++];

  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));
  
  free(sortSpace);
}


struct config loadConfig(const char *configFilePath){
  struct config settings;
  ifstream configFile;
  string line;
  bool errorInFile = false;
  
  settings.keepTopN = 100;
  settings.geneListFile = "";
  settings.expressionFile = "";
  settings.kickSize = 0;
  settings.threeSigma = -1;
  settings.twoSigma = -1;
  settings.oneSigma = -1;
  
  configFile.open(configFilePath);
  
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
      int readInt = atoi(readValue.c_str());
      if(0 > readInt || 255 < readInt){
        cerr << "Keep Top N value out of bounds (0,255]" << endl;
        exit(1);
      }
      settings.keepTopN = (u8) readInt;
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
      cerr << "setting threeSigma to " << readValue << endl;
      settings.threeSigma = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink2")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#") + 1;
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting twoSigma to " << readValue << endl;
      settings.twoSigma = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink3")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting oneSigma to " << readValue << endl;
      settings.oneSigma = atof(readValue.c_str());
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
  
  cerr << "topPick is " << (short) settings.keepTopN << endl;
  
  cerr << "kickSize is " << settings.kickSize << endl;
  
  if(settings.threeSigma < 0){
    cerr << "ERROR: threeSigma is less than 0" << endl;
    errorInFile = true;
  }else if(settings.threeSigma < settings.twoSigma){
    cerr << "ERROR: threeSigma is smaller than twoSigma" << endl;
    errorInFile = true;
  }else{
    cerr << "threeSigma is " << settings.threeSigma << endl;
  }
  if(settings.twoSigma < 0){
    cerr << "ERROR: twoSigma is less than 0" << endl;
    errorInFile = true;
  }else if(settings.twoSigma < settings.oneSigma){
    cerr << "ERROR: twoSigma is smaller than oneSigma" << endl;
    errorInFile = true;
  }else{
    cerr << "twoSigma is " << settings.twoSigma << endl;
  }
  if(settings.oneSigma < 0){
    cerr << "ERROR: oneSigma is smaller than 0" << endl;
    errorInFile = true;
  }else{
    cerr << "oneSigma is " << settings.oneSigma << endl;
  }
  
  if(errorInFile) exit(1);

  
  return settings;
}


void printClusters(queue< queue<size_t> > clusters, const vector<string> &names){
  for(size_t i = 0; !clusters.empty(); i++){
    cout << "cluster: " << (i+1) << endl;
    while(!clusters.front().empty()){
      cout << names[clusters.front().front()] << endl;
      clusters.front().pop();
    }
    clusters.pop();
  }
}


graph<geneData, f64>* constructGraph(const struct UDCorrelationMatrix &protoGraph, 
                                      cf64 threeSigma, cf64 oneSigma, cu8 maxNumEdges){
  graph<geneData, f64>* tr;
  size_t intermediateGraphItr;
  void *tmpPtr;
  pair<f64, size_t> *toSort, **intermediateGraph;
  size_t *columnSize;
  csize_t width = protoGraph.labels.size();
  
  cerr << "allocating preliminary memory..." << endl;
  intermediateGraphItr = 0;
  tmpPtr = malloc(sizeof(*intermediateGraph) * protoGraph.labels.size());
  intermediateGraph = (pair<f64, size_t>**) tmpPtr;
  tmpPtr = malloc(sizeof(*columnSize) * protoGraph.labels.size());
  columnSize = (size_t*) tmpPtr;
  
  //Construct an intermediate matrix to be entered into the graph.  This
  //is to avoid the high memory cost (and thus time cost) of the graph
  //over manipulating simpler data.
  for(size_t i = 0; i < protoGraph.labels.size(); i++){
    //cerr << "constructing intermediary matrix data entry..." << endl;
    size_t itr = 0;
    size_t baseOffset;
    tmpPtr = malloc(sizeof(*toSort) * (protoGraph.labels.size()-1));
    toSort = (pair<f64, size_t>*) tmpPtr;
    
    
    for(size_t j = 0; j < i; j++){
      baseOffset = (width*(width-1)/2) - (width-j)*((width-j)-1)/2 - j - 1;
      if(oneSigma <= protoGraph.UDMatrix[baseOffset + i])
        toSort[itr++] = pair<f64, size_t>(protoGraph.UDMatrix[baseOffset + i], j);
    }
    
    baseOffset = (width*(width-1)/2) - (width-i)*((width-i)-1)/2 - i - 1;
    for(size_t j = i+1; j < protoGraph.labels.size(); j++)
      if(oneSigma <= protoGraph.UDMatrix[baseOffset + j])
        toSort[itr++] = pair<f64, size_t>(protoGraph.UDMatrix[baseOffset + j], j);
    
    tmpPtr = realloc(toSort, sizeof(*toSort) * itr);
    toSort = (pair<f64, size_t>*) tmpPtr;
    
    if(0 == itr){
      columnSize[intermediateGraphItr] = 0;
      intermediateGraph[intermediateGraphItr++] = NULL;
      continue;
    }
    
    sortDoubleSizeTPairHighToLow(toSort, itr);
    
    size_t shrinkSize;
    if(toSort[0].first < threeSigma)  shrinkSize = 0;//Conditional jump or move depends on uninitialised value(s)
    else if(itr < maxNumEdges)  shrinkSize = itr;
    else                        shrinkSize = maxNumEdges;
    for(; shrinkSize != 0 && toSort[shrinkSize-1].first < oneSigma; shrinkSize--);
    columnSize[intermediateGraphItr] = shrinkSize;
    
    if(0 == shrinkSize){
      free(toSort);
      intermediateGraph[intermediateGraphItr++] = NULL;
    }else{
      tmpPtr = realloc(toSort, sizeof(*toSort) * shrinkSize);
      intermediateGraph[intermediateGraphItr++] = (pair<f64, size_t>*) tmpPtr;
    }
  }
  
  //Now that we don't need the very large UDMatrix in protoGraph, free
  //it.  This is to try very hard to stay within the memory envelope of
  //the original program.  Granted, I'm not sure the original program 
  //used double floating point values.
  free(protoGraph.UDMatrix);

  //now prepare the graph for all the data it is about to recieve, else
  //after the fact memory allocations can take minutes.
  cerr << "allocating for final graph..." << endl;
  tr = new graph<geneData, f64>();
  size_t totalVertexes, totalEdges;
  totalVertexes = totalEdges = 0;
  
  for(size_t i = 0; i < protoGraph.labels.size(); i++)
    if(columnSize[i]){
      totalVertexes++;
      totalEdges += columnSize[i];
    }
  
  tr->hintNumVertexes(totalVertexes);
  tr->hintNumEdges(totalEdges);
  
  for(size_t i = 0; i < protoGraph.labels.size(); i++){
    if(columnSize[i]){
      vertex<geneData, double> *target = tr->addVertex(geneData(i));
      target->hintNumEdges(columnSize[i]);
    }
  }
  
  for(size_t i = 0; i < protoGraph.labels.size(); i++){
    if(columnSize[i]){
      vertex<geneData, double> *leftVert;
      leftVert = tr->getVertexForValue(geneData(i));
      for(size_t j = 0; j < columnSize[i]; j++){
        vertex<geneData, double> *rightVert;
        rightVert = tr->getVertexForValue(geneData(intermediateGraph[i][j].second));
        if(NULL != rightVert){
          tr->addEdge(leftVert, rightVert, intermediateGraph[i][j].first);
        }
      }
    }
  }
  
  cerr << "minimizing memory in use by the graph..." << endl;
  for(size_t i = 0; i < protoGraph.labels.size(); i++){
    if(columnSize[i]){
      vertex<geneData, double> *target;
      target = tr->getVertexForValue(geneData(i));
      target->shrinkToFit();
    }
  }
  
  return tr;
}


void pruneGraph(graph<geneData, f64> *corrData, u8 keepTopN){
  edge<geneData, f64> **workSpace, **sortSpace;
  workSpace = (edge<geneData, f64>**) malloc(sizeof(*workSpace) * corrData->getNumEdges());
  sortSpace = (edge<geneData, f64>**) malloc(sizeof(*sortSpace) * corrData->getNumEdges());
  
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    queue<size_t> indiciesOfInterest;
    queue<size_t> IOISwap;
    edge<geneData, f64> *swap;
    size_t count[2] = {0, 0};
    
    vertex<geneData, f64> *target = corrData->getVertexes()[i];
    size_t size = target->getNumEdges();
    
    if(size <= keepTopN) continue;
    
    for(size_t j = 0; j < size; j++)
      workSpace[j] = target->getEdges()[j];
    
    //Quickmerge start//////////////////////////////////////////////////
    
    
    for(size_t j = 0; j < size-1; j++){
      if(workSpace[j]->weight < workSpace[j+1]->weight) count[0]++;
      else if(workSpace[j]->weight > workSpace[j+1]->weight) count[1]++;
    }
    
    if(count[0] > count[1])
      for(size_t j = 0; j <= size/2; j++){
        swap = workSpace[j];
        workSpace[j] = workSpace[(size-1) - j];
        workSpace[(size - 1) - j] = swap;
      }
    
    indiciesOfInterest.push(0);
    for(size_t j = 1; j < size; j++)
      if(workSpace[j-1]->weight < workSpace[j]->weight)
        indiciesOfInterest.push(j);
    indiciesOfInterest.push(size);
    
    
    while(indiciesOfInterest.size() > 2){
      while(indiciesOfInterest.size() > 2){
        size_t leftPtr, rightPtr, itr;
        csize_t leftStart = indiciesOfInterest.front(); indiciesOfInterest.pop();
        csize_t leftEnd = indiciesOfInterest.front();  indiciesOfInterest.pop();
        csize_t rightEnd = indiciesOfInterest.front();
        itr = 0;
        leftPtr = leftStart;
        rightPtr = leftEnd;
        csize_t copySize = (rightEnd - leftStart) * sizeof(*sortSpace);
        IOISwap.push(leftStart);
        
        while(leftPtr < leftEnd && rightPtr < rightEnd){
          if(workSpace[leftPtr]->weight > workSpace[rightPtr]->weight)
            sortSpace[itr++] = workSpace[leftPtr++];
          else
            sortSpace[itr++] = workSpace[rightPtr++];
        }
        while(leftPtr < leftEnd)   sortSpace[itr++] = workSpace[leftPtr++];
        while(rightPtr < rightEnd) sortSpace[itr++] = workSpace[rightPtr++];
        
        memcpy(&workSpace[leftStart], sortSpace, copySize);
      }
      if(indiciesOfInterest.size() == 2){
        IOISwap.push(indiciesOfInterest.front());
        indiciesOfInterest.pop();
      }
      IOISwap.push(size);
      
      indiciesOfInterest = IOISwap;
      while(!IOISwap.empty()) IOISwap.pop();
    }
    
    //Quickmerge end////////////////////////////////////////////////////
    
    for(size_t i = keepTopN; i < size; i++)
      corrData->removeEdge(workSpace[i]);
    
  }
  free(workSpace);
  free(sortSpace);
}


//quickMerge
void sortDoubleSizeTPairHighToLow(pair<f64, size_t> *toSort, csize_t size){
  size_t numRising;
  size_t i;
  void *tmpPtr;
  
  if(1 >= size) return;

  numRising = 0;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first < toSort[i+1].first) numRising++;
  }

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    pair<f64, size_t> tmp;
    for(i = 0; i < size/2; i++){
      tmp = toSort[i];
      toSort[i] = toSort[(size-1) - i];
      toSort[(size-1) - i] = tmp;
    }
  }

  size_t* indiciesOfInterest = new size_t[size];
  indiciesOfInterest[0] = 0;
  size_t IOISize = 1;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first < toSort[i+1].first)
      indiciesOfInterest[IOISize++] = (i+1);
  }
  indiciesOfInterest[IOISize++] = size;

  pair<f64, size_t> *sortSpace;
  tmpPtr = malloc(sizeof(*sortSpace) * size);
  sortSpace = (pair<f64, size_t>*) tmpPtr;
  
  size_t *newIndiciesOfInterest = new size_t[size];
  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){
      
      sortDoubleSizeTPairHighToLowHelper(toSort, indiciesOfInterest[i], 
                      indiciesOfInterest[i+1], indiciesOfInterest[i+2], 
                                                            sortSpace);

      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[i];
    }
    if(!(IOISize & 1)){
      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[IOISize-2];
    }
    newIndiciesOfInterest[NIOISize++] = size;
    memcpy(indiciesOfInterest, newIndiciesOfInterest, NIOISize * sizeof(*indiciesOfInterest));
    IOISize = NIOISize;
  }
  
  delete indiciesOfInterest;
  delete newIndiciesOfInterest;
  free(sortSpace);
  
}


void sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort, 
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex, 
                                          pair<f64, size_t> *sortSpace){
  size_t leftParser, rightParser, mergedParser;

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser].first > toSort[rightParser].first ?
        toSort[leftParser++] : toSort[rightParser++];

  while(leftParser < rightIndex)
    sortSpace[mergedParser++] = toSort[leftParser++];
  while(rightParser < endIndex)
    sortSpace[mergedParser++] = toSort[rightParser++];
  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));
}


void sortDoubleSizeTPairLowToHigh(pair<f64, size_t> *toSort, csize_t size){
  size_t numFalling;
  size_t i;
  void *tmpPtr;
  
  if(1 >= size) return;

  numFalling = 0;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first > toSort[i+1].first) numFalling++;
  }

  if(numFalling > (size >> 1)){
    //reverse so that more are in order
    pair<f64, size_t> tmp;
    for(i = 0; i < size/2; i++){
      tmp = toSort[i];
      toSort[i] = toSort[(size-1) - i];
      toSort[(size-1) - i] = tmp;
    }
  }

  size_t* indiciesOfInterest = new size_t[size];
  indiciesOfInterest[0] = 0;
  size_t IOISize = 1;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first > toSort[i+1].first)
      indiciesOfInterest[IOISize++] = (i+1);
  }
  indiciesOfInterest[IOISize++] = size;

  pair<f64, size_t> *sortSpace;
  tmpPtr = malloc(sizeof(*sortSpace) * size);
  sortSpace = (pair<f64, size_t>*) tmpPtr;
  
  size_t *newIndiciesOfInterest = new size_t[size];
  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){
      
      sortDoubleSizeTPairLowToHighHelper(toSort, indiciesOfInterest[i], 
                      indiciesOfInterest[i+1], indiciesOfInterest[i+2], 
                                                            sortSpace);

      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[i];
    }
    if(!(IOISize & 1)){
      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[IOISize-2];
    }
    newIndiciesOfInterest[NIOISize++] = size;
    memcpy(indiciesOfInterest, newIndiciesOfInterest, NIOISize * sizeof(*indiciesOfInterest));
    IOISize = NIOISize;
  }
  
  delete indiciesOfInterest;
  delete newIndiciesOfInterest;
  free(sortSpace);
  
}


void sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort, 
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex, 
                                          pair<f64, size_t> *sortSpace){
  size_t leftParser, rightParser, mergedParser;

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser].first < toSort[rightParser].first ?
        toSort[leftParser++] : toSort[rightParser++];

  while(leftParser < rightIndex)
    sortSpace[mergedParser++] = toSort[leftParser++];
  while(rightParser < endIndex)
    sortSpace[mergedParser++] = toSort[rightParser++];
  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));
}