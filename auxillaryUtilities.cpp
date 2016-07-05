#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <queue>
#include <string>
#include <thread>
#include <utility>


#include "auxillaryUtilities.hpp"
#include "graph.t.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"


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
using std::make_pair;


int verifyInput(int argc, char **argv, const string geneListFile){

  if( argc > 1 ) return EINVAL;

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


struct loadFromFileReturn loadFromFile(string geneListFile){
  
  vector<string> colHeaders, rowHeaders;
  vector< vector<f64> > corrMatrix;
  string tmp, headerGeneName, line;
  ifstream geneListReader, exprCoeffReader;
  struct loadFromFileReturn tr;
  
  //threads alloc this
  //toIgnore = (int*) malloc(sizeof(*toIgnore));
  
  geneListReader.open(geneListFile);
  
  while(geneListReader.good()){
    string geneName;
    geneListReader >> geneName;
    if(string("") == geneName)  continue;
    colHeaders.push_back(geneName);
  }
  geneListReader.close();
  
  
  
  return tr;
}



/*
void pruneGraph(graph<geneData, f64> *geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++){
    quickMergeEdges(&geneNetwork->getVertexes()[i], geneNetwork->getVertexes()[i]->getNumEdges());
  }

  for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++){
    vertex<geneData, double> *target = geneNetwork->getVertexes()[i];
    for(size_t j = target->getNumEdges()-1; j >= keepTopN; j--)
      geneNetwork->removeEdge(target->getEdges()[j]);
  }
}
*/


void quickMergeEdges(vertex<geneData, f64> *toPrune,
                                                    csize_t size){
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
  while(i < size-1)
    if(toPrune->getEdges()[i]->weight < toPrune->getEdges()[i+1]->weight)
      indiciesOfInterest.push_back(++i);
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


void removeLowEdges(graph<geneData, f64> *geneNetwork, cf64 &cutOff){
  for(size_t i = 0; i < geneNetwork->getNumEdges(); i++){
    edge<geneData, f64> *target = geneNetwork->getEdges()[i];
    if(target->weight < cutOff)
      geneNetwork->removeEdge(target);
  }
}


void removeWeakVerticies(graph<geneData, f64> *geneNetwork, cf64 TL1, 
                                                    cf64 TL2, cf64 TL3){
  //first, we need to remove all nodes which do not have 3 available
  //links.  It's technically possible, and technically possible for new
  //ones to occur after a given removal so this is expensive...
  bool disconnectedVerticiesFound;
  do{
    disconnectedVerticiesFound = false;
    for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++){
      vertex<geneData, f64> *target = geneNetwork->getVertexes()[i];
      if(3 > target->getNumEdges()){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
      }else{
        cf64 threeSigmaReq = target->getEdges()[0]->weight;
        cf64 twoSigmaReq   = target->getEdges()[1]->weight;
        cf64 oneSigmaReq   = target->getEdges()[2]->weight;
        if(TL1 > threeSigmaReq || TL2 > twoSigmaReq || TL3 > oneSigmaReq){
          geneNetwork->removeVertex(target);
          disconnectedVerticiesFound = true;
        }
      }
    }
  }while(disconnectedVerticiesFound);
}


struct config loadConfig(){
  struct config settings;
  ifstream configFile;
  string line;
  bool errorInFile = false;
  
  settings.keepTopN = 100;
  settings.geneListFile = "";
  settings.expressionFile = "";
  settings.kickSize = 0;
  settings.tripleLink1 = -1;
  settings.tripleLink2 = -1;
  settings.tripleLink3 = -1;
  
  configFile.open("SCCM_pipe.cfg");
  
  if(!configFile.is_open()){
    cout << "ERROR: cannot open required file \"SCCM_pipe.cfg\"" << endl;
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
      settings.keepTopN = atoi(readValue.c_str());
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
    cout << "ERROR: geneList file not set" << endl;
    errorInFile = true;
  }
  if(settings.expressionFile != ""){
    cout << "expression is in \"" << settings.expressionFile << "\"" << endl;
  }else{
    cout << "ERROR: expression file not set" << endl;
    errorInFile = true;
  }
  
  cout << "topPick is " << (short) settings.keepTopN << endl;
  
  cout << "kickSize is " << settings.kickSize << endl;
  
  if(settings.tripleLink1 < 0){
    cout << "ERROR: tripleLink1 is less than 0" << endl;
    errorInFile = true;
  }else if(settings.tripleLink1 < settings.tripleLink2){
    cout << "ERROR: tripleLink1 is smaller than tripleLink2" << endl;
    errorInFile = true;
  }else{
    cout << "tripleLink1 is " << settings.tripleLink1 << endl;
  }
  if(settings.tripleLink2 < 0){
    cout << "ERROR: tripleLink2 is less than 0" << endl;
    errorInFile = true;
  }else if(settings.tripleLink2 < settings.tripleLink3){
    cout << "ERROR: tripleLink2 is smaller than tripleLink3" << endl;
    errorInFile = true;
  }else{
    cout << "tripleLink2 is " << settings.tripleLink2 << endl;
  }
  if(settings.tripleLink3 < 0){
    cout << "ERROR: tripleLink3 is smaller than 0" << endl;
    errorInFile = true;
  }else{
    cout << "tripleLink3 is " << settings.tripleLink3 << endl;
  }
  
  if(errorInFile) exit(1);

  
  return settings;
}


void printClusters(vector< graph<geneData, f64>* > clusters){
  for(size_t i = 0; i < clusters.size(); i++){
    cout << "Cluster " << (i+1) << ": " << endl;
    for(size_t j = 0; j < clusters[i]->getNumVertexes(); j++){
      cout << clusters[i]->getVertexes()[j]->value.getName() << endl;
    }
    cout << endl;
  }
}


template <typename T> vector<T> range(const vector<T> &in, 
                                csize_t &start, csize_t &end){
  vector<T> tr;
  for(size_t i = start; i < end; i++) tr.push_back(in[i]);
  return tr;
}


struct correlationMatrix sortWeights(struct correlationMatrix &protoGraph, 
                                                          cu8 keepTopN){
  pthread_t *workers;
  struct addTopEdgesHelperStruct *toSend;
  
  cu32 numCPUs = thread::hardware_concurrency() < protoGraph.labels.size() ? 
                            thread::hardware_concurrency() : protoGraph.labels.size();
  
  //cu32 numCPUs = 1;
                            
  toSend = (struct addTopEdgesHelperStruct*) malloc(sizeof(*toSend) * numCPUs);
  
  workers = (pthread_t*) malloc(sizeof(*workers) * numCPUs);
  
  u32 anItr = 0;
  while(anItr < numCPUs){
    toSend[anItr].protoGraph = &protoGraph;
    toSend[anItr].startIndex = (anItr * protoGraph.labels.size())/numCPUs;
    toSend[anItr].endIndex =   ((anItr+1) * protoGraph.labels.size())/numCPUs;
    toSend[anItr].keepTopN = keepTopN;
    //cout << "Dispatching worker " << anItr << " to prune entries " 
    //     << toSend[anItr].startIndex << " to " << toSend[anItr].endIndex 
    //     << endl;
    pthread_create(&workers[anItr], NULL, addTopEdgesHelper, &toSend[anItr]);
    //addTopEdgesHelper(&toSend[anItr]);
    anItr++;
  }
  
  void *toIgnore;
  for(anItr = 0; anItr < numCPUs; anItr++)
    pthread_join(workers[anItr], &toIgnore);
  
  free(toSend);
  free(workers);
  
  return protoGraph;
}


void *addTopEdgesHelper(void *protoArgs){
  struct addTopEdgesHelperStruct *args;
  pair<size_t, f64> *workSpace, *sortSpace;
  
  args = (struct addTopEdgesHelperStruct*) protoArgs;
  
  struct correlationMatrix *protoGraph = args->protoGraph;
  csize_t startIndex = args->startIndex;
  csize_t endIndex = args->endIndex;
  cu8 keepTopN = args->keepTopN;
  
  sortSpace = (pair<size_t, f64>*) malloc(sizeof(*sortSpace) * protoGraph->labels.size());
  
  for(size_t i = startIndex; i < endIndex; i++){
    size_t count[2] = {0, 0};
    pair<size_t, f64> swap;
    queue<size_t> indiciesOfInterest;
    queue<size_t> IOISwap;
    size_t size;
    
    size = protoGraph->colSize[i];
    //printf("%lu copy data of size %lu\n", i, size);  fflush(stdout);
    
    //copy
    //TODO -- make this so that memcpy can be used
    workSpace = protoGraph->matrix[i];
    if(NULL == workSpace){
      continue;
    }
    
    //printf("%lu presort\n", i);  fflush(stdout);
    
    //sort
    //Quickmerge
    for(size_t j = 0; j < size-1; j++){
      if(workSpace[j].second < workSpace[j+1].second) count[0]++;
      else if(workSpace[j].second > workSpace[j+1].second) count[1]++;
    }
    
    if(count[0] > count[1])
      for(size_t j = 0; j <= size/2; j++){
        swap = workSpace[j];
        workSpace[j] = workSpace[(size-1) - j];
        workSpace[(size - 1) - j] = swap;
      }
    
    indiciesOfInterest.push(0);
    for(size_t j = 1; j < size; j++)
      if(workSpace[j-1].second < workSpace[j].second)
        indiciesOfInterest.push(j);
    indiciesOfInterest.push(size);
    

    
    //printf("%lu sort\n", i);  fflush(stdout);
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
        //printf("\t%lu", leftStart);
        
        bool shouldExit = false;
        if(leftPtr >= leftEnd){
          printf("ERROR: indexes shuffled!  leftPtr (%lu) >= leftEnd (%lu)\n", leftPtr, leftEnd);
          shouldExit = true;
        }
        if(rightPtr >= rightEnd){ 
          printf("ERROR: indexes shuffled!  rightPtr (%lu) >= rightEnd (%lu)\n", rightPtr, rightEnd);
          shouldExit = true;
        }
        if(shouldExit) free(NULL);
        
        while(leftPtr < leftEnd && rightPtr < rightEnd){
          if(workSpace[leftPtr].second > workSpace[rightPtr].second)
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
      //  printf("\t%lu", indiciesOfInterest.front());
        indiciesOfInterest.pop();
      }
      IOISwap.push(size);
      //printf("\t%lu", size);
      
      indiciesOfInterest = IOISwap;
      while(!IOISwap.empty()) IOISwap.pop();
    }
    
    void *tmpPtr =  realloc(protoGraph->matrix[i], keepTopN * sizeof(*protoGraph->matrix[i]));
    protoGraph->matrix[i] = (pair<size_t, f64>*) tmpPtr;
    
  }
  
  free(sortSpace);
  
  return NULL;
}


graph<geneData, f64>* constructGraph(const struct correlationMatrix &protoGraph){
  graph<geneData, f64>* tr;
  size_t edgeCount;
  
  tr = new graph<geneData, f64>();
  tr->hintNumVertexes(protoGraph.matrix.size());
  
  edgeCount = 0;
  for(size_t i = 0; i < protoGraph.matrix.size(); i++)
    edgeCount += protoGraph.colSize[i];
  tr->hintNumEdges(edgeCount);
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    geneData tmpL;
    tmpL.construct().setName(protoGraph.labels[i]);
    vertex<geneData, f64> *leftVert = tr->addVertex(tmpL);
    
    for(size_t j = 0; j < protoGraph.colSize[i]; j++){
      geneData tmpR;
      tmpR.construct().setName(protoGraph.labels[j]);
      tr->addEdge(leftVert , tr->addVertex(tmpR), protoGraph.matrix[i][j].second);
    }
  }
  
  return tr;
}
