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




void pruneGraph(graph<geneData, f64> *geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    quickMergeEdges(&geneNetwork->vertexArray[i], geneNetwork->vertexArray[i].numEdges);
  }

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    for(size_t j = geneNetwork->vertexArray[i].numEdges-1; j >= keepTopN; j--)
      geneNetwork->removeEdge(geneNetwork->vertexArray[i].edges[j]);
  }
}


void quickMergeEdges(vertex<geneData, f64> *toPrune,
                                                    csize_t size){
  size_t numRising;
  size_t i;

  numRising = 0;

  for(i = 0; i < size-1; i++){
    f64 prev, next;
    prev = toPrune->edges[i]->weight;
    next = toPrune->edges[i+1]->weight;
    if(prev < next) numRising++;
  }

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    edge<geneData, f64> *tmp;
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


void removeLowEdges(graph<geneData, f64> *geneNetwork, 
                                                  cf64&cutOff){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    if(geneNetwork->edgeArray[i].weight < cutOff)
      geneNetwork->removeEdge(&geneNetwork->edgeArray[i]);
}


void removeWeakVerticies(graph<geneData, f64> *geneNetwork){
  //first, we need to remove all nodes which do not have 3 available
  //links.  It's technically possible, and technically possible for new
  //ones to occur after a given removal so this is expensive...
  bool disconnectedVerticiesFound;
  do{
    disconnectedVerticiesFound = false;
    for(size_t i = 0; i < geneNetwork->numVertexes; i++){
      if(3 > geneNetwork->vertexArray[i].numEdges){
        geneNetwork->removeVertex(&geneNetwork->vertexArray[i]);
        disconnectedVerticiesFound = true;
      }else{
        f64 threeSigmaReq, twoSigmaReq, oneSigmaReq;
        threeSigmaReq = geneNetwork->vertexArray[i].edges[0]->weight;
        twoSigmaReq   = geneNetwork->vertexArray[i].edges[1]->weight;
        oneSigmaReq   = geneNetwork->vertexArray[i].edges[2]->weight;
        if(3.0 > threeSigmaReq || 2 > twoSigmaReq || 1 > oneSigmaReq){
          geneNetwork->removeVertex(&geneNetwork->vertexArray[i]);
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
    for(size_t j = 0; j < clusters[i]->numVertexes; j++){
      cout << clusters[i]->getVertexes()[j].value.name << endl;
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


struct upperDiagonalMatrix sortWeights(struct upperDiagonalMatrix &protoGraph, 
                                                          cu8 keepTopN){
  pthread_t *workers;
  struct addTopEdgesHelperStruct *toSend;
  
  //cu32 numCPUs = thread::hardware_concurrency() < protoGraph.labels.size() ? 
  //                          thread::hardware_concurrency() : protoGraph.labels.size();
  
  cu32 numCPUs = 1;
                            
  toSend = (struct addTopEdgesHelperStruct*) malloc(sizeof(*toSend) * numCPUs);
  
  workers = (pthread_t*) malloc(sizeof(*workers) * numCPUs);
  
  u32 anItr = 0;
  while(anItr < numCPUs){
    toSend[anItr].protoGraph = &protoGraph;
    toSend[anItr].startIndex = (anItr * protoGraph.labels.size())/numCPUs;
    toSend[anItr].endIndex =   ((anItr+1) * protoGraph.labels.size())/numCPUs;
    toSend[anItr].keepTopN = keepTopN;
    cout << "Dispatching worker " << anItr << " to prune entries " 
         << toSend[anItr].startIndex << " to " << toSend[anItr].endIndex 
         << endl;
    //pthread_create(&workers[anItr], NULL, addTopEdgesHelper, &toSend[anItr]);
    addTopEdgesHelper(&toSend[anItr]);
    anItr++;
  }
  
  //void *toIgnore;
  //for(anItr = 0; anItr < numCPUs; anItr++)
  //  pthread_join(workers[anItr], &toIgnore);
  
  free(toSend);
  free(workers);
  
  return protoGraph;
}


void *addTopEdgesHelper(void *protoArgs){
  struct addTopEdgesHelperStruct *args;
  pair<size_t, f64> *workSpace, *sortSpace;
  
  args = (struct addTopEdgesHelperStruct*) protoArgs;
  
  struct upperDiagonalMatrix *protoGraph = args->protoGraph;
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
    
    //printf("%lu presort\n", i);  fflush(stdout);
    
    //sort
    //Quickmerge
    for(size_t j = 0; j < size; j++){
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
      /*queue<size_t> DEBUG;
      while(indiciesOfInterest.size() > 1){
        size_t tmp = indiciesOfInterest.front();
        indiciesOfInterest.pop();
        if(tmp >= indiciesOfInterest.front()){
          printf("\non entry %lu, %lu >= %lu", i, tmp, indiciesOfInterest.front());
        }else{
          //printf("%lu < %lu\n", tmp, indiciesOfInterest.front());
          //printf(".");
          DEBUG.push(tmp);
        }
      }
      DEBUG.push(size);
      indiciesOfInterest = DEBUG;*/
      
      printf("\noutter\n");
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
        printf("\t%lu", leftStart);
        
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
        printf("\t%lu", indiciesOfInterest.front());
        indiciesOfInterest.pop();
      }
      IOISwap.push(size);
      printf("\t%lu", size);
      
      indiciesOfInterest = IOISwap;
      while(!IOISwap.empty()) IOISwap.pop();
    }
    
    void *tmpPtr =  realloc(protoGraph->matrix[i], keepTopN * sizeof(*protoGraph->matrix[i]));
    if(tmpPtr != protoGraph->matrix[i]) free(protoGraph->matrix[i]);
    protoGraph->matrix[i] = (pair<size_t, f64>*) tmpPtr;
    
  }
  
  free(sortSpace);
  
  return NULL;
}


struct upperDiagonalMatrix loadMatrix(cf64 &minValue){
  struct upperDiagonalMatrix tr;
  char left[128], right[128];
  string leftS, rightS;
  f64 weight;
  vector<size_t> columnSize;
  size_t leftIndex, rightIndex, numEntries;
  void *tmpPtr;
  
  numEntries = 0;
  
  while(!feof(stdin) && !ferror(stdin)){
    
    scanf("%127s %127s %lf", left, right, &weight);
    weight = std::abs(weight);//I think we care about absolute correlation
    if(weight < minValue) continue;
    
    leftS = string(left);
    rightS = string(right);
    
    if(tr.labelLookup.end() == tr.labelLookup.find(leftS)){
      leftIndex = tr.labels.size();
      tr.labelLookup.emplace(std::make_pair(leftS, leftIndex));
      tr.labels.push_back(leftS);
      numEntries++;
      
      leftIndex = tr.labels.size()-1;
      tr.matrix.push_back( (pair<size_t, f64>*) NULL);
      
      columnSize.push_back(0);
      tr.colSize.push_back(0);
    }else{
      leftIndex = tr.labelLookup[leftS];
    }
    
    if(tr.labelLookup.end() == tr.labelLookup.find(rightS)){
      rightIndex = tr.labels.size();
      tr.labelLookup.emplace(std::make_pair(rightS, rightIndex));
      tr.labels.push_back(rightS);
      numEntries++;
      
      rightIndex = tr.labels.size() -1;
      tr.matrix.push_back( (pair<size_t, f64>*) NULL );
      
      columnSize.push_back(0);
      tr.colSize.push_back(0);
    }else{
      rightIndex = tr.labelLookup[rightS];
    }
    
    
    
    if(columnSize[leftIndex] <= tr.colSize[leftIndex]){
      columnSize[leftIndex] = (tr.colSize[leftIndex]+1)*2;
      tmpPtr = realloc(tr.matrix[leftIndex], sizeof(pair<size_t, f64>) * columnSize[leftIndex] );
      if(NULL == tmpPtr) exit(-ENOMEM);
      tr.matrix[leftIndex] = (pair<size_t, f64>*) tmpPtr;
    }
    tr.matrix[leftIndex][tr.colSize[leftIndex]] = make_pair(rightIndex, weight);
    tr.colSize[leftIndex]++;
    
    if(columnSize[rightIndex] <= tr.colSize[rightIndex]){
      columnSize[rightIndex] = (tr.colSize[rightIndex]+1)*2;
      tmpPtr = realloc(tr.matrix[rightIndex], sizeof(pair<size_t, f64>) * columnSize[rightIndex]);
      if(NULL == tmpPtr) exit(-ENOMEM);
      tr.matrix[rightIndex] = (pair<size_t, f64>*) tmpPtr;
    }
    tr.matrix[rightIndex][tr.colSize[rightIndex]] = make_pair(leftIndex, weight);
    tr.colSize[rightIndex]++;
  }
  
  for(size_t i = 0; i < tr.colSize.size(); i++){
    tr.matrix[rightIndex] = (pair<size_t, f64>*) realloc(tr.matrix[i], tr.colSize[i] * sizeof(pair<size_t, f64>));
  }
  
  return tr;
}
