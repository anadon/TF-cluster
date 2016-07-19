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
  
  cerr << "topPick is " << (short) settings.keepTopN << endl;
  
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


void printClusters(queue< queue<size_t> > clusters, const vector<string> &names){
  for(size_t i = 0; !clusters.empty(); i++){
    cout << endl << "Cluster " << (i+1) << " : ";
    while(!clusters.front().empty()){
      cout << "\t" << names[clusters.front().front()];
      clusters.front().pop();
    }
    cout << endl;
    clusters.pop();
  }
}


/*template <typename T> vector<T> range(const vector<T> &in, 
                                csize_t &start, csize_t &end){
  vector<T> tr;
  for(size_t i = start; i < end; i++) tr.push_back(in[i]);
  return tr;
}*/


struct correlationMatrix sortWeights(struct correlationMatrix &protoGraph, 
                                                          cu8 keepTopN){
  pthread_t *workers;
  struct addTopEdgesHelperStruct *toSend;
  
  csize_t numCPUs = thread::hardware_concurrency() < protoGraph.labels.size() ? 
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
    pthread_create(&workers[anItr], NULL, addTopEdgesHelper, &toSend[anItr]);
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
    
    workSpace = protoGraph->matrix[i];
    if(NULL == workSpace) continue;
    
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
        indiciesOfInterest.pop();
      }
      IOISwap.push(size);
      
      indiciesOfInterest = IOISwap;
      while(!IOISwap.empty()) IOISwap.pop();
    }
    
    if(keepTopN < size){
      protoGraph->colSize[i] = keepTopN;
      void *tmpPtr =  realloc(protoGraph->matrix[i], protoGraph->colSize[i] * sizeof(*protoGraph->matrix[i]));
      protoGraph->matrix[i] = (pair<size_t, f64>*) tmpPtr;
    }
    
  }
  
  free(sortSpace);
  
  return NULL;
}


graph<geneData, f64>* constructGraph(const struct correlationMatrix &protoGraph){
  graph<geneData, f64>* tr;
  size_t edgeCount;
  
  tr = new graph<geneData, f64>();
  tr->hintNumVertexes(protoGraph.matrix.size());
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    vertex<geneData, double> *target = tr->addVertex(geneData(i));
    target->hintNumEdges(protoGraph.colSize[i]);
  }
  
  edgeCount = 0;
  for(size_t i = 0; i < protoGraph.matrix.size(); i++)
    edgeCount += protoGraph.colSize[i];
  tr->hintNumEdges(edgeCount);
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    vertex<geneData, f64> *leftVert = tr->getVertexForValue(geneData(i));
    for(size_t j = 0; j < protoGraph.colSize[i]; j++){
      vertex<geneData, f64> *rightVert = tr->getVertexForValue(geneData(protoGraph.matrix[i][j].first));
      if(leftVert == rightVert) raise(SIGABRT);
      tr->addEdge(leftVert , rightVert, protoGraph.matrix[i][j].second);
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
