#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
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
  
  
  exprCoeffReader.open(expressionFile);
  
  if(!exprCoeffReader.good()){
    cout << "could not access \"" << expressionFile << "\"" << endl;
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
        cout << "Error: Line " << i << " Token " << j << 
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
  vector< vector<double> > corrMatrix;
  string tmp, headerGeneName, line;
  ifstream geneListReader, exprCoeffReader;
  size_t iterations;
  double *results;
  unsigned int numCPUs;
  pthread_t *workers;
  struct corrHelpStruct *instructions;
  double **rawMatrix;
  double *tmpArr;
  int *toIgnore;
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
  
  
  exprCoeffReader.open(expressionFile);
  
  while(exprCoeffReader.good()){
    getline(exprCoeffReader, line);
    if(string("") != line){
      vector<string> toParse = tokenizeString(line);
      rowHeaders.push_back(toParse[0]);
      vector<double> row;
      for(size_t i = 1; i < toParse.size(); i++){
        row.push_back(atof(toParse[i].c_str()));
      }
      corrMatrix.push_back(row);
    }
  }
  exprCoeffReader.close();

  const size_t rowSize = corrMatrix[0].size();
  
  rawMatrix = (double**) malloc( sizeof(*rawMatrix) * corrMatrix.size());
  tmpArr = (double*) malloc( sizeof(**rawMatrix) * rowSize);
  
  for(size_t i = 0; i < corrMatrix.size(); i++){
    for(size_t j = 0; j < rowSize; j++)
      tmpArr[j] = corrMatrix[i][j];
    rawMatrix[i] = centerMean(tmpArr, rowSize);
  }
  free(tmpArr);

  //for(size_t i = 0; i < rowHeaders.size(); i++)
  //  geneNetwork->addVertex(geneData(rowHeaders[i]));
  
  iterations = ((corrMatrix.size() * corrMatrix.size()) - corrMatrix.size()) / 2;
  results = (double*) malloc(sizeof(*results) * iterations);
  //NOTE: this is for debugging -- remove once final
  memset(results, -1.0, sizeof(*results) * iterations);
  
  //TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  numCPUs = std::thread::hardware_concurrency();
  numCPUs = numCPUs < iterations ? numCPUs : iterations;
  
  workers = (pthread_t*) malloc(sizeof(*workers) * numCPUs);
  instructions = (struct corrHelpStruct*) malloc(sizeof(*instructions) * numCPUs);
  
  for(size_t i = 0; i < numCPUs; i++){
    instructions[i].side = corrMatrix.size();
    instructions[i].start = i * (iterations / numCPUs);
    instructions[i].end = (i+1) * (iterations / numCPUs);
    instructions[i].matrix = rawMatrix;
    instructions[i].vecLeng = rowSize;
    instructions[i].results = results;
    pthread_create(&workers[i], NULL, correlationHelper, &instructions[i]);
  }
  
  for(size_t i = 0; i < numCPUs; i++)
    pthread_join(workers[i], (void**) &toIgnore);
  
  free(workers);
  free(instructions);
  for(size_t i = 0; i < corrMatrix.size(); i++)
    free(rawMatrix[i]);
  free(rawMatrix);
  
  tr.UDMatrix = results;
  tr.genes = rowHeaders;
  
  return tr;
}


void *correlationHelper(void *protoArgs){
  size_t i, j, itr;
  const struct corrHelpStruct *args = (struct corrHelpStruct*) protoArgs;
  
  itr = args->start;
  
  //needed because of triangular addressing
  i = args->side - sqrtl((args->side * args->side) - 2*itr);//y
  j = itr - i*(args->side - i/2.0);//x
  
  //primer loop
  for(;j < args->side && itr < args->end; j++)
    args->results[itr++] = getCenteredCorrelation(args->matrix[i], args->matrix[j], args->vecLeng);
  
  
  //main loop
  for(i++;itr < args->end; i++){
    for(j = i+1; j < args->side && itr < args->end; j++){
      args->results[itr++] = getCenteredCorrelation(args->matrix[i], args->matrix[j], args->vecLeng);
    }
  }
  
  return NULL;
}


void pruneGraph(graph<geneData, double> *geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    quickMergeEdges(&geneNetwork->vertexArray[i], geneNetwork->vertexArray[i].numEdges);
  }

  for(size_t i = 0; i < geneNetwork->numVertexes; i++){
    for(size_t j = geneNetwork->vertexArray[i].numEdges-1; j >= keepTopN; j--)
      geneNetwork->removeEdge(geneNetwork->vertexArray[i].edges[j]);
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
  
  free(sortSpace);
}


//WARNING: This function assumes a mean of 0
double calculateCenteredSigma(const double *array, const size_t size){
  double sumOfSquareDifferences;

  sumOfSquareDifferences = 0;
  for(size_t i = 0; i < size; i++)
    sumOfSquareDifferences += (array[i] * array[i]);

  return sqrt(sumOfSquareDifferences / size);
}


/*double calculateSigma(const double *array, const size_t size){
  double mean;
  double meanOfSquareDifferences;

  mean = 0;
  for(size_t i = 0; i < size; i++)
    mean += array[i];
  mean /= numEdges;

  meanOfSquareDifferences = 0;
  for(size_t i = 0; i < numEdges; i++){
    double difference = links[i].weight - mean;
    meanOfSquareDifferences += (difference * difference);
  }
  meanOfSquareDifferences /= numEdges;

  return sqrt(meanOfSquareDifferences);
}*/


void convertCoeffToSigmaValue(struct loadFromFileReturn &fileData){
  double sigma;
  size_t size;
  
  size = fileData.genes.size();
  size = (size * (size-1))/2;
  
  sigma = calculateCenteredSigma(fileData.UDMatrix, size);
  
  for(size_t i = 0; i < size; i++)
    fileData.UDMatrix[i] /= sigma;
}


void removeLowEdges(graph<geneData, double> *geneNetwork, 
                                                  const double &cutOff){
  for(size_t i = 0; i < geneNetwork->numEdges; i++)
    if(geneNetwork->edgeArray[i].weight < cutOff)
      geneNetwork->removeEdge(&geneNetwork->edgeArray[i]);
}


void removeWeakVerticies(graph<geneData, double> *geneNetwork){
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
        double threeSigmaReq, twoSigmaReq, oneSigmaReq;
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
  
  settings.topPick = 100;
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
    cout << "ERROR: geneList file not set" << endl;
    errorInFile = true;
  }
  if(settings.expressionFile != ""){
    cout << "expression is in \"" << settings.expressionFile << "\"" << endl;
  }else{
    cout << "ERROR: expression file not set" << endl;
    errorInFile = true;
  }
  
  cout << "topPick is " << (short) settings.topPick << endl;
  
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


void printClusters(vector< graph<geneData, double>* > clusters){
  for(size_t i = 0; i < clusters.size(); i++){
    cout << "Cluster " << (i+1) << ": " << endl;
    for(size_t j = 0; j < clusters[i]->numVertexes; j++){
      cout << clusters[i]->getVertexes()[j].value.name << endl;
    }
    cout << endl;
  }
}


void printEdges(graph<geneData, double> *corrData){
  string left, right;
  double weight;
  
  cout << "Here are all the edges in the graph:" << endl;
  
  for(size_t i = 0; i < corrData->numEdges; i++){
    left = corrData->getEdges()[i].left->value.name;
    right = corrData->getEdges()[i].right->value.name;
    weight = corrData->getEdges()[i].weight;
    cout << left << " <--" << weight << "--> " << right << endl;
  }
  
  cout << endl;
}


template <typename T> vector<T> range(const vector<T> &in, 
                                const size_t &start, const size_t &end){
  vector<T> tr;
  for(size_t i = start; i < end; i++) tr.push_back(in[i]);
  return tr;
}


double covariance(const vector<double> &left, 
                    const vector<double> &right, const double &leftSum,
                                                const double &rightSum){
  double tr = 0, mixedSum = 0;
  const size_t size = left.size();
  if(size != right.size()){
    exit(-EINVAL);
  }
  
  for(size_t i = 0; i < size; i++)  mixedSum += (left[i] * right[i]);
  
  tr =  (mixedSum / size) / ((leftSum / size)*(rightSum / size));
  
  return tr;
}


double getMean(const double *array, const size_t &size){
  return getSum(array, size) / size;
}


double getStandardDeviation(const double *array, const size_t &size){
  double tr = 0;
  const double mean = getMean(array, size);
  double tmp;
  
  for(size_t i = 0; i < size; i++){
    tmp = array[i] - mean;
    tr += (tmp * tmp);
  }
  tr = sqrt(tr / size);
  
  return tr;
}


double getSum(const double *array, const size_t &size){
  double tr = 0;
  for(size_t i = 0; i < size; i++) tr += array[i];
  return tr;
}


double getSumOfSquares(const double *array, const size_t &size){
  return getSumOfMultipliedArrays(array, array, size);
}


double getSumOfMultipliedArrays(const double *left, const double *right, 
                                                    const size_t &size){
  double tr = 0;
  for(size_t i = 0; i < size; i++) tr += (left[i] * right[i]);
  return tr;
}


//WARNING: This functions assumes left and right have a mean of 0
double getCenteredCorrelation(const double *left, const double *right, 
                                                    const size_t &size){
  
  const double aSumOfSquares = getSumOfSquares(left, size);
  const double bSumOfSquares = getSumOfSquares(right, size);
  const double abCrossSum = getSumOfMultipliedArrays(left, right, size);
  
  return abCrossSum / sqrt(aSumOfSquares * bSumOfSquares);
}


double* centerMean(const double *array, const size_t &size){
  double *tr;
  
  const double mean = getMean(array, size);
  tr = (double*) malloc(sizeof(*array) * size);
  
  for(size_t i = 0; i < size; i++)
    tr[i] = array[i] - mean;
  
  return tr;
}


void addTopEdges(graph<geneData, double> *geneNetwork, 
                struct loadFromFileReturn &fileData, const u8 keepTopN){
  
  const u32 numCPUs = thread::hardware_concurrency() < fileData.genes.size() ? 
                            thread::hardware_concurrency() : fileData.genes.size();
                            
  struct addTopEdgesHelperStruct *toSend;
  toSend = (struct addTopEdgesHelperStruct*) malloc(sizeof(*toSend) * numCPUs);
  
  pthread_t *workers;
  workers = (pthread_t*) malloc(sizeof(*workers) * numCPUs);
  
  u32 anItr = 0;
  while(anItr < numCPUs){
    toSend[anItr].size = fileData.genes.size();
    toSend[anItr].geneNames = &fileData.genes;
    toSend[anItr].UDMatrix = fileData.UDMatrix;
    toSend[anItr].geneNetwork = geneNetwork;
    toSend[anItr].startIndex = (anItr * toSend[anItr].size)/numCPUs;
    toSend[anItr].endIndex =   (anItr * toSend[anItr].size)/numCPUs;
    toSend[anItr].keepTopN = keepTopN;
    pthread_create(&workers[anItr], NULL, addTopEdgesHelper, &toSend[anItr]);
    anItr++;
  }
  
  void *toIgnore;
  for(anItr = 0; anItr < numCPUs; anItr++)
    pthread_join(workers[anItr], &toIgnore);
  
  
}


void *addTopEdgesHelper(void *protoArgs){
  struct addTopEdgesHelperStruct *args;
  pair<string*, double> *workSpace, *sortSpace;
  
  args = (struct addTopEdgesHelperStruct*) protoArgs;
  
  const size_t startIndex = args->startIndex;
  const size_t endIndex = args->endIndex;
  const size_t size = args->size;
  workSpace = (pair<string*, double>*) malloc(sizeof(*workSpace) * (size - 1));
  sortSpace = (pair<string*, double>*) malloc(sizeof(*workSpace) * (size - 1));
  
  size_t itr;
  for(size_t i = startIndex; i < endIndex; i++){
    itr = 0;
    //copy
    for(size_t j = 0; j < i; j++){ //i-> x, j->y
      workSpace[itr].first = &(*args->geneNames)[j];
      workSpace[itr].second = args->UDMatrix[i + (j*size) - ((j*(j+1))/2)];
      itr++;
    }
    for(size_t j = i+1; i < size; j++){//i->y, j->x
      workSpace[itr].first = &(*args->geneNames)[j];
      workSpace[itr].second = args->UDMatrix[j + (i*size) - ((i*(i+1))/2)];
      itr++;
    }
    
    //sort
    //Quickmerge
    size_t count[2] = {0, 0};
    for(size_t j = 0; j < size-2; j++){
      if(workSpace[j].second < workSpace[j+1].second) count[0]++;
      else if(workSpace[j].second > workSpace[j+1].second) count[1]++;
    }
    
    string *swapTmpS;
    double  swapTmpD;
    if(count[0] > count[1])
      for(size_t j = 0; j < (size-1)/2; j++){
        swapTmpS = workSpace[j].first;
        swapTmpD = workSpace[j].second;
        workSpace[j].first = workSpace[(size - 2) - j].first;
        workSpace[j].second = workSpace[(size - 2) - j].second;
        workSpace[(size - 2) - j].first = swapTmpS;
        workSpace[(size - 2) - j].second = swapTmpD;
      }
    
    queue<size_t> indiciesOfInterest;
    queue<size_t> IOISwap;
    
    indiciesOfInterest.push(0);
    for(size_t j = 1; j < size - 1; j++)
      if(workSpace[j-1].second < workSpace[j].second)
        indiciesOfInterest.push(j);
    indiciesOfInterest.push(size - 1);
    
    while(indiciesOfInterest.size() > 2){
      while(indiciesOfInterest.size() > 1){
        size_t leftStart, leftPtr, leftEnd, rightPtr, rightEnd;
        itr = 0;
        leftStart = leftPtr = indiciesOfInterest.front(); indiciesOfInterest.pop();
        leftEnd = rightPtr = indiciesOfInterest.front();  indiciesOfInterest.pop();
        rightEnd = indiciesOfInterest.front();
        IOISwap.push(leftStart);
        
        while(leftPtr < leftEnd && rightPtr < rightEnd){
          if(workSpace[leftPtr].second > workSpace[rightPtr].second)
            sortSpace[itr++] = workSpace[leftPtr++];
          else
            sortSpace[itr++] = workSpace[rightPtr++];
        }
        while(leftPtr < leftEnd)   sortSpace[itr++] = workSpace[leftPtr++];
        while(rightPtr < rightEnd) sortSpace[itr++] = workSpace[rightPtr++];
        
        const size_t copySize = (rightEnd - leftStart) * sizeof(*sortSpace);
        memcpy(&workSpace[leftStart], sortSpace, copySize);
      }
      indiciesOfInterest = IOISwap;
      indiciesOfInterest.push(size - 1);
    }
    
    
    vertex<geneData, double> *leftVertex, *rightVertex;
    leftVertex = args->geneNetwork->getVertexForValue(geneData((*args->geneNames)[i]));
    itr = 0;
    for(u8 j = 0; j < args->keepTopN; j++){
      rightVertex = args->geneNetwork->getVertexForValue(geneData(*workSpace[j].first));
      args->geneNetwork->addEdge(leftVertex, rightVertex, workSpace[j].second);
    }
    
  }
  
  free(workSpace);
  free(sortSpace);
  return NULL;
}
