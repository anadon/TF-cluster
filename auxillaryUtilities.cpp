#include <string.h>
#include <math.h>


#include "auxillaryUtilities.hpp"
#include "graph.t.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"


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

  if( argc > 2 ) return EINVAL;

  if(argc == 2) freopen(argv[1], "r", stdin);

  return 0;
}


graph<struct geneData, double> loadFromFile(FILE *input){
  graph<struct geneData, double> toReturn;
  char *line = NULL;
  vector<string> tokens;

  getline(&line, NULL, stdin);
  tokens = tokenizeString(line);

  for(size_t i = 0; i < tokens.size(); i++) toReturn.addVertex({tokens[i], 0});

  while(feof(stdin)){

    getline(&line, NULL, stdin);
    tokens = tokenizeString(line);

    for(size_t i = 0; i < tokens.size(); i+=2){
      geneData toMatch = {tokens[i], 0};
      toReturn.addEdge(toReturn.vertexArray[i], toReturn.vertexArray[toReturn.geneNameToNodeID.find(toMatch)->second], atof(tokens[i+1].c_str()));
    }
  }

  return toReturn;
}


void pruneGraph(graph<struct geneData, double> geneNetwork, u8 keepTopN){

  for(size_t i = 0; i < geneNetwork.numNodes; i++){
    quickMergeEdges(geneNetwork.vertexArray[i], geneNetwork.vertexArray[i]->numEdges);
  }

  for(size_t i = 0; i < geneNetwork.numNodes; i++){
    for(size_t j = geneNetwork.vertexArray[i]->numEdges-1; j >= keepTopN; j--)
      geneNetwork.removeEdge(geneNetwork.vertexArray[i]->edges[j]);
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