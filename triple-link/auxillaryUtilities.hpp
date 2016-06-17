#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>
#include <utility>

#include "geneData.hpp"
#include "graph.hpp"

using std::hash;
using std::pair;
using std::size_t;
using std::string;
using std::vector;


typedef unsigned char u8;
typedef unsigned int  u32;

typedef double f64;

typedef const unsigned char cu8;
typedef const unsigned int  cu32;
typedef const std::size_t csize_t;
typedef const double cf64;


struct config{
  string geneListFile;
  string expressionFile;
  u8 keepTopN;
  size_t kickSize;  
  f64 tripleLink1;
  f64 tripleLink2;
  f64 tripleLink3;
};


struct corrHelpStruct{
  size_t side;
  size_t start;
  size_t end;
  f64**matrix;
  size_t vecLeng;
  f64*results;
};


struct loadFromFileReturn{
  f64*UDMatrix; // upper diagonal matrix
  vector<string> genes;
};


struct upperDiagonalMatrix{
  vector<string> labels;
  unordered_map<string, size_t> labelLookup;
  vector< vector< pair<string, f64> > > matrix;
};


struct addTopEdgesHelperStruct{
  struct upperDiagonalMatrix &protoGraph;
  u8 keepTopN;
  size_t startIndex;
  size_t endIndex;
};


vector<string> tokenizeString(string toTokenize,
                                    const char *delimitors = " \t,\n");

int verifyInput(int argc, char **argv, const string geneListFile);


struct config loadConfig();


struct loadFromFileReturn loadFromFile(const string geneListFile);


//TODO make multithreading safe
//Quick merge, keep upper values
//this is not multithreading safe, but can be made safe.
void prune(vertex<geneData, f64> *toPrune, u8 keepTopN = 100);


void quickMergeEdges(vertex<geneData, f64> *toPrune,
                                                    csize_t size);


void mergeHelper(edge<geneData, f64> **toSort, 
                        csize_t leftIndex, csize_t rightIndex, 
                                                csize_t endIndex);


void pruneGraph(graph<geneData, f64> *geneNetwork, u8 keepTopN);

void convertCoeffToSigmaValue(struct loadFromFileReturn &fileData);

void removeLowEdges(graph<geneData, f64> *geneNetwork, cf64 &cutOff);
                                                  
void removeWeakVerticies(graph<geneData, f64> *geneNetwork);

void printClusters(vector< graph<geneData, f64>* > clusters);

void printEdges(graph<geneData, f64> *corrData);

template <typename T> vector<T> range(const vector<T> &in, 
                                csize_t &start, csize_t &end);

struct upperDiagonalMatrix sortWeights(struct upperDiagonalMatrix &protoGraph, cu8 keepTopN);

void *addTopEdgesHelper(void *protoArgs);

struct upperDiagonalMatrix loadMatrix(cf64 &minValue);

#endif