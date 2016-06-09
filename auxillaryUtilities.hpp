#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>

#include "geneData.hpp"
#include "graph.hpp"

using std::vector;
using std::string;
using std::hash;
using std::size_t;


typedef unsigned char u8;
typedef unsigned int  u32;


struct config{
  string geneListFile;
  string expressionFile;
  u8 topPick;
  size_t kickSize;  
  double tripleLink1;
  double tripleLink2;
  double tripleLink3;
};


struct corrHelpStruct{
  size_t side;
  size_t start;
  size_t end;
  double **matrix;
  size_t vecLeng;
  double *results;
};


struct loadFromFileReturn{
  double *UDMatrix; // upper diagonal matrix
  vector<string> genes;
};


struct addTopEdgesHelperStruct{
  size_t startIndex;
  size_t endIndex;
  size_t size;
  vector<string> *geneNames;
  double *UDMatrix;
  u8 keepTopN;
  graph<geneData, double> *geneNetwork;
};


vector<string> tokenizeString(string toTokenize,
                                    const char *delimitors = " \t,\n");

int verifyInput(int argc, char **argv, string geneListFile, string expressionFile);


struct config loadConfig();


struct loadFromFileReturn loadFromFile(string geneListFile, 
                                                string expressionFile);


//TODO make multithreading safe
//Quick merge, keep upper values
//this is not multithreading safe, but can be made safe.
void prune(vertex<geneData, double> *toPrune, u8 keepTopN = 100);


void quickMergeEdges(vertex<geneData, double> *toPrune,
                                                    const size_t size);


void mergeHelper(edge<geneData, double> **toSort, 
                        const size_t leftIndex, const size_t rightIndex, 
                                                const size_t endIndex);


void pruneGraph(graph<geneData, double> *geneNetwork, u8 keepTopN);


double calculateCenteredSigma(const edge<geneData, double> * links, 
                                                const size_t numEdges);

void convertCoeffToSigmaValue(struct loadFromFileReturn &fileData);

void removeLowEdges(graph<geneData, double> *geneNetwork, 
                                                  const double &cutOff);
                                                  
void removeWeakVerticies(graph<geneData, double> *geneNetwork);

void printClusters(vector< graph<geneData, double>* > clusters);

void printEdges(graph<geneData, double> *corrData);

template <typename T> vector<T> range(const vector<T> &in, 
                                const size_t &start, const size_t &end);

double covariance(const vector<double> &left, 
                    const vector<double> &right, const double &leftSum, 
                                                const double &rightSum);

double getMean(const double *array, const size_t &size);

double getStandardDeviation(const double *array, const size_t &size);

double getSum(const double *array, const size_t &size);

double getSumOfSquares(const double *array, const size_t &size);

double getSumOfMultipliedArrays(const double *left, const double *right, 
                                                    const size_t &size);

double getCenteredCorrelation(const double *left, const double *right, 
                                                    const size_t &size);

double* centerMean(const double *array, const size_t &size);

void *correlationHelper(void *protoArgs);

void addTopEdges(graph<geneData, double> *geneNetwork, 
          struct loadFromFileReturn &fileData, const u8 keepTopN = 100);

void *addTopEdgesHelper(void *protoArgs);

#endif