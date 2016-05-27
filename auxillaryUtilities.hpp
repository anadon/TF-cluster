#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>

#include "geneData.hpp"
#include "graph.hpp"
//*
using std::vector;
using std::string;
using std::hash;
using std::size_t;
//*/

typedef unsigned char u8;




struct config{
  string geneListFile;
  string expressionFile;
  u8 topPick;
  size_t kickSize;  
  double tripleLink1;
  double tripleLink2;
  double tripleLink3;
};


vector<string> tokenizeString(char *input,
                                      const char *delimitors = " \t,");

int verifyInput(int argc, char **argv);


struct config loadConfig();


void loadFromFile(graph<struct geneData, double> *geneNetwork, 
                          string geneListFile, string expressionFile);


//TODO make multithreading safe
//Quick merge, keep upper values
//this is not multithreading safe, but can be made safe.
void prune(vertex<struct geneData, double> *toPrune, u8 keepTopN = 100);


void quickMergeEdges(vertex<struct geneData, double> *toPrune,
                                                    const size_t size);


void mergeHelper(edge<struct geneData, double> **toSort, 
                        const size_t leftIndex, const size_t rightIndex, 
                                                const size_t endIndex);


void pruneGraph(graph<struct geneData, double> *geneNetwork, u8 keepTopN);


double calculateSigma(const edge<struct geneData, double> ** links, 
                                                const size_t numEdges);

void convertCoeffToSigmaValue(graph<struct geneData, double> *geneNetwork, 
                                                          double sigma);

void removeLowEdges(graph<struct geneData, double> *geneNetwork, 
                                                  const double &cutOff);
                                                  
void removeWeakVerticies(graph<struct geneData, double> *geneNetwork);

void simpleError(const char *message);

#endif