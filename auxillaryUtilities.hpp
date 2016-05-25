#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>

#include "graph.hpp"


using std::vector;
using std::string;
using std::hash;
using std::size_t;

typedef unsigned char u8;

struct geneData{
  string name;
  size_t clusterGroup;

  bool operator==(const geneData &other) const {
    return (name == other.name && clusterGroup == other.clusterGroup);
  }
};


namespace std{
  template <> struct hash<geneData> {
    size_t operator()(const geneData& toHash) const {
      return (hash<string>()(toHash.name));
    }
  };
}


vector<string> tokenizeString(char *input,
                                      const char *delimitors = " \t,");

int verifyInput(int argc, char **argv);


graph<struct geneData, double> loadFromFile(FILE *input);


//TODO make multithreading safe
//Quick merge, keep upper values
//this is not multithreading safe, but can be made safe.
void prune(vertex<struct geneData, double> *toPrune, u8 keepTopN = 100);


void quickMergeEdges(vertex<struct geneData, double> *toPrune,
                                                    const size_t size);


void mergeHelper(edge<struct geneData, double> **toSort, const size_t leftIndex,
                        const size_t rightIndex, const size_t endIndex);


void pruneGraph(graph<struct geneData, double> geneNetwork, u8 keepTopN);


double calculateSigma(const edge<struct geneData, double> ** links, const size_t numEdges);

#endif