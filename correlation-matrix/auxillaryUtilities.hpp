#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

#include <vector>
#include <string>


#include "statistics.hpp"


using std::vector;
using std::string;
using std::size_t;


struct config{
  string geneListFile;
  string expressionFile;
  u8 topPick;
  size_t kickSize;  
  double tripleLink1;
  double tripleLink2;
  double tripleLink3;
};


struct loadFromFileReturn{
  f64 **matrix;
  size_t numRows, numCols;
  vector<string> genes;
};


struct correlationMatrix{
  f64 *UDMatrix; // 1D upper diagonal matrix, omit x=y diagonal
  size_t numElements;
};


struct corrHelpStruct{
  size_t side;
  size_t start;
  size_t end;
  double **matrix;
  size_t vecLeng;
  double *results;
};




vector<string> tokenizeString(string toTokenize,
                                    const char *delimitors = " \t,\n");

int verifyInput(int argc, char **argv, string geneListFile, string expressionFile);


struct config loadConfig();


struct loadFromFileReturn loadFromFile(string geneListFile, 
                                                string expressionFile);

void convertCoeffToSigmaValue(f64 *UDMatrix, csize_t bSize);

void printEdges(const struct loadFromFileReturn &fileData, 
                              const struct correlationMatrix &corrData);

template <typename T> vector<T> range(const vector<T> &in, 
                                csize_t &start, csize_t &end);

double covariance(const vector<double> &left, 
                    const vector<double> &right, const double &leftSum, 
                                                const double &rightSum);

double* centerMean(const double *array, csize_t &size);

void *correlationHelper(void *protoArgs);

struct correlationMatrix calculateCorrelationMatrix(const struct loadFromFileReturn &input);


#endif