#include "geneData.hpp"


geneData::geneData(const size_t newNameIndex):nameIndex(newNameIndex){
  //threeSigmaLink = twoSigmaLink = oneSigmaLink = false;
}


geneData geneData::operator=(geneData const &other){
  oneSigmaLink = other.oneSigmaLink;
  twoSigmaLink = other.twoSigmaLink;
  threeSigmaLink = other.threeSigmaLink;
  nameIndex = other.nameIndex;
  
  return *this;
}


bool operator==(const geneData &lhs, const geneData &rhs){
  return (lhs.nameIndex == rhs.nameIndex);
}