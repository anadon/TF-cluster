#include "geneData.hpp"


geneData::geneData(){
  name = "INVALID";
  threeSigmaLink = twoSigmaLink = oneSigmaLink = false;
}


geneData::geneData(string geneName){
  name = geneName;
  threeSigmaLink = twoSigmaLink = oneSigmaLink = false;
}


bool geneData::operator==(const geneData &other) const {
  return (name == other.name);
}