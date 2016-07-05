#include "geneData.hpp"


geneData::geneData(){
  threeSigmaLink = twoSigmaLink = oneSigmaLink = false;
  name = new string("");
}


geneData geneData::construct(){
  threeSigmaLink = twoSigmaLink = oneSigmaLink = false;
  name = new string("");
  return *this;
}


void geneData::destroy(){
  if(NULL == name) delete name;
}


geneData geneData::setName(const string &geneName){
  if(NULL == name){
    name = new string(geneName);
  }else{
    *name = geneName;
  }
  return *this;
}


std::string geneData::getName() const{
  return *name;
}


geneData geneData::operator=(geneData const &other){
  oneSigmaLink = other.oneSigmaLink;
  twoSigmaLink = other.twoSigmaLink;
  threeSigmaLink = other.threeSigmaLink;
  name = other.name;
  
  return *this;
}


bool operator==(const geneData &lhs, const geneData &rhs){
  return (lhs.getName() == rhs.getName());
}