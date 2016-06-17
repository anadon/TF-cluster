#ifndef GENEDATA_HPP
#define GENEDATA_HPP

#include <string>
#include <unordered_map>


using std::string;


class geneData{
  public:
  string name;
  bool threeSigmaLink;
  bool twoSigmaLink;
  bool oneSigmaLink;
  
  geneData();
  
  geneData(string geneName);

  bool operator==(const geneData &other) const;
};

namespace std{
  template <> struct hash<geneData> {
    std::size_t operator()(geneData const& toHash) const {
      return std::hash<std::string>()(toHash.name);
    }
  };
}

#endif