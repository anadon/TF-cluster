#ifndef GENEDATA_HPP
#define GENEDATA_HPP

#include <string>
#include <unordered_map>



class geneData{
  public:
  std::string name;
  bool threeSigmaLink;
  bool twoSigmaLink;
  bool oneSigmaLink;

  bool operator==(const geneData &other) const;
};

namespace std{
template <> struct hash<geneData> {
  size_t operator()(const geneData& toHash) const {
    return (hash<string>()(toHash.name));
  }
};
}

#endif