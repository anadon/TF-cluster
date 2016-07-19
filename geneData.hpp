#ifndef GENEDATA_HPP
#define GENEDATA_HPP

#include <string>
#include <unordered_map>


using std::string;


class geneData{
  
  geneData();
  
  public:
  size_t nameIndex;
  bool threeSigmaLink;
  bool twoSigmaLink;
  bool oneSigmaLink;
  
  geneData(const size_t newNameIndex);
  
  geneData operator=(geneData const &other);
};

bool operator==(const geneData &lhs, const geneData &rhs);

namespace std{
  template <> struct hash<geneData> {
    std::size_t operator()(geneData const& toHash) const {
      return toHash.nameIndex;
    }
  };
}

#endif