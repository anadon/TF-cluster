#ifndef GENEDATA_HPP
#define GENEDATA_HPP

#include <string>
#include <unordered_map>


using std::string;


class geneData{
  string *name;
  public:
  bool threeSigmaLink;
  bool twoSigmaLink;
  bool oneSigmaLink;
  
  geneData();
  
  geneData construct();
  
  void destroy();
  
  geneData setName(const string &geneName);
  
  string getName() const;
  
  geneData operator=(geneData const &other);
};

bool operator==(const geneData &lhs, const geneData &rhs);

namespace std{
  template <> struct hash<geneData> {
    std::size_t operator()(geneData const& toHash) const {
      return std::hash<std::string>()(toHash.getName());
    }
  };
}

#endif