#include "geneData.hpp"

bool geneData::operator==(const geneData &other) const {
  return (name == other.name);
}