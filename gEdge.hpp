#ifndef GEDGE_HPP
#define GEDGE_HPP

#include <unistd.h>

#include "auxillaryUtilities.hpp"


class gNode;

class gEdge{
  public:
  size_t leftEdgeIndex, rightEdgeIndex, edgeID;
  gNode *left, *right;
  union{
    double coeff;
    u8 sigma;
  };

  gEdge(gNode *newLeft, gNode *newRight, size_t myID);

  ~gEdge();
};

#endif