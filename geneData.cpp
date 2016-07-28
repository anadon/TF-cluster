/***********************************************************************
         FILE:  geneData.cpp

  DESCRIPTION:  Container for use with triple-link clustering

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include "geneData.hpp"

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

geneData::geneData(const size_t newNameIndex):nameIndex(newNameIndex){}


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

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////