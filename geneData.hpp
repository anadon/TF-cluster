/***********************************************************************
         FILE:  geneData.hpp

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
#ifndef GENEDATA_HPP
#define GENEDATA_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <unordered_map>

////////////////////////////////////////////////////////////////////////
//CLASS DEFINITION//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  
 **********************************************************************/
class geneData{

  geneData();
  
  public:
  size_t nameIndex;
  bool threeSigmaLink;
  bool twoSigmaLink;
  bool oneSigmaLink;
  
  
/*******************************************************************//**
 *  Make a new geneData with the index of the name it represents.
 **********************************************************************/
  geneData(const size_t newNameIndex);
  
  
/*******************************************************************//**
 *  Custom copy operator.
 * 
 * @param[in] other geneData to copy nameIndex from.
 **********************************************************************/
  geneData operator=(geneData const &other);
};


/*******************************************************************//**
 *  Custom equality operator; compared based only on nameIndex.
 * 
 * @param[in] lhs Left of equality comparison
 * @param[in] rhs Right of equality comparison.
 **********************************************************************/
bool operator==(const geneData &lhs, const geneData &rhs);

namespace std{
  template <> struct hash<geneData> {
    
/*******************************************************************//**
 *  Add a hash specialization to handle geneData in the standard 
 * namespace.  Hash is given as the name index, since it can be assumed
 * this is unique within a given graph.
 * 
 * @param[in] toHash geneData that will give it's nameIndex as a hash.
 **********************************************************************/
    std::size_t operator()(geneData const& toHash) const {
      return toHash.nameIndex;
    }
  };
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif