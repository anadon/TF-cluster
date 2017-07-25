/*Copyright 2016-2017 Josh Marshall************************************/

/***********************************************************************
    This file is part of TF-Cluster.

    TF-Cluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TF-Cluster.  If not, see <http://www.gnu.org/licenses/>.
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
