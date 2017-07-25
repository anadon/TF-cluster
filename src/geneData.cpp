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
