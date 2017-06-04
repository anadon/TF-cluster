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

#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <graph.tpp>
#include <vector>

#include "geneData.hpp"

////////////////////////////////////////////////////////////////////////
//PUBLIC FUNCTION DECLARATION///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Clustering algorithm described in "TF-Cluster: A pipeline for
 * identifying functionally coordinated transcription factors via
 * network decomposition of the shared coexpression connectivity matrix
 * (SCCM)".  Returns a queue of queues, with each subqueue holding label
 * indexes that represent a particular cluster.
 *
 * @param[in,out] geneNetwork Graph of genes which are parsed with the
 *                            triple-link algorithm.  The graph is
 *                            modified by this operation.
 * @param[in] threeSigma High value used for strong edges in
 *                       triple-link.
 * @param[in] twoSigma Medium value used for edges in triple-link.
 **********************************************************************/
queue< queue<size_t> > tripleLink(graph<geneData, unsigned char> *geneNetwork,
                    const struct config &settings);

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
