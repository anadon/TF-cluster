/*******************************************************************//**
         FILE:  tripleLink.hpp

  DESCRIPTION:  Public interface for using tripleLink clustering

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/
#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <vector>

#include "geneData.hpp"
#include "graph.hpp"

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
queue< queue<size_t> > tripleLink(graph<geneData, double> *geneNetwork,
                    cf64 threeSigma, cf64 twoSigma);

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif