
/*******************************************************************//**
         FILE:  auxillaryUtilities.cpp

  DESCRIPTION:  Miscelaneous functions used in TF-cluser

         BUGS:  Correlation values are larger than perl version
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

#include <iostream>

#include "diagnostics.hpp"
#include "edge.t.hpp"
#include "graph.t.hpp"
#include "vertex.t.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::cout;
using std::cerr;
using std::endl;
using std::queue;
using std::string;
using std::vector;

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void printClusters(queue< queue<size_t> > clusters,
                                            const vector<string> &TFs){
  for(size_t i = 0; !clusters.empty(); i++){
    cout << "cluster: " << (i+1) << endl;
    while(!clusters.front().empty()){
      cout << TFs[clusters.front().front()] << endl;
      clusters.front().pop();
    }
    clusters.pop();
  }
}


void printCoincidenceMatrix(UpperDiagonalSquareMatrix<u8> matrix, 
                                cu8 maxMatch, const vector<string> TFs){
  f64 **mtr;
  size_t n = matrix.getSideLength();
  mtr = (f64**) malloc(sizeof(*mtr) * n);
  for(size_t i = 0; i < n; i++)
    mtr[i] = (f64*) malloc(sizeof(**mtr) * n);

  for(size_t i = 0; i < n; i++){\
    mtr[i][i] = maxMatch;
    for(size_t j = i+1; j < n; j++){
      mtr[i][j] = mtr[j][i] = matrix.getValueAtIndex(i, j);
    }
  }

  for(size_t i = 0; i < n; i++){
    cout << TFs[i] << "\t";
    for(size_t j = 0; j < n; j++){
      cout << mtr[i][j] << "\t";
    }
    cout << endl;
    free(mtr[i]);
  }
  free(mtr);
  cout << endl << endl;
}


void printCorrelationMatrix(const CMF &protoGraph){
  //vector< pair<size_t, double>* > matrix;
  //unordered_map<string, size_t> labelLookup;
  //vector<string> labels;
  //vector<size_t> colSize;
  printf("\t");
  for(size_t i = 0; i < protoGraph.GeneLabels.size();i++)
    printf("%s\t", protoGraph.GeneLabels[i].c_str());


  for(size_t i = 0; i < protoGraph.TFLabels.size(); i++){
    printf("%s\t", protoGraph.TFLabels[i].c_str());
    for(size_t j = 0; j < protoGraph.GeneLabels.size(); j++){
      printf("%lf\t", protoGraph.fullMatrix[i][j]);
    }
    printf("\n");
  }

}


void graphvizRepresentation(graph<geneData, f64> *corrData, 
                                                vector<string> labels){
  cout << "digraph rep{" << endl;
  cout << endl;
  cout << "bgcolor=\"transparent\"" << endl;
  cout << "fontcolor = black" << endl;
  cout << endl;
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    cout << labels[corrData->getVertexes()[i]->value.nameIndex] << endl;
  }
  
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    vertex<geneData, f64> *targetV;
    targetV = corrData->getVertexes()[i];
    for(size_t j = 0; j < targetV->getNumEdges(); j++){
      edge<geneData, f64> *targetE;
      targetE = corrData->getVertexes()[i]->getEdges()[j];
      if(targetV == targetE->left){
        cout << labels[targetE->left->value.nameIndex] << " -> ";
        cout << labels[targetE->right->value.nameIndex] << " ;" << endl;
      }
    }
  }
  
  cout << "}";
}


void printEdgeWeights(graph<geneData, u8> *corrData){
  for(size_t i = 0; i < corrData->getNumEdges(); i++)
    cout << (int) corrData->getEdges()[i]->weight << endl;
}


void printProtoGraph(const CMF &toPrint){
  for(size_t i = 0; i < toPrint.numRows(); i++){
    for(size_t j = 0; j < toPrint.numCols(); j++){
      fprintf(stdout, "%s\t%s\t%lf\n", toPrint.TFLabels[i].c_str(), 
                      toPrint.GeneLabels[j].c_str(), toPrint.fullMatrix[i][j]);
    }
  }
}


void printGraph(graph<geneData, f64> *toPrint,
                                          const vector<string> &labels){
  for(size_t i = 0; i < toPrint->getNumVertexes(); i++){
    for(size_t j = 0; j < toPrint->getVertexes()[i]->getNumEdges();
                                                                  j++){
      fprintf(stderr, "%s\t%s\t%lf\n",
          labels[toPrint->getVertexes()[i]->value.nameIndex].c_str(),
labels[toPrint->getVertexes()[i]->getEdges()[j]->other(toPrint->getVertexes()[i])->value.nameIndex].c_str(),
                    toPrint->getVertexes()[i]->getEdges()[j]->weight);
    }
    fprintf(stderr, " ]\n\n"
"========================================================================"
"\n\n");
  }
}


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
