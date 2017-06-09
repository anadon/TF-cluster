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


void printCoincidenceMatrix(UpperDiagonalSquareMatrix<u8> *matrix, 
                                cu8 maxMatch, const vector<string> TFs){
  f64 **mtr;
  size_t n = matrix->getSideLength();
  mtr = (f64**) malloc(sizeof(*mtr) * n);
  for(size_t i = 0; i < n; i++)
    mtr[i] = (f64*) malloc(sizeof(**mtr) * n);

  for(size_t i = 0; i < n; i++){\
    mtr[i][i] = maxMatch;
    for(size_t j = i+1; j < n; j++){
      mtr[i][j] = mtr[j][i] = matrix->getValueAtIndex(i, j);
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
