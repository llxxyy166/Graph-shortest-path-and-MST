//
//  Algorithm.h
//  free_scale_graph
//
//  Created by xinye lei on 15/11/14.
//  Copyright © 2015年 xinye lei. All rights reserved.
//

#ifndef Algorithm_h
#define Algorithm_h

#include <stdio.h>
#include "Graph.h"


typedef struct path{  //structure to store the path
    int pre;          //previous node in the shortest path, used by dij algorithm
    int index;        //node ID, used by all algorithm
    double weight;    //distance, used by all algorithm
    int mark;         //mark that wether this node is processed, used by dij algorithm
    int nextV;        //next vertex on the path, used by all algorithm
}path;


path* shortpath(Graph *g, int source); //Dijkstra algorithm
path** AllPairDis(Graph *g);           //Dij algorithm find all pair distance
/*
 * Heap based Dijkstra algorithm. First function calculate from single source,
 * Second call the first n times to compute all pair distance.
 * Only used for graph without negative weight.
 */


path* shortpath_BF(Graph *g, int source);
path** AllPairDis_BF(Graph *g);
/*
 * Bellman-Ford algorithm. First compute shorteset distance from single source,
 * second call the first function compute all pairs. Can be used for graph with
 * negative weight. 
 */


path** Floyd(Graph *g);
/*
 * Floyd algorithm directly compute all pair distance.Can find the path by visit
 * nextV in each node.
 */


path** MatirxMulti(Graph *g);
/*
 * Matrix multiplicationalgorithm directly compute all pair distance. Only get the
 * distance, can not find the path.
 */


int OutputMatrix(path **matrix, int size, int cluster, char *fileName);
/*
 * Print the matrix into a text file so that partner can analyse.
 */

int check(Graph *g); //check whether 2 algorithms provide same result.


double Prim(Graph *g);
double Kruskal(Graph *g);
/*
 * MST algorithm, find the 
 */

int FreeMatrix(path **matrix, size_t size);

#endif /* Algorithm_h */
