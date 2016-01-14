//
//  main.c
//  free_scale_graph
//
//  Created by xinye lei on 15/10/30.
//  Copyright © 2015年 xinye lei. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "Graph.h"
#include "PrioQ.h"
#include "Algorithm.h"
#include <pthread.h>
#define MAXIMUM 99999999

typedef enum method {
    BF,DIJ,FLY,MATRIX,PRIM,KRUSKAL
}method;

double Test(Graph *g, size_t size, method method) {
    clock_t ts, tf;
    double duration;
    int expTimes = 5;
    path **container[expTimes];
    ts = clock();
    for (int i = 0; i < expTimes; i++) {
        path **result = NULL;
        if (method == PRIM) {
            Prim(g);
        }
        if (method == KRUSKAL) {
            Kruskal(g);
        }
        if (method == BF) {
            result = AllPairDis_BF(g);
        }
        if (method == DIJ) {
            result = AllPairDis(g);
        }
        if (method == FLY) {
            result = Floyd(g);
        }
        if (method == MATRIX) {
            result = MatirxMulti(g);
        }
        container[i] = result;
    }
    tf = clock();
    duration =  (tf - ts) / (double)CLOCKS_PER_SEC / expTimes;
    for (int i = 0; i < expTimes; i++) {
        FreeMatrix(container[i], size);
    }
    return duration;
}


int main(int argc, const char * argv[]) {
    if (argc!=2)
    {
        printf("See README for usage\n");
        return 1;
    }
    const char* filename=argv[1];
    Graph *g=readGraph(filename);
    size_t size = getNumNodes(g);
    
    double t1 = Test(g, size, PRIM);
    double t2 = Test(g, size, KRUSKAL);
    double t3 = Test(g, size, DIJ);
    double t4 = Test(g, size, BF);
    double t5 = Test(g, size, FLY);
    double t6 = Test(g, size, MATRIX);
    
    printf("PRIM: %f\n",t1);
    printf("KRUSKAL: %f\n",t2);
    printf("DIJ: %f\n",t3);
    printf("BELLMAN-FORD: %f\n",t4);
    printf("FLOYD: %f\n",t5);
    printf("MATRIX MULTIPLICATION: %f\n",t6);
    freeGraph(g);

    return 0;
}
