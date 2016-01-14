//
//  Algorithm.c
//  free_scale_graph
//
//  Created by xinye lei on 15/11/14.
//  Copyright © 2015年 xinye lei. All rights reserved.
//
#include <stdlib.h>
#include "Algorithm.h"
#include "PrioQ.h"

#define MAXIMUM 99999999

/* ****************************Dijsktra************************************************************ */
path* shortpath(Graph *g, int source){
    PQueue pq = PriorityQueueCreate();    //working queue
    path *dist = malloc(sizeof(path)*getNumNodes(g));  //create result array
    for (int i=0; i<getNumNodes(g); i++) {            //initialize result array
        dist[i].mark=0;
        dist[i].index=i;
        dist[i].weight=MAXIMUM;
        dist[i].pre=-1;
        if (i==source) {
            dist[i].pre=i;
            dist[i].weight=0;
        }
    }
    PriorityQueueAddItem(pq, dist[source].weight, dist+source);
    while (PriorityQueueSize(pq)) {
        path *temp=PriorityQueueRemoveMinItem(pq);
        if (temp->mark==1) { //a node might be pushed into queue several time thus need to check mark
            continue;
        }
        temp->mark=1;        //mark the node
        
        /*
         check all its neighbors. Modify the result array if necessary. If an element(dist[i]) is modified, repush it into the queue. Since it has smaller weight than before, it will be populated ealier. This means same address (dist[i]) might be pushed into queue several times. But we always get it with smalleset weight(the latest version). And then if we meet this element again, we just ignore it.
         */
        for (unsigned neigh = neigh_first(g, temp->index); !neigh_done(g); neigh = neigh_next(g)) {
            double w=get_weight(g);
            if (dist[neigh].weight>w+dist[temp->index].weight) {
                dist[neigh].weight=w+dist[temp->index].weight;
                dist[neigh].pre=temp->index;
                PriorityQueueAddItem(pq, dist[neigh].weight, dist+neigh);
            }
        }
    }
    PriorityQueueDestroy(pq);
    return dist;
}

path** AllPairDis(Graph *g) { 
    path **result = (path**)malloc(sizeof(path*)*getNumNodes(g));
    for (int i = 0; i < getNumNodes(g); i++) {
        result[i] = shortpath(g, i);
    }
    return result;
}


/* ****************************************************************************************** */

/* ****************************Bellman Ford************************************************************ */
path* shortpath_BF(Graph *g, int source) {
    size_t size =getNumNodes(g);
    path *dist = malloc(sizeof(path)*getNumNodes(g));
    for (int i = 0; i < size; i++) {       //init distance
        dist[i].index = i;
        dist[i].weight = MAXIMUM;
        dist[i].pre = -1;
        if (i==source) {
            dist[i].pre=i;
            dist[i].weight=0;
        }
    }
    for (int i = 0; i < size; i++) {        //n times relax
        for (int j = 0; j < size; j++) {    //relax all edges
            for (unsigned neigh = neigh_first(g, j); !neigh_done(g); neigh = neigh_next(g)) {
                double weight = get_weight(g);
                if (dist[neigh].weight > weight + dist[j].weight) {
                    dist[neigh].weight = weight + dist[j].weight;
                    dist[neigh].pre = j;
                }
            }
        }
    }
    return dist;
}

path** AllPairDis_BF(Graph *g) {
    path **result = (path**)malloc(sizeof(path*)*getNumNodes(g));
    for (int i = 0; i < getNumNodes(g); i++) {
        result[i] = shortpath_BF(g, i);
    }
    return result;
}

/* ****************************************************************************************** */

/* ****************************Matrix Multiplication****************************************** */

path** MatirxMulti(Graph *g) {
    size_t size = getNumNodes(g);
    path **m1 = (path**)malloc(sizeof(path*) * size);    //alloc 2 matrices
    path **m2 = (path**)malloc(sizeof(path*) * size);
    for (int i = 0; i < size; i++) {
        m1[i] = malloc(sizeof(path) * size);
        m2[i] = malloc(sizeof(path) * size);
    }
    for (int i = 0; i < size; i++) {                     //init 2  matrices
        for (int j = 0; j < size; j++) {
            m1[i][j].weight = MAXIMUM;
            m2[i][j].weight = MAXIMUM;
            if (i == j) {
                m1[i][j].weight = 0;
                m2[i][j].weight = 0;
            }
        }
    }
    for (int i = 0; i < size; i++) {                     //fill the matrix with edge
        for (unsigned neigh = neigh_first(g, i); !neigh_done(g); neigh = neigh_next(g)) {
            m1[i][neigh].weight = get_weight(g);
            m1[neigh][i].weight = get_weight(g);
        }
    }
    
    int power = 1;                   //current power
    int flag = 1;                    //flag that detemine which matrix is orginal data
    //Here we only have 2 matrices. 1 for original data and 1 to hold the result.
    while (power < size) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double min = MAXIMUM;
                for (int k = 0; k < size; k++) {
                    if (flag == 1) {
                        min = m1[i][k].weight + m1[k][j].weight < min ? m1[i][k].weight + m1[k][j].weight : min;
                    }
                    if (flag == 2) {
                        min = m2[i][k].weight + m2[k][j].weight < min ? m2[i][k].weight + m2[k][j].weight : min;
                    }
                }
                if (flag == 1) {
                    m2[i][j].weight = min;
                }
                if (flag == 2) {
                    m1[i][j].weight = min;
                }
            }
        }
        flag = flag == 1 ? 2 : 1;    //modify flag and power after each iteration
        power *= 2;
    }
    
    //free the matrix that will not return
    if (flag == 1) {
        for (int i = 0; i < size; i++) {
            free(m1[i]);
        }
        free(m1);
        return m2;
    }
    else {
        for (int i = 0; i < size; i++) {
            free(m2[i]);
        }
        free(m2);
        return m1;
    }
}

/* ****************************Floyd************************************************************ */

path** Floyd(Graph *g) {
    path** matrix = (path**)malloc(sizeof(path*) * getNumNodes(g));     //alloc
    for (int i = 0; i < getNumNodes(g); i++) {
        matrix[i] = (path*)malloc(sizeof(path) * getNumNodes(g));
    }
    for (int i = 0; i < getNumNodes(g); i++) {                         //init
        for (int j = 0; j < getNumNodes(g); j++) {
            matrix[i][j].weight = MAXIMUM;
            if (i == j) {
                matrix[i][j].weight = 0;
            }
            matrix[i][j].nextV = -1;
        }
    }
    for (int i = 0; i < getNumNodes(g); i++) {                       //fill edges
        for (unsigned neigh = neigh_first(g, i); !neigh_done(g); neigh = neigh_next(g)) {
            matrix[i][neigh].weight = get_weight(g);
            matrix[i][neigh].nextV = neigh;
            matrix[neigh][i].weight = get_weight(g);
            matrix[neigh][i].nextV = i;
        }
    }

    for (int i = 0; i < getNumNodes(g); i++) {      //add each node, consider whether the distance decrease
        for (int j = 0; j < getNumNodes(g); j++) {
            for (int k = 0; k < getNumNodes(g); k++) {
                if (matrix[j][k].weight > matrix[j][i].weight + matrix[i][k].weight) {
                    matrix[j][k].weight = matrix[j][i].weight + matrix[i][k].weight;
                    matrix[k][j].weight = matrix[j][i].weight + matrix[i][k].weight;
                    matrix[j][k].nextV = i;
                    matrix[k][j].nextV = i;
                }
            }
        }
    }
    return matrix;
}
/* ********************************************************************************************** */

int OutputMatrix(path **matrix, int size, int cluster, char* fileName) {
    FILE *fp = fopen(fileName, "w");
    fprintf(fp, "%d %d\n",size,cluster);
    if (!fp) {
        printf("File does not exist\n");
        return 1;
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            fprintf(fp, "%f ",matrix[i][j].weight);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
    return 0;
}

int check(Graph *g) {          //check 2 algorithm's result
    path **f = MatirxMulti(g);
    path **p = AllPairDis_BF(g);
    for (int i = 0; i < getNumNodes(g); i++) {
        for (int j = 0; j < getNumNodes(g); j++) {
            if (j == i) {
                continue;
            }
            if (p[i][j].weight - f[i][j].weight > 0.000001) {
                return 1;
            }
        }    }
    for (int i = 0; i < getNumNodes(g); i++) {
        free(f[i]);
        free(p[i]);
    }
    free(f);
    free(p);
    return 0;
}

/* *****************************************Prim********************************************* */

typedef struct Mst {          //edge structure
    int from;
    int to;
    double weight;
}Mst;

double Prim(Graph *g) {
    PQueue pq = PriorityQueueCreate();   //working queue
    size_t size = getNumNodes(g);
    int mark[size];                      //mark whether an edge is existed
    for (int i = 0; i < size; i++) {     //init mark
        mark[i] = 0;
    }
    mark[0] = 1;      //init queue by pushing edge start from 0 into it
    for (unsigned neigh = neigh_first(g, 0); !neigh_done(g); neigh = neigh_next(g)) {
        Mst *mst = malloc(sizeof(Mst));
        mst->weight = get_weight(g);
        mst->from = 0;
        mst->to = neigh;
        PriorityQueueAddItem(pq, mst->weight, mst);
    }

    double mstSum = 0;
    while (PriorityQueueSize(pq)) {
        Mst *temp = PriorityQueueRemoveMinItem(pq);
        if (mark[temp->from] && mark[temp->to]) {    //edge already exist
            continue;
        }
        mstSum += temp->weight;
        int node = mark[temp->from] ? temp->to : temp->from;  //new node
        mark[node] = 1;
        for (unsigned neigh = neigh_first(g, node); !neigh_done(g); neigh = neigh_next(g)) {
            double w = get_weight(g);      //push all the edge start from new node into queue
            Mst *mst = malloc(sizeof(Mst));
            mst->from = node;
            mst->to =neigh;
            mst->weight = w;
            PriorityQueueAddItem(pq, w, mst);
        }
        free(temp);
    }
    PriorityQueueDestroy(pq);
    return mstSum;
}



/* *****************************************Kruskal********************************************* */

int *set;          //record the component a node belongs to

int find(int index) {  //find which component a node belongs to
    return index == set[index] ? index : (set[index] = find(set[index]));
}

double Kruskal(Graph *g) {
    PQueue pq = PriorityQueueCreate();
    size_t size = getNumNodes(g);
    set = malloc(sizeof(int) * size);
    for (int i = 0; i < size; i++) {   //init array hold
        set[i] = i;
    }
    for (int i = 0; i < size; i++) {
        for (unsigned neigh = neigh_first(g, i); !neigh_done(g); neigh = neigh_next(g)) {
            double w = get_weight(g);
            Mst *temp = malloc(sizeof(Mst));
            temp->from = i;
            temp->to = neigh;
            temp->weight = w;
            PriorityQueueAddItem(pq, w, temp);
        }
    }
    double mstSum = 0;
    while (PriorityQueueSize(pq)) {
        Mst *mst = PriorityQueueRemoveMinItem(pq);
        if (find(mst->from) != find(mst->to)) {
            mstSum += mst->weight;
            set[find(mst->from)] = find(mst->to);
        }
        free(mst);
    }
    free(set);
    PriorityQueueDestroy(pq);
    return mstSum;
}

int FreeMatrix(path **matrix, size_t size) {
    if (!matrix) {
        return 1;
    }
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
    return 0;
}











