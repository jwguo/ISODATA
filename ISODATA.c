/* 9913525
 *
 * ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

// ==========need to change to command line options============================
#define DESIRED_CLUSTER_NUM 6
#define CLUSTER_SIZE_MIN 3 // m0
#define STD_DEVIATION_THRESH 6 // for splitting
#define SPLIT_FRACTION 0.8 // belong to (0, -1]
#define LUMP_THRESH 12 // for lumping unit: distance
#define LUMP_PAIR_MAX_PER_IT 2
#define ITERATION_MAX 20
#define DATA_PATH "./data.txt"
// ============================================================================

#define LINE_BUF_SIZE 128
#define TOKEN_BUF_SIZE 32
#define INIT_FEATURE_VECTOR_NUM 16

typedef struct
{
    double w;
    double x;
    double y;
    double z;
    int    cluster;
    int    core;
} FVector, Center;

void readData(FILE *fp, Center **pCenters, int *centerNdx, FVector **pFVectors, int *fVectorNdx);
void _readFVector(FVector **pFVector, int *index, int *total, char *token);
void didayDynamicClusterMethod(FVector *pFVectors,
                               int      fVectorTotal,
                               int      clusterCnt,
                               int     *pCoreCnt,
                               int      iteration);
double twoNorm(FVector a, FVector b);
double distCE(FVector  *pFVectors,
              int       fVectorTotal,
              int       clusterCnt,
              FVector **pCores,
              int      *coresNdx);
void randCore(FVector  *pFVectors,
              int       fVectorTotal,
              int       clusterCnt,
              int      *pCoreCnt,
              int     **pCluster,
              int      *clusterNdx,
              FVector **pCores,
              int      *coresNdx);

int main(int argc, char **argv)
{
    int  i;
    FILE *fp;

    fp = fopen(DATA_PATH, "r");
    if (NULL == fp) {
        fprintf(stderr, "error opening file %s\n", DATA_PATH);
        exit(EXIT_FAILURE);
    }

    Center *pCenters = calloc(INIT_FEATURE_VECTOR_NUM, sizeof(Center));
    int centerNdx = 0;
    FVector *pFVectors = calloc(INIT_FEATURE_VECTOR_NUM, sizeof(FVector));
    int fVectorNdx = 0;
    readData(fp, &pCenters, &centerNdx, &pFVectors, &fVectorNdx);

    /*int coreCnt[DESIRED_CLUSTER_NUM];
    for (i = 0; i < DESIRED_CLUSTER_NUM; i++) {
        if (0 == i)
            coreCnt[i] = N0;
        if (1 == i)
            coreCnt[i] = N1;
        if (2 == i)
            coreCnt[i] = N2;
        if (3 == i)
            coreCnt[i] = N3;
    }
    didayDynamicClusterMethod(pFVectors, fVectorNdx, DESIRED_CLUSTER_NUM, coreCnt, ITERATION_MAX);*/

    return 0;
}

void readData(FILE *fp, Center **pCenters, int *centerNdx, FVector **pFVectors, int *fVectorNdx)
{
    char  lineBuf[LINE_BUF_SIZE];
    char *token;
    char *delim = ",";
    int   centerTotal  = INIT_FEATURE_VECTOR_NUM;
    int   fVectorTotal = INIT_FEATURE_VECTOR_NUM;
    int   flag = 0;

    while (fgets(lineBuf, LINE_BUF_SIZE, fp)) {
        if (0 == strcmp(lineBuf, "\n")) {
            // change to read sample feature vectors
            flag = 1;
            continue;
        }
        token = strtok(lineBuf, delim);
        if (0 == flag) {
            // read initial cluster centers
            _readFVector(*pCenters, *centerNdx, &centerTotal, token);
        } else if (1 == flag) {
            // read sample feature vectors
            _readFVector(*pFVectors, *fVectorNdx, &fVectorTotal, token);
        }
    }
    int i;
    for (i = 0; i < *fVectorNdx; i++) {
        printf("(%d, %d, %d, %d)\n", (int)(*pFVectors)[i].w,
                                     (int)(*pFVectors)[i].x,
                                     (int)(*pFVectors)[i].y,
                                     (int)(*pFVectors)[i].z);
    }
}

void _readFVector(FVector **pFVector, int *index, int *total, char *token)
{
    char *delim = ",";
    char  compo = 'w';

    while (NULL != token) {
        if ('w' == compo) {
            compo = 'x';
            (*pFVector)[*index].w = atof(token);
        } else if ('x' == compo) {
            compo = 'y';
            (*pFVector)[*index].x = atof(token);
        } else if ('y' == compo) {
            compo = 'z';
            (*pFVector)[*index].y = atof(token);
        } else if ('z' == compo) {
            compo = 'w';
            (*pFVector)[*index].cluster = -1;
            (*pFVector)[*index].core  = -1;
            (*pFVector)[(*index)++].z = atof(token);
            if (*total <= *index) {
                (*total) += 2 * INIT_FEATURE_VECTOR_NUM;
                *pFVector = realloc(*pFVector, sizeof(FVector) * (*total));
            }
        }
        token = strtok(NULL, delim);
    }
}

void didayDynamicClusterMethod(FVector *pFVectors,
                               int      fVectorTotal,
                               int      clusterCnt,
                               int     *pCoreCnt,
                               int      iteration)
{
    srand(time(NULL));
    int i, j, k, m;
    FVector  initCenters[clusterCnt];
    FVector *pCores[clusterCnt];
    int      coresNdx[clusterCnt];
    int     *pCluster[clusterCnt];
    int      clusterNdx[clusterCnt];
    double   minDist;
    double   tmp;
    FVector *pCandidate;

    for (i = 0; i < clusterCnt; i++) {
        pCores[i]     = calloc(pCoreCnt[i], sizeof(FVector));
        pCluster[i]   = calloc(fVectorTotal, sizeof(int));
        coresNdx[i]   = 0;
        clusterNdx[i] = 0;
    }

    // assign the first clusterCnt feature vectors as initCenters
    for (i = 0; i < clusterCnt; i++) {
        initCenters[i] = pFVectors[i];
        initCenters[i].cluster = i;
    }

    // assign feature vectors to clusters by min. distance
    for (i = 0; i < fVectorTotal; i++) {
        minDist = DBL_MAX;
        for (j = 0; j < clusterCnt; j++) {
            tmp = twoNorm(initCenters[j], pFVectors[i]);
            if (tmp < minDist) {
                pCandidate = &(initCenters[j]);
                minDist = tmp;
            }
        }
        pFVectors[i].cluster = pCandidate->cluster;
        pCluster[pCandidate->cluster][clusterNdx[pCandidate->cluster]++] = i;
    }

    // randomly select multicenter cores for all clusters
    randCore(pFVectors,
             fVectorTotal,
             clusterCnt,
             pCoreCnt,
             pCluster,
             clusterNdx,
             pCores,
             coresNdx);

    // find the best cluster sets and core sets iteratively
    int bestCoreSet[fVectorTotal];
    int bestClusterSet[fVectorTotal];
    double minDistCE = DBL_MAX;
    for (i = 0; i < iteration; i++) {
        for (j = 0; j < clusterCnt; j++)
            clusterNdx[j] = 0;

        for (j = 0; j < fVectorTotal; j++) {
            // assign feature vectors to clusters via min. dist. to cores
            minDist = DBL_MAX;
            for (k = 0; k < clusterCnt; k++) {
                for (m = 0; m < coresNdx[k]; m++) {
                    tmp = twoNorm(pCores[k][m], pFVectors[j]);
                    if (tmp < minDist) {
                        pCandidate = &(pCores[k][m]);
                        minDist = tmp;
                    }
                }
            }
            pFVectors[j].cluster = pCandidate->cluster;
            pCluster[pCandidate->cluster][clusterNdx[pCandidate->cluster]++] = j;
        }
        // compute D(k) and D(C, E)
        tmp = distCE(pFVectors, fVectorTotal, clusterCnt, pCores, coresNdx);
        if (tmp < minDistCE) {
            // candidate
            minDistCE = tmp;
            for (j = 0; j < fVectorTotal; j++) {
                bestCoreSet[j] = pFVectors[j].core;
                bestClusterSet[j] = pFVectors[j].cluster;
            }
        }
        // reselect random new core sets
        randCore(pFVectors,
                 fVectorTotal,
                 clusterCnt,
                 pCoreCnt,
                 pCluster,
                 clusterNdx,
                 pCores,
                 coresNdx);
    }

    // print final results to standard output
    int memberCnts[clusterCnt];
    for (i = 0; i < clusterCnt; i++)
        memberCnts[i] = 0;
    for (i = 0; i < clusterCnt; i++) {
        printf("\ncluster #%d:\n", i);
        for (j = 0; j < fVectorTotal; j++) {
            if (pFVectors[j].cluster == i) {
                memberCnts[i]++;
                printf("(%d, %d, %d)",
                       (int)pFVectors[j].x,
                       (int)pFVectors[j].y,
                       (int)pFVectors[j].z);
                if (pFVectors[j].core == i)
                    printf(" <- core of cluster #%d", i);
                printf("\n");
            }
        }
    }
    printf("\n");
    for (i = 0; i < clusterCnt; i++)
        printf("# of cluster %d: %d\n", i, memberCnts[i]);

    printf("\nmin( D(C,E) ): %f\n", minDistCE);
}

double twoNorm(FVector a, FVector b)
{
    return (pow(pow(a.x - b.x, 2) +
                pow(a.y - b.y, 2) +
                pow(a.z - b.z, 2), 0.5));
}

double distCE(FVector  *pFVectors,
              int       fVectorTotal,
              int       clusterCnt,
              FVector **pCores,
              int      *coresNdx)
{
    int i, j;
    int clusterNum;
    double result = 0;

    for (i = 0; i < fVectorTotal; i++) {
        // traverse feature vectors
        clusterNum = pFVectors[i].cluster;
        for (j = 0; j < coresNdx[clusterNum]; j++) {
            result += twoNorm(pCores[clusterNum][j], pFVectors[i]);
        }
    }

    return result;
}

void randCore(FVector  *pFVectors,
              int       fVectorTotal,
              int       clusterCnt,
              int      *pCoreCnt,
              int     **pCluster,
              int      *clusterNdx,
              FVector **pCores,
              int      *coresNdx)
{
    int i, j;
    int randNdx;

    for (i = 0; i < clusterCnt; i++)
        coresNdx[i] = 0;

    for (i = 0; i < fVectorTotal; i++)
        pFVectors[i].core = -1;

    for (i = 0; i < clusterCnt; i++) {
        // traverse each cluster
        if (clusterNdx[i] <= pCoreCnt[i]) {
            // cluster member count is smaller than default core count
            for (j = 0; j < clusterNdx[i]; j++) {
                // each cluster member is core
                pFVectors[pCluster[i][j]].core = i;
                pCores[i][coresNdx[i]++] = pFVectors[pCluster[i][j]];
            }
        } else {
            // randomly select cores
            for (j = 0; j < pCoreCnt[i]; ) {
                randNdx = rand() % clusterNdx[i];
                if (pFVectors[pCluster[i][randNdx]].core != i) {
                    pFVectors[pCluster[i][randNdx]].core = i;
                    pCores[i][coresNdx[i]++] = pFVectors[pCluster[i][randNdx]];
                    j++;
                }
            }
        }
    }
}

