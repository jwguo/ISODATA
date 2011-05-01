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
#define EPSILON 0.001 // threshold for determining changeFlag
#define DATA_PATH "./data.txt"
// ============================================================================

#define LINE_BUF_SIZE 128
#define TOKEN_BUF_SIZE 32
#define INIT_FEATURE_VECTOR_NUM 16

typedef struct FVector
{
    double          w;
    double          x;
    double          y;
    double          z;
    int             cluster;
    int             core;
    struct FVector *pNextMember;
    struct FVector *pNextCenter;
    int             totalMembers;
} FVector, Center;

void readData(FILE     *fp,
              Center  **pCenters,
              int      *centerNdx,
              FVector **pFVectors,
              int      *fVectorNdx);
void _readFVector(FVector **pFVector, int *index, int *total, char *token);
void isodata(Center  *pInitCenters,
             int      initCenterNdx,
             FVector *pFVectors,
             int      fVectorNdx,
             int      desiredClusterNum,
             int      clusterSizeMin,
             double   stdDeviationThresh,
             double   splitFraction,
             double   lumpThresh,
             int      lumpPairMaxPerIt,
             int      iterationMax,
             double   epsilon);
void classify(FVector **pFVectors,
              int       fVectorNdx,
              Center   *pNewCenters);
//void didayDynamicClusterMethod(FVector *pFVectors,
//                               int      fVectorTotal,
//                               int      clusterCnt,
//                               int     *pCoreCnt,
//                               int      iteration);
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

    isodata(pCenters,
            centerNdx,
            pFVectors,
            fVectorNdx,
            DESIRED_CLUSTER_NUM,
            CLUSTER_SIZE_MIN,
            STD_DEVIATION_THRESH,
            SPLIT_FRACTION,
            LUMP_THRESH,
            LUMP_PAIR_MAX_PER_IT,
            ITERATION_MAX,
            EPSILON);

    return 0;
}

void readData(FILE     *fp,
              Center  **pCenters,
              int      *centerNdx,
              FVector **pFVectors,
              int      *fVectorNdx)
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
        if (0 == flag)
            // read initial cluster centers
            _readFVector(&*pCenters, &*centerNdx, &centerTotal, token);
        else if (1 == flag)
            // read sample feature vectors
            _readFVector(&*pFVectors, &*fVectorNdx, &fVectorTotal, token);
    }
    int i;
    printf("Centers:\n");
    for (i = 0; i < *centerNdx; i++) {
        printf("(%d, %d, %d, %d)\n", (int)(*pCenters)[i].w,
                                     (int)(*pCenters)[i].x,
                                     (int)(*pCenters)[i].y,
                                     (int)(*pCenters)[i].z);
    }
    printf("\nSample vectors:\n");
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
            (*pFVector)[*index].cluster     = -1;
            (*pFVector)[*index].core        = -1;
            (*pFVector)[*index].pNextMember = NULL;
            (*pFVector)[*index].pNextCenter = NULL;
            (*pFVector)[*index].totalMembers       = 0;
            (*pFVector)[(*index)++].z       = atof(token);
            if (*total <= *index) {
                (*total) += 2 * INIT_FEATURE_VECTOR_NUM;
                *pFVector = realloc(*pFVector, sizeof(FVector) * (*total));
            }
        }
        token = strtok(NULL, delim);
    }
}

void isodata(Center  *pInitCenters,
             int      initCenterNdx,
             FVector *pFVectors,
             int      fVectorNdx,
             int      desiredClusterNum,
             int      clusterSizeMin,
             double   stdDeviationThresh,
             double   splitFraction,
             double   lumpThresh,
             int      lumpPairMaxPerIt,
             int      iterationMax,
             double   epsilon)
{
    int      i, j;
    int      changeFlag;
    int      clusterCntChangeFlag = 0;
    int      realClusterNum;
    int      splitFlag[iterationMax];
    int      lumpFlag[iterationMax];
    Center  *pNewCenters  = calloc(initCenterNdx, sizeof(Center));
    int      newCenterNdx = initCenterNdx;
    FVector *pMember;
    int      memberCnt;

    // step 1: initialize {split,lump}Flag
    for (i = 0; i < iterationMax; i++)
        splitFlag[i] = lumpFlag[i] = 2;

    // step 2: apply new centers, and classify sample feature vectors
    realClusterNum = initCenterNdx;
    changeFlag     = 1;
    pNewCenters[0] = pInitCenters[0];
    for (i = 1; i < newCenterNdx; i++) {
        pNewCenters[i] = pInitCenters[i];
        pNewCenters[i-1].pNextCenter = &(pNewCenters[i]);
    }

    FVector* classifiedFVectors[fVectorNdx];
    for (i = 0; i < fVectorNdx; i++)
        classifiedFVectors[i] = &(pFVectors[i]);
    classify(classifiedFVectors, fVectorNdx, pNewCenters);

    // step 3: remove clusters that have not enough members
    for (i = 0; i < newCenterNdx; i++) {
        if (pNewCenters[i].totalMembers < clusterSizeMin) {
            // members not enough, reclassify members in  this cluster
            clusterCntChangeFlag = 1;
            memberCnt = pNewCenters[i].totalMembers;
            realClusterNum--;
            if (0 == i) {
                pNewCenters[i] = pNewCenters[i + 1];

            } else
                pNewCenters[i - 1].pNextCenter = pNewCenters[i].pNextCenter;

            FVector* pReclassifiedFVectors[memberCnt];
            pMember = pNewCenters[i].pNextMember;
            for (j = 0; j < memberCnt; j++) {
                pReclassifiedFVectors[j] = pMember;
                pMember = pMember->pNextMember;
            }
            classify(pReclassifiedFVectors, memberCnt, pNewCenters);
        }
    }

    // step 4: calculate new centers and decide changeFlag
    Center *pCenter = &(pNewCenters[0]);
    Center newCenter[realClusterNum];
    double totalDiff = 0;
    for (i = 0; i < realClusterNum; i++) {
        pMember = pCenter->pNextMember;
        while (NULL != pMember) {
            newCenter[i].w += pMember->w;
            newCenter[i].x += pMember->x;
            newCenter[i].y += pMember->y;
            newCenter[i].z += pMember->z;

            pMember = pMember->pNextMember;
        }
        memberCnt = pCenter->totalMembers;
        newCenter[i].w /= memberCnt;
        newCenter[i].x /= memberCnt;
        newCenter[i].y /= memberCnt;
        newCenter[i].z /= memberCnt;

        totalDiff += twoNorm(newCenter[i], *pCenter);
        pCenter = pCenter->pNextCenter;
    }
    if (0 == clusterCntChangeFlag && totalDiff < epsilon)
        changeFlag = 0;

    // step 5:
}

void classify(FVector **pFVectors,
              int       fVectorNdx,
              Center   *pNewCenters)
{
    int      i;
    double   tmpDouble;
    double   minDist;
    FVector *pCandidate;
    FVector *child;
    Center  *pCenter;

    for (i = 0; i < fVectorNdx; i++) {
        minDist = MAXFLOAT;
        pCenter = &(pNewCenters[0]);
        while (NULL != pCenter) {
            tmpDouble = twoNorm(*(pFVectors[i]), *pCenter);
            if (tmpDouble > minDist) {
                minDist      = tmpDouble;
                pCandidate   = pCenter;
            }
            pCenter = pCenter->pNextCenter;
        }
        child = pCandidate->pNextMember;
        pCandidate->pNextMember = pFVectors[i];
        pFVectors[i]->pNextMember = child;
        (pCandidate->totalMembers)++;
    }
}

/*void didayDynamicClusterMethod(FVector *pFVectors,
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
}*/

double twoNorm(FVector a, FVector b)
{
    return (pow(pow(a.w - b.w, 2) +
                pow(a.x - b.x, 2) +
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

