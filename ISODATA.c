/* 9913525. Jing-Wei, Guo
 *
 * The dimension of feature vectors should be made programmable
 *
 * ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

// ===========need to change to command line options============================
#define DIMENSION 4
#define DESIRED_CLUSTER_NUM 6
#define CLUSTER_SIZE_MIN 3 // m0
#define STD_DEVIATION_THRESH 6 // for splitting
#define SPLIT_FRACTION 0.8 // belong to (0, -1]
#define LUMP_THRESH 12 // for lumping unit: distance
#define LUMP_PAIR_MAX_PER_IT 2
#define ITERATION_MAX 20
#define EPSILON 0.001 // threshold for determining changeFlag
#define DATA_PATH "./data.txt"
// =============================================================================

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
    int             totalCenters;
} FVector, Center;

void readData(FILE     *fp,
              Center  **pCenters,
              int      *centerNdx,
              FVector **pFVectors,
              int      *fVectorNdx);
void _readFVector(FVector **pFVector, int *index, int *total, char *token);
int // 0: end looping, 1: keep looping
isodata(Center  *pCenters,
        FVector *pFVectors,
        int      totalFVectors,
        int      desiredClusterNum,
        int      clusterSizeMin,
        double   stdDeviationThresh,
        double   splitFraction,
        double   lumpThresh,
        int      lumpPairMaxPerIt,
        int      iterationMax,
        double   epsilon,
        int     *pSplitFlags,
        int     *pLumpFlags);
void classify(FVector **pFVectors,
              int       totalFVectors,
              Center   *pCenters);
int // 0: next iteration, 1: do lumping, 2: end
_split(int     changeFlag,
       int    *pSplitFlag,
       int    *pLastLumpFlag,
       Center *pCenters,
       int    *realClusterNum,
       double  stdDeviationThresh,
       int     desiredClusterNum,
       double  avgDistToCenters,
       double *pAvgDistToCenter,
       int     clusterSizeMin,
       int     loopCnt,
       double  splitFraction);
void __split(int      index,
             double   sigma,
             FVector *pSplitCenter,
             FVector *pCenters,
             double   splitFraction);
void   _lump(int *pLumpFlag, Center *pCenter, int realClusterNum);
double twoNorm(FVector a, FVector b);

int main(int argc, char **argv)
{
    int   i;
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

    int splitFlag[ITERATION_MAX];
    int lumpFlag[ITERATION_MAX];
    /* flag meaning: 0: splitting or lumping starts
     *               1: splitting or lumping is completed successfully
     *               2: flag's initial state after the i'th iteration
     * ************************************************************************/
    // step 1: initialize {split,lump}Flag
    for (i = 0; i < ITERATION_MAX; i++)
        splitFlag[i] = lumpFlag[i] = 2;
    int loopFlag = 1;
    while (1 == loopFlag)
        loopFlag = isodata(pCenters,
                           pFVectors,
                           fVectorNdx,
                           DESIRED_CLUSTER_NUM,
                           CLUSTER_SIZE_MIN,
                           STD_DEVIATION_THRESH,
                           SPLIT_FRACTION,
                           LUMP_THRESH,
                           LUMP_PAIR_MAX_PER_IT,
                           ITERATION_MAX,
                           EPSILON,
                           splitFlag,
                           lumpFlag);
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
        if (0 == flag) {
            // read initial cluster centers
            _readFVector(&*pCenters, &*centerNdx, &centerTotal, token);
        }
        else if (1 == flag)
            // read sample feature vectors
            _readFVector(&*pFVectors, &*fVectorNdx, &fVectorTotal, token);
    }
    int i;
    for (i = 0; i < *centerNdx - 1; i++)
        (*pCenters)[i].pNextCenter = &((*pCenters)[i + 1]);
    (*pCenters)[0].totalCenters = *centerNdx;

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
            (*pFVector)[*index].cluster      = -1;
            (*pFVector)[*index].core         = -1;
            (*pFVector)[*index].pNextMember  = NULL;
            (*pFVector)[*index].pNextCenter  = NULL;
            (*pFVector)[*index].totalMembers = 0;
            (*pFVector)[*index].totalCenters = 0;
            (*pFVector)[(*index)++].z        = atof(token);
            if (*total <= *index) {
                (*total) += 2 * INIT_FEATURE_VECTOR_NUM;
                *pFVector = realloc(*pFVector, sizeof(FVector) * (*total));
            }
        }
        token = strtok(NULL, delim);
    }
}

int isodata(Center  *pCenters,
        FVector *pFVectors,
        int      totalFVectors,
        int      desiredClusterNum,
        int      clusterSizeMin,
        double   stdDeviationThresh,
        double   splitFraction,
        double   lumpThresh,
        int      lumpPairMaxPerIt,
        int      iterationMax,
        double   epsilon,
        int     *pSplitFlags,
        int     *pLumpFlags)
{
    static int loopCnt = 0;
    int      i, j, rc;
    int      changeFlag;
    int      clusterCntChangeFlag = 0;
    int      realClusterNum;
    //Center  *pNewCenters  = calloc(centerNdx, sizeof(Center));
    //int      newCenterNdx = centerNdx;
    FVector *pMember;
    int      memberCnt;
    Center  *pCenter;

    //realClusterNum = centerNdx;
    //pNewCenters[0] = pCenters[0];
    //for (i = 1; i < newCenterNdx; i++) {
    //    pNewCenters[i] = pCenters[i];
    //    pNewCenters[i-1].pNextCenter = &(pNewCenters[i]); // !! free
    //}

    // step 2: classify sample feature vectors
    changeFlag = 1;
    FVector* classifiedFVectors[totalFVectors];
    for (i = 0; i < totalFVectors; i++)
        classifiedFVectors[i] = &(pFVectors[i]);
    classify(classifiedFVectors, totalFVectors, pCenters);

    // step 3: remove clusters that have not enough members
    realClusterNum = pCenters[0].totalCenters;
    pCenter = &(pCenters[0]);
    Center *pLastCenter = NULL;
    while (NULL != pCenter) {
        if (pCenter->totalMembers < clusterSizeMin) {
            // members not enough, reclassify members in  this cluster
            clusterCntChangeFlag = 1;
            memberCnt = pCenter->totalMembers;
            realClusterNum--;

            if (NULL == pLastCenter)
                *pCenter = *(pCenter->pNextCenter);
            else
                pLastCenter->pNextCenter = pCenter->pNextCenter;

            FVector* pReclassifiedFVectors[memberCnt];
            pMember = pCenter->pNextMember;
            for (j = 0; j < memberCnt; j++) {
                pReclassifiedFVectors[j] = pMember;
                pMember = pMember->pNextMember;
            }
            classify(pReclassifiedFVectors, memberCnt, pCenters);
        }
        pLastCenter = pCenter;
        pCenter     = pCenter->pNextCenter;
    }
    pCenters[0].totalCenters = realClusterNum;

    // step 4: calculate new centers and decide changeFlag
    pCenter = &(pCenters[0]);
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

        pCenter->w = newCenter[i].w;
        pCenter->x = newCenter[i].x;
        pCenter->y = newCenter[i].y;
        pCenter->z = newCenter[i].z;

        pCenter = pCenter->pNextCenter;
    }
    if (0 == clusterCntChangeFlag && totalDiff < epsilon)
        changeFlag = 0;

    // step 5: calculate average distance to centers
    double avgDistToCenters = 0;
    double avgDistToCenter[realClusterNum];
    double dist;
    for (i = 0; i < realClusterNum; i++)
        avgDistToCenter[i] = 0;

    pCenter = &(pCenters[0]);
    for (i = 0; i < realClusterNum; i++) {
        pMember = pCenter->pNextMember;
        while (NULL != pMember) {
            dist = twoNorm(*pMember, *pCenter);
            avgDistToCenter[i] += dist;
            avgDistToCenters   += dist;
            pMember = pMember->pNextMember;
        }
        avgDistToCenter[i] /= pCenter->totalMembers;
        pCenter = pCenter->pNextCenter;
    }
    avgDistToCenters /= totalFVectors;

    // step 6:
    loopCnt++;
    if (loopCnt >= iterationMax) // end
        return 0;
    else if (realClusterNum <= (desiredClusterNum + 1) / 2 ||
             (realClusterNum < 2 * desiredClusterNum &&
              1 == loopCnt % 2)) {
        // split condition fulfilled
        rc = _split(changeFlag,
                    &(pSplitFlags[loopCnt]),
                    &(pLumpFlags[loopCnt - 1]),
                    pCenters,
                    &realClusterNum,
                    stdDeviationThresh,
                    desiredClusterNum,
                    avgDistToCenters,
                    avgDistToCenter,
                    clusterSizeMin,
                    loopCnt,
                    splitFraction);
        if (0 == rc) { // next iteration
            return 1;
        } else if (1 == rc) { // do lumping
            //rc = _lump(&(pSplitFlags[loopCnt]), pNewCenters, realClusterNum);
            if (0 == rc) { // next iteration
                return 1;
            } else // end
                return 0;
        } else // end
            return 0;
    } else {
        // lump condition fulfilled
        //rc = _lump(&(pSplitFlags[loopCnt]), pNewCenters, realClusterNum);
        if (0 == rc) { // next iteration
            return 1;
        } else // end
            return 0;
    }
}

void classify(FVector **pFVectors,
              int       totalFVectors,
              Center   *pCenters)
{
    int      i;
    double   tmpDouble;
    double   minDist;
    FVector *pCandidate;
    Center  *pCenter;

    for (i = 0; i < totalFVectors; i++) {
        minDist = DBL_MAX;
        pCenter = &(pCenters[0]);
        while (NULL != pCenter) {
            tmpDouble = twoNorm(*(pFVectors[i]), *pCenter);
            if (tmpDouble > minDist) {
                minDist    = tmpDouble;
                pCandidate = pCenter;
            }
            pCenter = pCenter->pNextCenter;
        }
        pFVectors[i]->pNextMember = pCandidate->pNextMember;
        pCandidate->pNextMember   = pFVectors[i];
        (pCandidate->totalMembers)++;
    }
}

int _split(int     changeFlag,
           int    *pSplitFlag,
           int    *pLastLumpFlag,
           Center *pCenters,
           int    *realClusterNum,
           double  stdDeviationThresh,
           int     desiredClusterNum,
           double  avgDistToCenters,
           double *pAvgDistToCenter,
           int     clusterSizeMin,
           int     loopCnt,
           double  splitFraction)
{
    // step 7:
    printf("SPLIT\n");
    int      i, j;
    double   sigma[*realClusterNum][DIMENSION];
    Center  *pCenter;
    FVector *pMember;

    for (i = 0; i < *realClusterNum; i++)
        for (j = 0; j < DIMENSION; j++)
            sigma[i][j] = 0;

    *pSplitFlag = 0;
    pCenter = &(pCenters[0]);
    for (i = 0; i < *realClusterNum; i++) {
        pMember = pCenter->pNextMember;
        while (NULL != pMember) {
            for (j = 0; j < DIMENSION; j++) {
                if (0 == j)
                    sigma[i][j] += pow(pMember->w - pCenter->w, 2); // !!
                else if (1 == j)
                    sigma[i][j] += pow(pMember->x - pCenter->x, 2);
                else if (2 == j)
                    sigma[i][j] += pow(pMember->y - pCenter->y, 2);
                else if (3 == j)
                    sigma[i][j] += pow(pMember->z - pCenter->z, 2);
            }
            pMember = pMember->pNextMember;
        }
        for (j = 0; j < DIMENSION; j++)
            sigma[i][j] = pow(sigma[i][j] / pCenter->totalMembers, 0.5);
        pCenter = pCenter->pNextCenter;
    }

    // step 8
    pCenter = &(pCenters[0]);
    for (i = 0; i < *realClusterNum; i++) {
        for (j = 0; j < DIMENSION; j++)
            if (sigma[i][j] > stdDeviationThresh &&
                (*realClusterNum <= (desiredClusterNum + 1) / 2 ||
                 pAvgDistToCenter[i] > avgDistToCenters ||
                 pCenter->totalMembers >= 2 * clusterSizeMin)) {
                // split
                __split(j, sigma[i][j], pCenter, pCenters, splitFraction);
                *pSplitFlag = 1;
                (*realClusterNum)++;
                break;
            }
        pCenter = pCenter->pNextCenter;
    }

    if (1 == *pSplitFlag)
        return 0; // next iteration
    else if (1 >= loopCnt)
        return 1; // do lumping
    else if (0 != *pLastLumpFlag)
        return 1; // do lumping
    else if (1 == changeFlag)
        return 0; // next iteration
    else
        return 2; // end
}

void __split(int      index,
             double   sigma,
             FVector *pSplitCenter,
             FVector *pCenters,
             double   splitFraction)
{
    int     i;
    int     oldCenterCnt = pCenters[0].totalCenters;
    Center *pCenter;
    Center *pLastCenter;
    Center *pNewCenters = calloc(++(pCenters[0].totalCenters), sizeof(Center));

    pCenter     = &(pCenters[0]);
    pLastCenter = NULL;
    for (i = 0; i < oldCenterCnt; i++) {
        if (pCenter != pSplitCenter)
            pNewCenters[i] = *pCenter;
        else if (NULL != pLastCenter)
            pLastCenter->pNextCenter = pCenter->pNextCenter;

        pLastCenter = pCenter;
        pCenter     = pCenter->pNextCenter;
    }

    pSplitCenter->pNextCenter  = NULL;
    pSplitCenter->pNextMember  = NULL;
    pSplitCenter->totalMembers = 0;

    Center *pSplit1 = &(pNewCenters[oldCenterCnt - 1]);
    Center *pSplit2 = &(pNewCenters[oldCenterCnt]);

    pNewCenters[oldCenterCnt - 2].pNextCenter = pSplit1;
    *pSplit1 = *pSplit2 = *pSplitCenter;
    pSplit1->pNextCenter = pSplit2;

    if (0 == index) {
        pSplit1->w += splitFraction * sigma;
        pSplit2->w -= splitFraction * sigma;
    } else if (1 == index) {
        pSplit1->x += splitFraction * sigma;
        pSplit2->x -= splitFraction * sigma;
    } else if (2 == index) {
        pSplit1->y += splitFraction * sigma;
        pSplit2->y -= splitFraction * sigma;
    } else if (3 == index) {
        pSplit1->z += splitFraction * sigma;
        pSplit2->z -= splitFraction * sigma;
    }

    free(pCenters);
}

void _lump(int *pLumpFlag, Center *pCenter, int realClusterNum)
{
    printf("LUMP\n");
}

double twoNorm(FVector a, FVector b)
{
    return (pow(pow(a.w - b.w, 2) +
                pow(a.x - b.x, 2) +
                pow(a.y - b.y, 2) +
                pow(a.z - b.z, 2), 0.5));
}

