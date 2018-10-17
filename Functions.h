#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#define K 4
#define MAX_HAMMING_SIZE 70
#define MAX_DIST 99999999
#define S1 5
#define S2 2
#define UPDATE 1

double hammingDistance(double *item1, double *item2, int dimension);
double eucledianDistance(double *item1, double *item2, int dimension);
double manhattanDistance(double *item1, double *item2, int dimension);
double cosineDistance(double *item1, double *item2, int dimension);
int* k_medoids(double** array, char distance_metric, int dimensions, int plithos_items, int K);
int* concentrate(double** array, char metriki, char type, int dimensions, int plithos_items, int K);
int* PAM_assignment(double** array, char distance_metric, int dimensions, int plithos_items, int *medoids, int K);
double* findCentroid(double** array, int plithos_items, int dimensions);
int* loydsNEW(double** array, int* assignments, int plithos_items, int dimensions, char distance_metric, int K, int* best_meds);
int *selectMoutofN(int m, int n);
int* CLARA(double** array, int n, int plithos_items, int dimensions, char metriki, char type, char method, int* best_meds, int K);
int *CLARANS(double** array, int plithos_items, int dimensions, char metriki, char type, char method, int* best_meds, int K, int S, int fraction);
double* Silhouette(double** array, int plithos_items, int dimensions, char metriki, int* medoids, int* ass, int K);
void read_configfile(char* config, int* K, int* s, int* fraction);
double findDistance(double** array, int item1, int item2, int dimensions, char metriki);
