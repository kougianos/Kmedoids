#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Functions.h"
//#define K 4
#define MAX_HAMMING_SIZE 70
#define MAX_DIST 99999999
#define S1 5
#define S2 2
#define UPDATE 1

int* k_medoids(double** array, char distance_metric, int dimensions, int plithos_items, int K)
{
	int i, *kmeds, thesi_k, j, int_sd;
	double sum_dist, d_r, probability, dist, *min_distances;
	
	kmeds = (int*)malloc(K*sizeof(int));
	min_distances = (double*)malloc(plithos_items*sizeof(double));
	
	kmeds[0] = rand()%plithos_items;
	thesi_k=1;
	
	for(thesi_k=1; thesi_k<K; thesi_k++)
	{
		sum_dist = 0.0;
		for(i=0; i<plithos_items; i++)
		{	 
			min_distances[i] = MAX_DIST;
			for(j=0; j<thesi_k; j++)
			{
				if(distance_metric == 'h')
					dist = hammingDistance(array[kmeds[j]], array[i], dimensions);
				else if(distance_metric == 't')
					dist = array[kmeds[j]][i];
				else if(distance_metric == 'e')
					dist = eucledianDistance(array[kmeds[j]], array[i], dimensions);
				else if(distance_metric == 'c')
					dist = cosineDistance(array[kmeds[j]], array[i], dimensions);
				else// manhattan
					dist = manhattanDistance(array[kmeds[j]], array[i], dimensions);
					
				if(dist<min_distances[i])	
					min_distances[i] = dist;
			}
			sum_dist = sum_dist + min_distances[i]*min_distances[i];
		}
		
		int_sd = (int)(sum_dist);
		d_r = rand()%int_sd;

		probability = 0.0;
		for(i=0; i<plithos_items; i++)
		{
			probability = probability + min_distances[i]*min_distances[i];
			if(probability >= d_r)
			{
				kmeds[thesi_k] = i;
				break;
			}	
		}		
	}
	
	free(min_distances);
	return kmeds;
}

int* concentrate(double** array, char metriki, char type, int dimensions, int plithos_items, int K) 
{
	int i, j, t, *medoids, thesi_min;
	double *vi, *sumDist, **distances, sum, sum2, min;
	distances = (double**)malloc(plithos_items*sizeof(double));
	for(i=0;i<plithos_items;i++)
		distances[i]=(double*)malloc(plithos_items*sizeof(double));
		
	vi = (double*)malloc(plithos_items*sizeof(double));
	sumDist = (double*)malloc(plithos_items*sizeof(double));
	medoids = (int*)malloc(K*sizeof(int));
	
	if(type == 'h')
	{
		for(i=0;i<plithos_items;i++)
		{
			for(j=0;j<plithos_items;j++)
				distances[i][j] = hammingDistance(array[i], array[j], dimensions);
		}
	}
	else if(type == 'e')
	{
		if(metriki == 'e')
		{
			for(i=0;i<plithos_items;i++)
			{
				for(j=0;j<plithos_items;j++)
					distances[i][j] = eucledianDistance(array[i], array[j], dimensions);
			}
		}
		else if(metriki == 'c')
		{
			for(i=0;i<plithos_items;i++)
			{
				for(j=0;j<plithos_items;j++)
					distances[i][j] = cosineDistance(array[i], array[j], dimensions);
			}
		}
		else if(metriki == 'm')
		{
			for(i=0;i<plithos_items;i++)
			{
				for(j=0;j<plithos_items;j++)
					distances[i][j] = manhattanDistance(array[i], array[j], dimensions);
			}
		}
		else {printf("dosate lathos metriki\n"); return NULL;}
	}
	else if(type == 'm')
	{
		for(i=0;i<plithos_items;i++)
			{
				for(j=0;j<plithos_items;j++)
					distances[i][j] = array[i][j];
			}
	}
	else {printf("dosate lathos type\n");}
	
	for(i=0; i<plithos_items; i++)
	{
		sumDist[i] = 0.0;
		for(j=0; j<plithos_items; j++)
			sumDist[i] = sumDist[i] + distances[i][j];
	}
	
	for(i=0;i<plithos_items;i++)
	{
		vi[i] = 0.0;
		for(j=0;j<plithos_items;j++)
		{
			vi[i] = vi[i] + distances[i][j]/sumDist[j];
		}
	}
	
	for(i=0; i<K; i++)
	{
		min = MAX_DIST;
		thesi_min = -1;
		for(j=0; j<plithos_items; j++)
		{
			if(vi[j] < min)
			{
				min = vi[j];
				thesi_min = j;
			}
		}
		
		medoids[i] = thesi_min;
		vi[thesi_min] = MAX_DIST;
	}
	return medoids;	
}

double hammingDistance(double *item1, double *item2, int dimension)
{
	int i;
	double sum = 0;
	
	for(i=0; i<dimension; i++)
			sum = sum + fabs(item1[i] - item2[i]);	
	
	return sum;
}

double eucledianDistance(double *item1, double *item2, int dimension)
{
	int i;
	double sum = 0;
	for(i=0; i<dimension; i++)
			sum = sum + (item1[i]-item2[i])*(item1[i]-item2[i]);	
	
	return sqrt(sum);
}

double manhattanDistance(double *item1, double *item2, int dimension)
{
	int i;
	double sum = 0;
	for(i=0; i<dimension; i++)
			sum = sum + fabs(item1[i]-item2[i]);	
	
	return sum;
}

double cosineDistance(double *item1, double *item2, int dimension)
{
	int i = 0;
	double ar = 0.0, par1 = 0.0, par2 = 0.0;
	for(i=0; i<dimension; i++)
	{
		ar = ar + item1[i]*item2[i];
		par1 = par1 + item1[i]*item1[i];
		par2 = par2 + item2[i]*item2[i];
	} 	
	return (ar / (sqrt(par1)*sqrt(par2)) );
}

double* findCentroid(double** array, int plithos_items, int dimensions)
{
	int i, j;
	double *centroid, *sumdim;
	centroid = (double*)malloc(dimensions*sizeof(double));
	sumdim = (double*)malloc(dimensions*sizeof(double));
	for(i=0;i<dimensions;i++)
		sumdim[i] = 0.0;
	for(i=0;i<plithos_items;i++)
	{
		for(j=0; j<dimensions;j++)
		{
			sumdim[j] += array[i][j];
		}
	}
	for(i=0;i<dimensions;i++)
		centroid[i] = sumdim[i] / plithos_items;
	return centroid;
}

int* PAM_assignment(double** array, char distance_metric, int dimensions, int plithos_items, int *medoids, int K)
{
	int i, j, *assignment;
	double min_dist, dist;
	
	assignment = (int*)malloc(plithos_items*sizeof(int));
	
	for(i=0; i<plithos_items; i++)
	{
		assignment[i] = -1;	min_dist = MAX_DIST;
		for(j=0; j<K; j++)
		{
			if(distance_metric == 'h')
				dist = hammingDistance(array[medoids[j]], array[i], dimensions);
			else if(distance_metric == 't')
				dist = array[medoids[j]][i];
			else if(distance_metric == 'e')
				dist = eucledianDistance(array[medoids[j]], array[i], dimensions);
			else if(distance_metric == 'c')
				dist = cosineDistance(array[medoids[j]], array[i], dimensions);
			else// manhattan
				dist = manhattanDistance(array[medoids[j]], array[i], dimensions);
		
			if(dist<min_dist)
			{
				min_dist = dist;
				assignment[i] = j;
			}
		}
	}

	return assignment;
}

int* loydsNEW(double** array, int* assignments, int plithos_items, int dimensions, char distance_metric, int K, int* best_meds)
{
	int i,j,k, count, different_ass;
	double *centroid, min_dist, dist, sum;
	int *medoids, *new_ass;
	medoids = (int*)malloc(K*sizeof(int));
	centroid = (double*)malloc(dimensions*sizeof(double));
	new_ass = (int*)malloc(plithos_items*sizeof(int));
		
	while(1)
	{
		
		if(distance_metric != 't')
		{
			for(i=0;i<K;i++)
			{
				count = 0;
				for(j=0;j<dimensions;j++)
					centroid[j] = 0.0;
					
				for(j=0;j<plithos_items;j++)
				{
					if(assignments[j] == i) //an anikei to item sto cluster
					{
						count++;
						for(k=0;k<dimensions;k++)
							centroid[k] += array[j][k];
					}
				}	
			
				for(j=0;j<dimensions;j++)
					centroid[j] /= count;
			
					
				medoids[i] = -1;  min_dist = MAX_DIST;
				for(j=0; j<plithos_items; j++)
				{
					if(assignments[j] == i)
					{
						if(distance_metric == 'h') dist = hammingDistance(centroid, array[j], dimensions);
						else if(distance_metric == 'e') dist = eucledianDistance(centroid, array[j], dimensions);
						else if(distance_metric == 'c') dist = cosineDistance(centroid, array[j], dimensions);
						else dist = manhattanDistance(centroid, array[j], dimensions); //manhattan
						
						if(dist<min_dist)
						{
							min_dist = dist;
							medoids[i] = j;
						}
					}
				}
			}
		}
		else //if(distance_metric == 't')
		{
			for(k=0;k<K;k++)
			{	
				min_dist = MAX_DIST; medoids[k] = -1;
				for(i=0;i<plithos_items;i++)
				{
					if(assignments[i] == k) //an anikei sto cluster
					{
						sum = 0;
						for(j=0;j<plithos_items;j++)
							if(assignments[j] == k)
								sum+=array[i][j];
								
						if(sum<min_dist)
						{
							min_dist = sum;
							medoids[k] = i;
						}
					}	
				}
			}
		}

		new_ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, medoids, K);
		
		different_ass = 0;
		for(i=0;i<plithos_items;i++)
			if(new_ass[i] != assignments[i])
				different_ass = 1;
						
		if(different_ass == 0)
		{
			for(i=0;i<K;i++)
				best_meds[i]=medoids[i];
			free(assignments);
			return new_ass;
		}
		
		for(i=0;i<plithos_items;i++)
			assignments[i] = new_ass[i];
	}
}

int *selectMoutofN(int m, int n)
{
	int *p, i, j, exists;
	p = (int*)malloc(m*sizeof(int));	
	if(p==NULL)
	{
		printf("\nInside function selectMoutofN, malloc failed");
		exit(-1);	
	}
	
	for(i=0; i<m;)
	{
		p[i] = rand()%n;
		exists = 0;
		for(j=0; j<i; j++)
			if(p[j] == p[i])
				exists = 1;
				
		if(exists==0)
			i++;
	}
	
	return p;	
}

int* CLARA(double** array, int n, int plithos_items, int dimensions, char metriki, char type, char method, int* best_meds, int K)     //method='k' gia k_medoids, 'c' gia concentrate 
{
	int i, j, k, z, *moutofn, *medoids, **ass, *small_ass, min, count;
	double **small_array, sum, min_cost, dist, min_dist, *centroid;
	centroid = (double*)malloc(dimensions*sizeof(double));
	small_ass=(int*)malloc(n*sizeof(int));
	ass=(int**)malloc(S1*sizeof(int*));
	for(i=0;i<S1;i++)
		ass[i]=(int*)malloc(plithos_items*sizeof(int));
		
	int **meds=(int**)malloc(S1*sizeof(int*));
	for(i=0;i<S1;i++)
		meds[i]=(int*)malloc(K*sizeof(int));
		
	if(metriki!='t')
	{
		small_array=(double**)malloc(n*sizeof(double*));
		for(i=0;i<n;i++)
			small_array[i]=(double*)malloc(dimensions*sizeof(double));
	}
	else
	{
		small_array=(double**)malloc(n*sizeof(double*));
		for(i=0;i<n;i++)
			small_array[i]=(double*)malloc(n*sizeof(double));
	}
	
	moutofn=(int*)malloc(n*sizeof(int));
	medoids=(int*)malloc(K*sizeof(int));
	for(i=0;i<S1;i++)
	{
		moutofn = selectMoutofN(n, plithos_items);
		//for(j=0; j<n; j++)
		//	printf("\nmoutofn[%d]= %d", j, moutofn[j]);
		
		//upologismos small_array
		if(metriki!='t')
		{
			for(j=0;j<n;j++)
			{
				for(k=0;k<dimensions;k++)
					small_array[j][k] = array[moutofn[j]][k];
			}
		}
		else
		{
			for(j=0;j<n;j++)
			{
				for(k=0;k<n;k++)
				{
					//printf("\n %d,%d <- %d, %d", j, k, moutofn[j], moutofn[k]);
					//printf(" : %lf", array[moutofn[j]][moutofn[k]]);
					small_array[j][k] = array[moutofn[j]][moutofn[k]];
					//printf("\n %d,%d <- %d, %d", j, k, moutofn[j], moutofn[k]);
				}
			}
		}
		
		//end upologismos small_array
		
		//upologismos arxikwn medoids me small_array
		if(method=='k')
			medoids = k_medoids(small_array, metriki, dimensions, n, K);
		else if(method=='c')
			medoids = concentrate(small_array, metriki, type, dimensions, n, K);
		else {printf("\nDosate lathos methodo!\n"); exit(-1);}
		//end upologismos arxikwn medoids
		
		
		//PAM ston small_array
		small_ass = PAM_assignment(small_array, metriki, dimensions, n, medoids, K);  
		
		//upologismos newn medoids
		if(metriki!='t')
		{	
			for(j=0;j<K;j++)
			{
				count=0;
				for(k=0;k<dimensions;k++) {centroid[k] = 0.0;}
				for(k=0;k<n;k++)
				{
					if(small_ass[k] == j)//an aniksei sto cluster
					{
						count++;
						for(z=0;z<dimensions;z++)
							centroid[z]+=small_array[k][z];
					}
				} 
				for(k=0;k<dimensions;k++)
					centroid[k] /= count;
					
				medoids[j] = -1;  min_dist = MAX_DIST;
				for(k=0; k<n; k++)
					{
						if(small_ass[k] == j)
						{
							if(metriki == 'h') dist = hammingDistance(centroid, small_array[k], dimensions);
							else if(metriki == 'e') dist = eucledianDistance(centroid, small_array[k], dimensions);
							else if(metriki == 'c') dist = cosineDistance(centroid, small_array[k], dimensions);
							else dist = manhattanDistance(centroid, small_array[k], dimensions); //manhattan
							
							if(dist<min_dist)
							{
								min_dist = dist;
								medoids[j] = k;
							}
						}
					}	
			}
		}  
		else //if(distance_metric == 't')
		{
			for(j=0;j<K;j++)
			{
				min_dist = MAX_DIST; medoids[j] = -1;
				for(k=0;k<n;k++)
				{
					if(small_ass[k] == j) //an anikei sto cluster
					{
						sum = 0;
						for(z=0;z<n;z++)
							if(small_ass[z] == j)
								sum+=small_array[k][z];
								
						if(sum<min_dist)
						{
							min_dist = sum;
							medoids[j] = k;
						}
					}	
				} 
			} //for(j=0;j<K;j++) printf("\n MEDOID %d = %d",j,medoids[j]);
		}
		//end upologismos newn medoids 
		
		
		for(j=0; j<K; j++)
			medoids[j] = moutofn[medoids[j]];
		
		//PAM assignment 
		ass[i] = PAM_assignment(array, metriki, dimensions, plithos_items, medoids, K);
		
		for(j=0;j<K;j++)
		{
			meds[i][j] = medoids[j];
		}
		
		//upologismos kostous
		sum=0.0;
		min_cost = MAX_DIST;
		if(metriki!='t')
		{	
			for(j=0;j<K;j++)
			{ 
				for(k=0;k<plithos_items;k++)
				{
					if(ass[i][k] == j) //an anikei to item sto cluster
					{
						if(metriki == 'h'){ dist = hammingDistance(array[medoids[j]], array[k], dimensions); sum+=dist;}
						else if(metriki == 'e'){ dist = eucledianDistance(array[medoids[j]], array[k], dimensions);  sum+=dist;}
						else if(metriki == 'c'){ dist = cosineDistance(array[medoids[j]], array[k], dimensions);  sum+=dist;}
						else {dist = manhattanDistance(array[medoids[j]], array[k], dimensions);  sum+=dist;} //manhattan
					}
				}
			}
			if(sum<min_cost)
			{
				min_cost=sum;
				min=i;
			}	
		}
		else //metriki == 't'
		{
			for(j=0;j<K;j++)
			{
				for(k=0;k<plithos_items;k++)
				{
					if(ass[i][k] == j)
						sum += array[k][medoids[j]];
				}
			}
			if(sum<min_cost)
			{
				min_cost=sum;
				min=i;
			}
		}	
				
		//end upologismos kostous
//		printf("\n KOSTOS %lf", sum);
	}
//	printf("\n %lf, %d", min_cost, min);
	for(i=0; i<K; i++)
		best_meds[i] = medoids[i];
	
	return ass[min];
}
		
	
int *CLARANS(double** array, int plithos_items, int dimensions, char metriki, char type, char method, int* best_meds, int K, int S, int fraction)   //method='k' gia k_medoids, 'c' gia concentrate 
{
	int i, j, k, num, *medoids, exists, *ass, *best_ass, Q, q, becomes_m, becomes_t, thesi_becomes_t, epanalipseis;
	double last_cost, current_cost, min_cost, dist;
	
	ass=(int*)malloc(plithos_items*sizeof(int));
	best_ass=(int*)malloc(plithos_items*sizeof(int));	
	medoids=(int*)malloc(K*sizeof(int));
	if(S==0)
		epanalipseis=2;
	else
		epanalipseis=S;
	if(fraction!=0)
		Q=fraction;
	else
	{
		
		Q = 250;
		if(0.12*K*(plithos_items-K) > Q)
			Q = 0.12*K*(plithos_items-K);
	}
	//Q = 10;
	printf("\n Q = %d", Q);

	min_cost = MAX_DIST;
	
	for(i=0; i<epanalipseis; i++)
	{
		//upologismos arxikwn medoids
		if(method=='k')
			medoids = k_medoids(array, metriki, dimensions, plithos_items, K);
		else if(method=='c')
			medoids = concentrate(array, metriki, type, dimensions, plithos_items, K);
		else {printf("\nDosate lathos methodo!\n"); exit(-1);}
		//end upologismos arxikwn medoids
		
		ass = PAM_assignment(array, metriki, dimensions, plithos_items, medoids, K);
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		current_cost=0.0;
		for(j=0;j<K;j++)
		{ 
			for(k=0;k<plithos_items;k++)
			{
				if(ass[k] == j) //an anikei to item sto cluster
				{
					if(metriki == 't'){ dist = array[medoids[j]][k], current_cost+=dist;}
					else if(metriki == 'h'){ dist = hammingDistance(array[medoids[j]], array[k], dimensions); current_cost+=dist;}
					else if(metriki == 'e'){ dist = eucledianDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;}
					else if(metriki == 'c'){ dist = cosineDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;}
					else {dist = manhattanDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;} //manhattan
				}
			}
		}
		//printf("\n current_cost = %lf, ", current_cost);
	//	for(j=0; j<K; j++) printf(" %4d ", medoids[j]);
		if(current_cost<min_cost)
		{
			min_cost = current_cost;
			for(j=0; j<plithos_items; j++)
				best_ass[j] = ass[j];
			
			for(j=0; j<K; j++)
				best_meds[j] = medoids[j];
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		last_cost = current_cost;
		
		for(q=0; q<Q; q++)
		{	
			//swap m<->t
			thesi_becomes_t = rand()%K;
			becomes_t = medoids[thesi_becomes_t];
			while(1)
			{
				becomes_m = rand()%plithos_items;
				exists = 0;
				for(j=0; j<K; j++)
					if(medoids[j] == becomes_m)
						exists = 1;
						
				if(exists == 0)
				{
					medoids[thesi_becomes_t] = becomes_m;
					break;
				}	
			}
			
			
			ass = PAM_assignment(array, metriki, dimensions, plithos_items, medoids, K);
		
			//////////////////////////////////////////////////////////////////////////////////////////////////
			current_cost=0.0;
			for(j=0;j<K;j++)
			{ 
				for(k=0;k<plithos_items;k++)
				{
					if(ass[k] == j) //an anikei to item sto cluster
					{
						if(metriki == 't'){ dist = array[medoids[j]][k], current_cost+=dist;}
						else if(metriki == 'h'){ dist = hammingDistance(array[medoids[j]], array[k], dimensions); current_cost+=dist;}
						else if(metriki == 'e'){ dist = eucledianDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;}
						else if(metriki == 'c'){ dist = cosineDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;}
						else {dist = manhattanDistance(array[medoids[j]], array[k], dimensions);  current_cost+=dist;} //manhattan
					}
				}
			}
		//	printf("\n current_cost = %lf, ", current_cost);
		//	for(j=0; j<K; j++) printf(" %4d ", medoids[j]);
			
			if(current_cost<min_cost)
			{
				min_cost = current_cost;
				for(j=0; j<plithos_items; j++)
					best_ass[j] = ass[j];
				
				for(j=0; j<K; j++)
					best_meds[j] = medoids[j];
			}
			
			if(current_cost<last_cost)
				last_cost = current_cost;
		
			else
			{
				// epanaferw to swap
				medoids[thesi_becomes_t] = becomes_t;
			}
				
		}
	}
	return best_ass;	
}

double* Silhouette(double** array, int plithos_items, int dimensions, char metriki, int* medoids, int* ass, int K)
{
	int i, j, k, *itemsInCluster, nCluster; //nCluster = neighborCluster
	double min_dist, dist, *a, *b, *s, max; //to max einai gia sugrisi anamesa se a(i) kai b(i)
	a=(double*)malloc(plithos_items*sizeof(double));
	b=(double*)malloc(plithos_items*sizeof(double));
	s=(double*)malloc(plithos_items*sizeof(double));
	itemsInCluster=(int*)malloc(K*sizeof(int));
	
	
	//upologismos pinaka itemsInCluster
	for(i=0;i<K;i++)
	{
		itemsInCluster[i] = 0;
		for(j=0;j<plithos_items;j++)
		{
			if(ass[j] == i)
				itemsInCluster[i]++;
		}
	}
	//end upologismos pinaka itemsInCluster
				
			
	for(i=0;i<K;i++)
	{
		for(j=0;j<plithos_items;j++) 
		{
			//upologismos a gia to j item
			dist=0.0;
			if(ass[j] == i)
			{
				for(k=0;k<plithos_items;k++)
				{
					if(ass[k] == i)
						dist+=findDistance(array, j, k, dimensions, metriki);
				}
				a[j]=dist/itemsInCluster[i];
			}//end upologismos a gia to j item
			
			//upologismos nCluster
			min_dist=MAX_DIST; nCluster=-1; dist=0.0;
			for(k=0;k<K;k++)
			{
				if(ass[j] != k)
				{
					dist=findDistance(array, j, medoids[k], dimensions, metriki);
					if(dist<min_dist)
					{
						min_dist=dist;
						nCluster=k;
					}
				}
			}//end upologismos nCluster
				
			//upologismos b gia to j item
			dist=0.0;
			if(ass[j] == i)
			{
				for(k=0;k<plithos_items;k++)
				{
					if(ass[k] == nCluster)
						dist+=findDistance(array, j, k, dimensions, metriki);
				}
				b[j]=dist/itemsInCluster[nCluster];
			}//end upologismos b gia to j item
			
			max=a[j];
			if(a[j] < b[j])
				max=b[j];
				
			s[j] = (b[j] - a[j]) / max;		
		}
	}
	
	free(a); free(b); free(itemsInCluster);
	return s;
}

		

void read_configfile(char* config, int* K, int* s, int* fraction)
{
	char *dummy;
    int data;
    
    FILE *configfile;
    
    dummy=(char*)malloc(50*sizeof(char));
    
    configfile = fopen(config, "r");
    if(configfile == NULL)
    { printf("Failed to open configuration file\n"); exit(-1); }
    
	while(!feof(configfile))
    {
        fscanf(configfile, "%s", dummy);
        if(strcmp(dummy, "number_of_clusters:") == 0)
        {
        	fscanf(configfile, "%d", &data);
        	if(data > 0)
                *K = data;
            else
            	printf("LATHOS ARITHMOS CLUSTER");
		}   
		if(strcmp(dummy, "clarans_set_fraction:") == 0)
        {
            fscanf(configfile, "%d", &data);
            *fraction = data;

        }
        if(strcmp(dummy, "clarans_iterations:") == 0)
        {
            fscanf(configfile, "%d", &data);
            *s = data;
        }
    }
    
    free(dummy);
    fclose(configfile);
}

double findDistance(double** array, int item1, int item2, int dimensions, char metriki)
{
	double dist;
	if(metriki == 't'){ dist = array[item1][item2];}
	else if(metriki == 'h'){ dist = hammingDistance(array[item1], array[item2], dimensions);}
	else if(metriki == 'e'){ dist = eucledianDistance(array[item1], array[item2], dimensions);}
	else if(metriki == 'c'){ dist = cosineDistance(array[item1], array[item2], dimensions);}
	else {dist = manhattanDistance(array[item1], array[item2], dimensions);} //manhattan
	return dist;
}

