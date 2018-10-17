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

int main(int argc, const char * argv[])
{
	FILE *infile, *qfile, *outfile;
    char dummy[30], dataType[30], distance_metric, xar, ham[MAX_HAMMING_SIZE], **itemNames, inFile[50], qFile[50], outFile[50];
    int dimensions, j, plithos_items, i,  vec, m, *KM, *CON, *ass, *MoutofN, **meds, *clarans_ass, K, k, s, fraction, soumaGAMW, *best_meds, complete=0; 
    double distance_matrix, **array, dist, cost, sil;
    int *new_ass;
	double *siloueta, *soumaAnaKlaster, oliki=0.0, olikio, olikoSUM;
	siloueta=(double*)malloc(plithos_items*sizeof(double));
	soumaAnaKlaster=(double*)malloc(K*sizeof(double));
	best_meds = (int*)malloc(K*sizeof(int));
	new_ass = (int*)malloc(plithos_items*sizeof(int));
	
    srand(time(NULL));
	
	if(argc <= 1)
    { printf("Lathos plithos orismatwn\n");  exit(-1);}
    else
    {
    	if(argc==7)
    	{
    		for(i=1; i<argc; i+=2)
       		{
	            if(argv[i][1] == 'd')
				   strcpy(inFile, argv[i+1]);
	            else if(argv[i][1] == 'c')
	                strcpy(qFile, argv[i+1]);
	            else if(argv[i][1] == 'o')
	                strcpy(outFile, argv[i+1]);
	            else
	            { printf("Lathos orisma %c\n", argv[i][1]); exit(-1); }
        	}
        }
    	else if(argc==8)	
		{
			complete=1;
			for(i=1; i<argc-1; i+=2)
       		{
	            if(argv[i][1] == 'd')
				   strcpy(inFile, argv[i+1]);
	            else if(argv[i][1] == 'c')
	                strcpy(qFile, argv[i+1]);
	            else if(argv[i][1] == 'o')
	                strcpy(outFile, argv[i+1]);
	            else
	            { printf("Lathos orisma %c\n", argv[i][1]); exit(-1); }
        	}
		} 
    }
	read_configfile(qFile, &K, &s, &fraction);
 	
	int sum[K];
	
	infile = fopen(inFile, "r");
	if(infile == NULL){ printf("Failed to open input file\n"); exit(-1); }
	
	/*----------------------- Reading from input file ----------------------*/
	
	fscanf(infile,"%s %s", dummy, dataType);
 	printf("\n %s  %s", dummy, dataType);
	
	if(dataType[0] == 'v' || dataType[0] == 'e') // vector
 	{
 		fscanf(infile,"%s %s", dummy, dataType);
 		printf("\n %s  %s", dummy, dataType);
 		distance_metric = dataType[0]; // 'e' or 'm' or 'c'
 		// vriskw diastaseis pinaka (plithos items kai dimensions)
 		
 		 fscanf(infile, " %c", &xar); // me keno, gia na min parei to '\n'
		 for(dimensions = 0; xar!='\n';)
		 {
		 	if(xar == '\t') dimensions++; //oi diastaseis xorizontai metaksi tous me tabs, opote auksanoume dimensions gia kathe /t pou sinantame
		 	fscanf(infile, "%c", &xar);
		 }
 		printf("\n DIMENSIONS = %d", dimensions);
 		
 		for(plithos_items = 0; !feof(infile); plithos_items++)
 		{
 			for(j=0; j<dimensions+1; j++) // onoma item + diastaseis tou
 				fscanf(infile, "%s", dummy);
		}
 		printf("\n plithos_items = %d", plithos_items);
	
 		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia disdiastato pinaka pou tha periexei tis times twn items
 		for(i=0; i<plithos_items; i++)
		 	array[i] = (double*)malloc(dimensions*sizeof(double));
		
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)
		 	itemNames[i] = (char*)malloc(20*sizeof(char));		
 		
 		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
 		fscanf(infile,"%s %s %s %s", dummy, dummy, dummy, dummy); //diavazoume tis protes 4 sumvoloseires pou den mas endiaferoun
 		for(i=0; i<plithos_items; i++)
 		{
 			fscanf(infile, "%s", itemNames[i]);	//diavazoume to onoma tou item kai to vazoume stin i thesi tou pinaka itemNames
 			for(j=0; j<dimensions; j++)
 				fscanf(infile, "%lf", &array[i][j]); //diavazoume tin diastasi tou ekastote i item kai tin vazoume stin j thesi tou pinaka array
		}			
 	}
	else if(dataType[0] == 'h') // hamming
	{
		distance_metric = 'h';
		for(plithos_items = -1; !feof(infile); plithos_items++) //-1 gia to prwto enter, sto telos tis 1is grammis
			fscanf(infile,"%s %s", dummy, ham);
		
		dimensions = strlen(ham); //i deuteri simvoloseira pou diavazoume se kathe grammi einai oi diastaseis tou item
			
		printf("\n plithos_items = %d, dimensions = %d", plithos_items, dimensions);
	
		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia disdiastato pinaka pou tha periexei tis times twn items
 		for(i=0; i<plithos_items; i++)
		 	array[i] = (double*)malloc(dimensions*sizeof(double));
		
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)
		 	itemNames[i] = (char*)malloc(20*sizeof(char));
		 	
		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
 		fscanf(infile,"%s %s", dummy, dummy);
 		for(i=0; i<plithos_items; i++)
 		{
 			fscanf(infile, "%s %s", itemNames[i], ham);	//se kathe grammi i proti sumvoloseira einai to onoma tou item kai i deuteri i timi tou
 			for(j=0; j<dimensions; j++)
 				array[i][j] = ham[j] - '0';
		}	
	}
	else if(dataType[0] == 'm') // matrix
	{
		distance_metric = 't'; // maTrix
		plithos_items = 0;
		
		for(i=0; !feof(infile); i++)
		{
			fscanf(infile, " %c", &xar);
			if(xar == ',') plithos_items++; //kathe fora pou vriskoume komma auksanoume to plithos items kata 1
			if(xar=='\n') break;
		}
		plithos_items++; //giati to teleutaio item den exei komma
		
		printf("\nDM plithos_items = %d\n", plithos_items);
		array = (double**)malloc(plithos_items*sizeof(double*)); //pianoume xoro gia  disdiastato pinaka pou tha periexei tis apostaseis twn items
 		for(i=0; i<plithos_items; i++) array[i] = (double*)malloc(plithos_items*sizeof(double));
		itemNames = (char**)malloc(plithos_items*sizeof(char*)); //pianoume xoro gia pinaka pou tha periexei ta onomata twn items
 		for(i=0; i<plithos_items; i++)itemNames[i] = (char*)malloc(20*sizeof(char));
 		
 		fseek(infile, 0, SEEK_SET); // pame to deikti tou arxeiou apostasi 0 apo tin arxi
		fscanf(infile,"%s %s %s", dummy, dummy, dummy); //oi treis protes sumvoloseires den mas endiaferoun
		
		for(i=0, j=0; i<plithos_items-1;)
		{
			fscanf(infile, " %c", &xar);
			if(xar==',') //an diavasoume komma, simainei oti exoume diavasei ena olokliro itenName
			{
				dummy[j]='\0'; //opote vazoume ton teleutaio xaraktira tou itemName na einai \0 gia na exoume simvoloseira
				strcpy(itemNames[i], dummy); //antigrafi tou itemName pou exoume diavasei mexri tora stin i thesi tou pinaka itemNames
				j = 0; //vazoume to j=0 giati tha ksekinisoume na diavazoume neo onoma
				i++; //pame stin epomeni thesi tou pinaka itemNames
			}
			else //an den diavasoume komma, simainei oti diavazoume xaraktira pou anikei sto onoma tou item
			{
				dummy[j] = xar;
				j++;
			}
		}
		fscanf(infile,"%s", itemNames[plithos_items-1]); //sto teleutaio item den sunantame komma, opote to kanoume me to xeri	
		printf("\n LAST ITEM= %s\n", itemNames[plithos_items-1]);
		
		for(i=0; i<plithos_items; i++)
 		{
 			for(j=0; j<plithos_items; j++)
 				fscanf(infile, "%lf", &array[i][j]);	//gemizoume ton disdiastato pinaka me tis apostaseis ton items
		}	
	}
	else // wrong fileType
	{ printf("\nWrong input file type\n Exiting program :-( \n"); exit(-1); }
	
	/*------------------ End of reading from input file --------------------*/	
	
	/*============================================================================*/
	outfile = fopen(outFile, "w");
	if(outfile == NULL)
   	{ printf("Failed to open output file\n"); exit(-1); }
    
	fprintf(outfile, "============ I1A1U1 ============\n");
	KM = k_medoids(array, distance_metric, dimensions, plithos_items, K);
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, KM, K);
	new_ass = loydsNEW( array,  ass, plithos_items, dimensions, distance_metric, K, best_meds);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
	for(i=0; i<plithos_items; i++)
		sum[new_ass[i]] += 1;
	for(j=0; j<K; j++) 
	{
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	}
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}

	
	fprintf(outfile, "\n\n============ I1A1U2 ============\n");
	KM = k_medoids(array, distance_metric, dimensions, plithos_items, K);
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, KM, K);
	new_ass = CLARA(array, 40+2*K, plithos_items, dimensions,  distance_metric, dataType[0], 'k', best_meds, K);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
		for(i=0; i<plithos_items; i++)
			sum[new_ass[i]] += 1;
			
			
	for(j=0; j<K; j++) 
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}
	
	fprintf(outfile, "\n\n============ I1A1U3 ============\n");
	KM = k_medoids(array, distance_metric, dimensions, plithos_items, K);
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, KM, K);
	new_ass = CLARANS(array, plithos_items, dimensions,  distance_metric, dataType[0], 'k', best_meds, K, s, fraction);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
		for(i=0; i<plithos_items; i++)
			sum[new_ass[i]] += 1;
			
			
	for(j=0; j<K; j++) 
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}
	
	
	fprintf(outfile, "\n\n============ I2A1U1 ============\n");
	CON = concentrate(array, distance_metric, dataType[0], dimensions, plithos_items, K);				
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, CON, K);
	new_ass = loydsNEW( array,  ass, plithos_items, dimensions, distance_metric, K, best_meds);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
		for(i=0; i<plithos_items; i++)
			sum[new_ass[i]] += 1;
			
			
	for(j=0; j<K; j++) 
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}
	
	fprintf(outfile, "\n\n============ I2A1U2 ============\n");
	CON = concentrate(array, distance_metric, dataType[0], dimensions, plithos_items, K);				
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, CON, K);
	new_ass = CLARA(array, 40+2*K, plithos_items, dimensions,  distance_metric, dataType[0], 'c', best_meds, K);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
		for(i=0; i<plithos_items; i++)
			sum[new_ass[i]] += 1;
					
	for(j=0; j<K; j++) 
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}
	
	fprintf(outfile, "\n\n============ I2A1U3 ============\n");
	CON = concentrate(array, distance_metric, dataType[0], dimensions, plithos_items, K);				
	ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, CON, K);
	new_ass = CLARANS(array, plithos_items, dimensions,  distance_metric, dataType[0], 'k', best_meds, K, s, fraction);
	siloueta = Silhouette(array, plithos_items, dimensions, distance_metric, best_meds, new_ass, K);
	olikoSUM = 0.0;
	for(i=0; i<K; i++)
			sum[i] = 0;
		for(i=0; i<plithos_items; i++)
			sum[new_ass[i]] += 1;
			
			
	for(j=0; j<K; j++) 
		fprintf(outfile, "\nCLUSTER-%d {size: %d, medoid: %d}", j, sum[j], best_meds[j]);
	
	for(i=0;i<K;i++)
	{
		soumaGAMW=0;
		for(j=0;j<plithos_items;j++)
		{
			if(new_ass[j] == i)
			{
				soumaGAMW++;
				soumaAnaKlaster[i]+=siloueta[j];
			}
		}
		soumaAnaKlaster[i]/=soumaGAMW;
	}
	fprintf(outfile, "\nSilhouette per cluster: [");
	for(i=0;i<K;i++)
		fprintf(outfile, "%lf, ", soumaAnaKlaster[i]);
	fprintf(outfile, "]");
	for(i=0;i<K;i++)
		olikoSUM+=soumaAnaKlaster[i];
	fprintf(outfile, " Sum: %lf", olikoSUM);
	
	if(complete==1)
	{
		fprintf(outfile, "\n");
		for(j=0;j<K;j++)
		{ 
			fprintf(outfile, "\n");
			fprintf(outfile, "\nCLUSTER-%d {", j);
			for(k=0;k<plithos_items;k++)
			{
				if(new_ass[k] == j) //an anikei to item sto cluster
					fprintf(outfile, " item%d", k);
			}
			fprintf(outfile, "}", j);
		}
	}

	fclose(outfile);
	
	/*-------------- Deallocate memory -------------*/
	free(ass); free(CON); free(KM); free(new_ass); free(siloueta);
	for(i=0; i<plithos_items; i++)                //free itemNames, array
	{
		free(itemNames[i]); free(array[i]);
		itemNames[i] = NULL; array[i] = NULL;
	}
	free(itemNames); free(array);
	itemNames = NULL; array = NULL;	
		
	return 0;  
}

















	
