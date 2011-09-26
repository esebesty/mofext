#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <math.h>*/
#include "matrix.h"
#include "mystrings.h"
#include "similar.h"
#include "extend.h"
#include <ctype.h>

/******************************************************/
void    basefreq(char *seq, double *bases){ /* This is an undocumented feature calculate the base frequency */
unsigned int i,k;
int a=0,t=0,g=0,c=0;

k = strlen(seq);

for (i=0; i < k; i++){/*FIXME This is working only with DNA matrix */
  if (toupper(seq[i]) =='A') a++;
  if (toupper(seq[i]) =='G') g++;
  if (toupper(seq[i]) =='T') t++;
  if (toupper(seq[i]) =='C') c++;
  if (seq[i] == 'R') {a++;g++;}
  if (seq[i] == 'Y') {c++;t++;}
  if (seq[i] == 'S') {g++;c++;}
  if (seq[i] == 'W') {a++;t++;}
  if (seq[i] == 'K') {g++;t++;}
  if (seq[i] == 'M') {a++;c++;}
  if (seq[i] == 'N') {a++;c++;t++;g++;}
}

bases[0] = (double)a / (double)k;
bases[1] = (double)c / (double)k;
bases[2] = (double)g / (double)k;
bases[3] = (double)t / (double)k;
return;
}


/******************************************************/
double pvalue(char *seq, double *bases){ /* Undocumented feature. Calculate the probability. Need more specificity */
unsigned int    i,k;
double          p=1.0;

k = strlen(seq);

for (i=0; i < k; i++){ /*FIXME  This is working only with DNA matrix */
	if (seq[i] == 'A') p*=  bases[0];
	if (seq[i] == 'C') p*=  bases[1];
	if (seq[i] == 'G') p*=  bases[2];
	if (seq[i] == 'T') p*=  bases[3];
	if (seq[i] == 'a') p*=  bases[0] * 0.8;
	if (seq[i] == 'c') p*=  bases[1] * 0.8;
	if (seq[i] == 'g') p*=  bases[2] * 0.8;
	if (seq[i] == 't') p*=  bases[3] * 0.8;
	if (seq[i] == 'R') p*= (bases[0] + bases[2]);
	if (seq[i] == 'Y') p*= (bases[1] + bases[3]);
	if (seq[i] == 'S') p*= (bases[2] + bases[1]);
	if (seq[i] == 'W') p*= (bases[0] + bases[3]);
	if (seq[i] == 'K') p*= (bases[2] + bases[3]);
	if (seq[i] == 'M') p*= (bases[0] + bases[1]);
}

return(p);
}


/******************************************************/
double sumscore(char *seqA, char *seqB, double *percent){
  int      i,k;
  double   sum=0, self=0;

  k = strlen(seqB);

  for(i=0; i < k; i++){
         sum  += getscore(seqA[i],seqB[i]);
         self += getscore(seqB[i],seqB[i]);
  }
  *percent = sum / self;
  return(sum);
}

/******************************************************/
void similarity(char *id, char *seqA, char *seqB, int wordsize, double limit, double *bases, char *output){
  int m,n;
  int i,j,k;
  char *seqA_f, *seqB_f, *overA, *overB, *subseq, *subpatt;
  double score,extscore;        /* The seqA_f <-> seqB_f score
                                 * and the extended score  */
  double percent;               /* The self/score percentage */

  int startA,startB,lenA,lenB;  /* This variables need by GetCorrectPos */
  int comm_start, comm_len;
  double p;

m = strlen(seqA) - wordsize + 1;
n = strlen(seqB) - wordsize + 1;

seqA_f  = malloc(wordsize+1);
seqB_f  = malloc(wordsize+1);
overA   = malloc(strlen(seqA)+1);
overB   = malloc(strlen(seqB)+1);
subseq  = malloc(strlen(seqA)+1); 
subpatt = malloc(strlen(seqB)+1);

if((seqA_f == NULL) || (seqB_f == NULL) ||
   (overA  == NULL) || (overB  == NULL) ||
   (subseq == NULL) || (subpatt== NULL)) {
	fprintf(stderr,"Memory is not enought\n");
	exit(2);
}

for (i=0; i<m; i++){
  substr(seqA,i,wordsize,seqA_f);
  for (j=0; j<n; j++){

      substr(seqB,j,wordsize,seqB_f);
      score = sumscore(seqA_f,seqB_f,&percent);

      if (percent >= limit){
	startA = i+1;
	startB = j+1;
	lenA   = strlen(seqA);
	lenB   = strlen(seqB);

	GetCorrectPos(&startA,&startB,&lenA,&lenB);

	substr(seqA,startA,lenA,overA);
	substr(seqB,startB,lenB,overB);
	extscore = ExtScore(overA, overB, wordsize, &comm_start, &comm_len);

	if (extscore != score){
		substr(overA,comm_start,comm_len,subseq);
		substr(overB,comm_start,comm_len,subpatt);
		startA = startA+comm_start+1;
		startB = startB+comm_start+1;
		i+=comm_len;
		j+=comm_len;
	}
	else{
		substr(seqA,i,wordsize,subseq);  
		substr(seqB,j,wordsize,subpatt); 
		startA = i+1;
		startB = j+1;
	}
	p = pvalue(subpatt,bases);

	for (k=0; k<12; k++){
		switch(output[k]){
		case 'i':
			printf("%s ",id);
			break;
		case 's':
			printf("%.2f ",score);
			break;
		case 'e':
			printf("%.2f ",extscore);
			break;
		case 'p': /* Print out the undocumented probability */
			printf("%f ",p);
			break;
		case 'd':
			printf("%s ",subseq);
			break;
		case 'D':
			printf("%s ",seqA);
			break;
		case 'q':
			printf("%s ",subpatt);
			break;
		case 'Q':
			printf("%s ",seqB);
			break;
		case 'f':
			printf("%d ",startB);
			break;
		case 'F':
			printf("%d ",startA);
			break;
		case '%':
			printf("%.2f ",percent);
		};
	}
	printf("\n");
      }
  }
}

free(subpatt);
free(subseq);
free(overB);
free(overA);
free(seqA_f);
free(seqB_f);
return;
}/* similarity */
