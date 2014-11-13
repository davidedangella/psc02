
#include <stdio.h>
#include <stdlib.h>

//void mmalloc(void const ** arr, const size_t size, char const * const name){
//	 if ( (*arr = malloc(size)) == NULL ) {
//		fprintf(stderr, "malloc failed to allocate first dimension of %s", name);
//		exit(1);
//	}
//}



int allocate_local_arrays(int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
        double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
        double **BP, double **SU, int* points_count, int*** points, int** elems, int **local_to_global){

	int i;

	 // allocating LCC
	    if ( (*LCC = (int**) malloc((*NINTCF + 1) * sizeof(int*))) == NULL ) {
	        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
	        return -1;
	    }

	    for ( i = 0; i < *NINTCF + 1; i++ ) {
	        if ( ((*LCC)[i] = (int *) malloc(6 * sizeof(int))) == NULL ) {
	            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
	            return -1;
	        }
	    }
	    // allocate other arrays
	    if ( (*BS = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BS) failed\n");
	        return -1;
	    }

	    if ( (*BE = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BE) failed\n");
	        return -1;
	    }

	    if ( (*BN = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BN) failed\n");
	        return -1;
	    }

	    if ( (*BW = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BW) failed\n");
	        return -1;
	    }

	    if ( (*BL = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BL) failed\n");
	        return -1;
	    }

	    if ( (*BH = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BH) failed\n");
	        return -1;
	    }

	    if ( (*BP = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(BP) failed\n");
	        return -1;
	    }

	    if ( (*SU = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "malloc(SU) failed\n");
	        return -1;
	    }

	    if ( (*elems = (int*) malloc((*NINTCF + 1) * 8 * sizeof(int))) == NULL ) {
	           fprintf(stderr, "malloc failed to allocate elems");
	           return -1;
	   }

	    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL ) {
	            fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
	            return -1;
		}

		for ( i = 0; i < *points_count; i++ ) {
			if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
				fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
				return -1;
			}
		}

	    if ( (*local_to_global = (int *) malloc((*NEXTCF + 1) * sizeof(int))) == NULL ) {
	        fprintf(stderr, "malloc(local_to_global) failed\n");
	        return -1;
	    }

	    return 0;
}


int realloc_local_array(int* nextcf, double** bs, double** be, double** bn, double** bw,
        double** bl, double** bh, double** bp, double** su, int** local_to_global){

	 if ( (*bs = (double *) realloc(*bs, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BS) failed\n");
	        return -1;
	    }

	    if ( (*be = (double *) realloc(*be, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BE) failed\n");
	        return -1;
	    }

	    if ( (*bn = (double *) realloc(*bn, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BN) failed\n");
	        return -1;
	    }

	    if ( (*bw = (double *) realloc(*bw, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BW) failed\n");
	        return -1;
	    }

	    if ( (*bl = (double *) realloc(*bl, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BL) failed\n");
	        return -1;
	    }

	    if ( (*bh = (double *) realloc(*bh, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BH) failed\n");
	        return -1;
	    }

	    if ( (*bp = (double *) realloc(*bp, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(BP) failed\n");
	        return -1;
	    }

	    if ( (*su = (double *) realloc(*su, (*nextcf + 1) * sizeof(double))) == NULL ) {
	        fprintf(stderr, "realloc(SU) failed\n");
	        return -1;
	    }

	    if ( (*local_to_global = (int *) realloc(*su, (*nextcf + 1) * sizeof(int))) == NULL ) {
	        fprintf(stderr, "realloc(local_to_global) failed\n");
	        return -1;
	    }



	    return 0;
}





