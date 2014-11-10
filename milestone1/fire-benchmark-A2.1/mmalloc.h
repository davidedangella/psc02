/*
 * mmalloc.h
 *
 *  Created on: Nov 7, 2014
 *      Author: davide
 */

#ifndef MMALLOC_H_
#define MMALLOC_H_

int allocate_local_arrays(int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
        double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
        double **BP, double **SU, int* points_count, int*** points, int** elems);

int realloc_local_array(int* nextcf, double** bs, double** be, double** bn, double** bw,
        double** bl, double** bh, double** bp, double** su);


#endif /* MMALLOC_H_ */
