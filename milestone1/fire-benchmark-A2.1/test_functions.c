/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "util_read_files.h"
#include "util_write_files.h"
#include "test_functions.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars) {
    
    // global sized variables, for reading the input file
    int nintci_m, nintcf_m;  
    int nextci_m, nextcf_m;
    int **lcc_m; 
    double *bs_m, *be_m, *bn_m, *bw_m, *bh_m, *bl_m;
    double *bp_m;  
    double *su_m;  
    int points_count_m;  
    int** points_m;  
    int* elems_m;  
    int i;

    // read the entire file
    int f_status = read_binary_geo( file_in, &nintci_m, &nintcf_m, &nextci_m,
            &nextcf_m, &lcc_m, &bs_m, &be_m, &bn_m, &bw_m, &bl_m, &bh_m, &bp_m,
            &su_m, &points_count_m, &points_m, &elems_m );
    if ( f_status != 0 ) {
        printf( "Error in reading input file \n" );
        return -1;
    }

    // allocate distribution vector
    double *distr;
    if ( ( distr = (double *) malloc( ( nintcf_m + 1 ) * sizeof(double) ) )
            == NULL ) {
        printf( "malloc failed to allocate distr array" );
        return -1;
    }
    for ( i = nintci_m; i < ( nintcf_m + 1 ); i++ ) {
        distr[i] = 0.0;
    }

    // copy the local values using the generated map
    for ( i = 0; i < local_num_elems; i++ ) {
        distr[local_global_index[i]] = scalars[i];
    }
    // write vtk file
    vtk_write_unstr_grid_header( file_in, file_vtk_out, nintci_m, nintcf_m,
            points_count_m, points_m, elems_m );
    vtk_append_double( file_vtk_out, "SCALARS", nintci_m, nintcf_m, distr );
    printf( "Distribution VTK file succesfully generated! \n" );

    // free the allocated memory
    free( su_m );
    free( bp_m );
    free( bh_m );
    free( bl_m );
    free( bw_m );
    free( bn_m );
    free( be_m );
    free( bs_m );
    free( elems_m );

    for ( i = 0; i < nintcf_m + 1; i++ ) {
        free( lcc_m[i] );
    }
    free( lcc_m );

    for ( i = 0; i < points_count_m; i++ ) {
        free( points_m[i] );
    }
    free( points_m );
    free( distr );

    return 0;
}



void compute_local_elements_nodes( int nintci, int nintcf, int* global_elems, int** global_points, int global_points_count, int** local_points, int* local_elems, int* local_points_count, int* local_to_global){
	int i, j;//, k;
	int global_elem;
	int global_node;
	int local_node_ix = 0;

	int * points_global_to_local = (int*) malloc((global_points_count)*sizeof(int));
	// points_to_local: will have index from 0 to max_index of the point contained in the domain
	// the rest will be -1 --> point with global index [i] not in the local domain

	for(i=0;i<global_points_count;i++)
		points_global_to_local[i]=-1;


	for(i=nintci; i<nintcf+1; i++){

		global_elem = local_to_global[i]*8;

		for(j=0; j<8; j++){
			global_node = global_elems[ global_elem + j ];

			if( points_global_to_local[  global_node   ] == -1 ) {
				points_global_to_local[  global_node   ] = local_node_ix;

					memcpy(local_points[local_node_ix], global_points[ global_node ], 3*sizeof(int));

				local_node_ix++;
			}

			local_elems[i*8+j] = points_global_to_local[ global_node ];
		}

	}

	*local_points_count=local_node_ix;

	free(points_global_to_local);
}


int test_distribution3(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars) {

    // global sized variables, for reading the input file
    int nintci_m, nintcf_m;
    int nextci_m, nextcf_m;
    int **lcc_m;
    double *bs_m, *be_m, *bn_m, *bw_m, *bh_m, *bl_m;
    double *bp_m;
    double *su_m;
    int points_count_m;
    int** points_m;
    int* elems_m;
    int i;

    // read the entire file
    int f_status = read_binary_geo( file_in, &nintci_m, &nintcf_m, &nextci_m,
            &nextcf_m, &lcc_m, &bs_m, &be_m, &bn_m, &bw_m, &bl_m, &bh_m, &bp_m,
            &su_m, &points_count_m, &points_m, &elems_m );
    if ( f_status != 0 ) {
        printf( "Error in reading input file \n" );
        return -1;
    }

    // allocate distribution vector
    double *distr;
    if ( ( distr = (double *) malloc( ( nintcf_m + 1 ) * sizeof(double) ) )
            == NULL ) {
        printf( "malloc failed to allocate distr array" );
        return -1;
    }
    for ( i = nintci_m; i < ( nintcf_m + 1 ); i++ ) {
        distr[i] = 0.0;
    }

    // copy the local values using the generated map
    for ( i = 0; i < local_num_elems; i++ ) {
        distr[local_global_index[i]] = scalars[i];
    }

    int ** local_points = malloc(points_count_m*sizeof(int*));
    for(i=0;i<points_count_m;i++)
    	local_points[i] = malloc(3*sizeof(int));

    int * local_elems = malloc(local_num_elems*8*sizeof(int));
    int local_points_count;


    compute_local_elements_nodes( 0, local_num_elems-1, elems_m, points_m, points_count_m, local_points, local_elems, &local_points_count, local_global_index);

    // write vtk file
    vtk_write_unstr_grid_header( file_in, file_vtk_out, 0, local_num_elems-1,
    		local_points_count, local_points, local_elems );
    vtk_append_double( file_vtk_out, "SCALARS", 0, local_num_elems-1, scalars );
    printf( "Distribution VTK file succesfully generated! \n" );

    // free the allocated memory
    free( su_m );
    free( bp_m );
    free( bh_m );
    free( bl_m );
    free( bw_m );
    free( bn_m );
    free( be_m );
    free( bs_m );
    free( elems_m );

    for ( i = 0; i < nintcf_m + 1; i++ ) {
        free( lcc_m[i] );
    }
    free( lcc_m );

    for ( i = 0; i < points_count_m; i++ ) {
        free( points_m[i] );
    }
    free( points_m );
    free( distr );

    return 0;
}


int test_distribution2(char *file_in, int nintci_m, int nintcf_m, int nextci_m
		, int nextcf_m, int** lcc_m,
		double *bs_m,double *be_m,double *bn_m,double *bw_m,double *bh_m,double *bl_m,
        double *bp_m, double *su_m, int points_count_m, int** points_m, int* elems_m,
        char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars) {

    // global sized variables, for reading the input file
    /*int nintci_m, nintcf_m;
    int nextci_m, nextcf_m;
    int **lcc_m;
    double *bs_m, *be_m, *bn_m, *bw_m, *bh_m, *bl_m;
    double *bp_m;
    double *su_m;
    int points_count_m;
    int** points_m;
    int* elems_m;
    int i;*/


    // allocate distribution vector
    /*double *distr;
    if ( ( distr = (double *) malloc( ( nintcf_m + 1 ) * sizeof(double) ) )
            == NULL ) {
        printf( "malloc failed to allocate distr array" );
        return -1;
    }
    for ( i = nintci_m; i < ( nintcf_m + 1 ); i++ ) {
        distr[i] = 0.0;
    }

    // copy the local values using the generated map
    for ( i = 0; i < local_num_elems; i++ ) {
        distr[local_global_index[i]] = scalars[i];
    }*/
    // write vtk file
    vtk_write_unstr_grid_header( file_in, file_vtk_out, nintci_m, nintcf_m,
            points_count_m, points_m, elems_m );
    vtk_append_double( file_vtk_out, "SCALARS", nintci_m, nintcf_m, scalars );
    printf( "Distribution VTK file succesfully generated! \n" );

    // free the allocated memory
    /*free( su_m );
    free( bp_m );
    free( bh_m );
    free( bl_m );
    free( bw_m );
    free( bn_m );
    free( be_m );
    free( bs_m );
    free( elems_m );

    for ( i = 0; i < nintcf_m + 1; i++ ) {
        free( lcc_m[i] );
    }
    free( lcc_m );

    for ( i = 0; i < points_count_m; i++ ) {
        free( points_m[i] );
    }
    free( points_m );
    free( distr );
*/

    return 0;
}

/* Write statistics to pstats.dat 
 * @param input_key: 1 - tjunc
 					 2 - drall
 					 3 - pent
 					 4 - cojack
 * @param part_key: 1 - classic
 					2 - dual
 					3 - nodal
 * @param read_key: 1 - oneread
 				    2 - allread
 * @param my_rank: current process rank
 * @param time_usec: execution time in microseconds
 * @return
 */

int write_pstats_exectime( int input_key,
				           int part_key,
				           int read_key,
				           int my_rank,
				           long long time_usec ){

	// append to existing file
    FILE *fp = fopen( "pstats.dat", "a" );
    if( fp == NULL ){
        printf( "Error opening file pstats.dat for writing\n" );
        return -1;
    }


    fprintf( fp, "EXECTIME %d %d %d %d %lld\n", input_key, part_key, read_key, my_rank, time_usec );
    fclose( fp );

    return 0;
}

/* Write statistics to pstats.dat 
 * @param input_key: 1 - tjunc
 					 2 - drall
 					 3 - pent
 					 4 - cojack
 * @param part_key: 1 - classic
 					2 - dual
 					3 - nodal
 * @param my_rank: current process rank
 * @param local_intc: number of local internal cells
 * @param local_extc: number of local external cells
 * @return
 */

int write_pstats_partition( int input_key,
				            int part_key,
				            int my_rank,
				            int local_intc,
				            int local_extc ){

	// append to existing file
    FILE *fp = fopen( "pstats.dat", "a" );
    if( fp == NULL ){
        printf( "Error opening file pstats.dat for writing\n" );
        return -1;
    }


    fprintf( fp, "PARTITION %d %d %d %d %d\n", input_key, part_key, my_rank, local_intc, local_extc );
    fclose( fp );

    return 0;
}
