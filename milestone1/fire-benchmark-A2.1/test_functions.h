/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#ifndef TEST_FUNCTIONS_H_
#define TEST_FUNCTIONS_H_

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars);

int write_pstats_exectime( int input_key,
				           int part_key,
				           int read_key,
				           int my_rank,
				           long long time_usec );

int write_pstats_partition( int input_key,
				            int part_key,
				            int my_rank,
				            int local_intc,
				            int local_extc );

int test_distribution2(char *file_in, int nintci_m, int nintcf_m, int nextci_m
		, int nextcf_m, int** lcc_m,
		double *bs_m,double *be_m,double *bn_m,double *bw_m,double *bh_m,double *bl_m,
        double *bp_m, double *su_m, int points_count_m, int** points_m, int* elems_m,
        char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars);

int test_distribution3(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars);


#endif /* TEST_FUNCTIONS_H_ */

