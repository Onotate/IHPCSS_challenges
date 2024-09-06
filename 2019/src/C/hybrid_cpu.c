/**
 * @file hybrid_cpu.c
 * @brief Contains the MPI + OpenMP version of Laplace.
 * @note This code was originaly written by John Urbanic for PSC 2014, later modified by Ludovic Capelli.
 * @author John Urbanic
 * @author Ludovic Capelli
 **/

#include <stdio.h> // printf
#include <stdlib.h> // EXIT_FAILURE
#include <math.h> // fabs
#include <mpi.h> // MPI_*
#include <string.h> // strcmp
#include "util.h"  

/**
 * @brief Runs the experiment.
 * @pre The macro 'ROWS' contains the number of rows (excluding boundaries) per MPI process. It is a define passed as a compilation flag, see makefile.
 * @pre The macro 'COLUMNS' contains the number of columns (excluding boundaries). It is a define passed as a compilation flag, see makefile.
 **/
int main(int argc, char *argv[])
{
	// Temperature grid.
	double temperature[ROWS+2][COLUMNS+2];
	// Temperature grid from last iteration
	double temperature_last[ROWS+2][COLUMNS+2]; 
	// Current iteration.
    int iteration = 0;
    // Temperature change for our MPI process
    double dt;
    // Temperature change across all MPI processes
    double dt_global = 100;
    // The number of MPI processes in total
    int comm_size;
    // The rank of my MPI process
    int my_rank;
    // Status returned by MPI calls
    MPI_Status status;
    // Sender and receiver rank for Sendrecv operation
    int from, to;

    double timer_simulation;

    // The usual MPI startup routines
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if(provided < MPI_THREAD_FUNNELED)
    {
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if(strcmp(VERSION_RUN, "hybrid_small") == 0 && comm_size != 2)
    {
        printf("The small version is meant to be run with 2 MPI processes, not %d.\n", comm_size);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else if(strcmp(VERSION_RUN, "hybrid_big") == 0 && comm_size != 8)
    {
        printf("The big version is meant to be run with 8 MPI processes, not %d.\n", comm_size);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if(my_rank == 0)
    {
        printf("Running on %d MPI processes\n\n", comm_size);
    }

    // Initialise temperatures and temperature_last including boundary conditions
    initialise_temperatures(temperature, temperature_last);

    ///////////////////////////////////
    // -- Code from here is timed -- //
    ///////////////////////////////////
    if(my_rank == 0)
    {
        start_timer(&timer_simulation);
        #pragma omp parallel
        {
            printf("Application run using %d OpenMP threads.\n", omp_get_num_threads());
        } // End of OpenMP parallel region
    }

    while(dt_global > MAX_TEMP_ERROR && iteration <= MAX_NUMBER_OF_ITERATIONS)
    {
        iteration++;

        // Main calculation: average my four neighbours
	    #pragma omp parallel for default(none) shared(temperature, temperature_last)
		for(unsigned int i = 1; i <= ROWS; i++)
		{
			for(unsigned int j = 1; j <= COLUMNS; j++)
			{
				temperature[i][j] = 0.25 * (temperature_last[i-1][j  ] +	
											temperature_last[i  ][j+1] +
											temperature_last[i  ][j-1] +
                      temperature_last[i+1][j  ]);
			}
		}

        //////////////////////
        // HALO SWAP PHASE //
        ////////////////////

        // Send down, receive up phase 
        to = my_rank == comm_size-1 ? MPI_PROC_NULL : my_rank + 1;
        from = my_rank == 0 ? MPI_PROC_NULL : my_rank - 1;
        MPI_Sendrecv(&temperature[ROWS][1], COLUMNS, MPI_DOUBLE, to, 0,
                     &temperature_last[0][1], COLUMNS, MPI_DOUBLE, from, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);

        // Send up, receive down phase
        to = my_rank == 0 ? MPI_PROC_NULL : my_rank - 1;
        from = my_rank == comm_size-1 ? MPI_PROC_NULL : my_rank + 1;
        MPI_Sendrecv(&temperature[1][1], COLUMNS, MPI_DOUBLE, to, 0,
                     &temperature_last[ROWS+1][1], COLUMNS, MPI_DOUBLE, from, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);

        //////////////////////////////////////
        // FIND MAXIMAL TEMPERATURE CHANGE //
        ////////////////////////////////////
        dt = 0.0;

		#pragma omp parallel for default(none) shared(temperature, temperature_last) reduction(max:dt)
        for(unsigned int i = 1; i <= ROWS; i++)
        {
            for(unsigned int j = 1; j <= COLUMNS; j++)
            {
    	        dt = fmax(fabs(temperature[i][j]-temperature_last[i][j]), dt);
    	        temperature_last[i][j] = temperature[i][j];
            }
        }

        // We know our temperature delta, we now need to sum it with that of other MPI processes
        // TODO: investigate
        MPI_Reduce(&dt, &dt_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Periodically print test values
        if((iteration % PRINT_FREQUENCY) == 0)
        {
            if(my_rank == comm_size - 1)
            {
                track_progress(iteration, temperature);
    	    }
        }
    }

    // Slightly more accurate timing and cleaner output 
    MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////////////////////////////
    // -- Code from here is no longer timed -- //
    /////////////////////////////////////////////
    if(my_rank == 0)
    {
        stop_timer(&timer_simulation);
        print_summary(iteration, dt_global, timer_simulation);
    }

	// Print the halo swap verification cell value 
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == comm_size - 2)
	{
		printf("Value of halo swap verification cell [%d][%d] is %.18f\n", ROWS_GLOBAL - ROWS - 1, COLUMNS - 1, temperature[ROWS][COLUMNS]);
	}

    MPI_Finalize();
}
