/**
 * @file main.f08
 * @brief This file provides you with the original implementation of pagerank.
 * Your challenge is to optimise it using OpenMP and/or MPI.
 * @author Ludovic Capelli (l.capelli@epcc.ed.ac.uk)
 **/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>

/// The number of vertices in the graph.
#define GRAPH_ORDER 1000
/// Parameters used in pagerank convergence, do not change.
#define DAMPING_FACTOR 0.85
/// The number of seconds to not exceed forthe calculation loop.
#define MAX_TIME 10

/**
 * @brief Indicates which vertices are connected.
 * @details If an edge links vertex A to vertex B, then adjacency_matrix[A][B]
 * will be 1.0. The absence of edge is represented with value 0.0.
 * Redundant edges are still represented with value 1.0.
 */
short int adjacency_matrix[GRAPH_ORDER][GRAPH_ORDER];
double max_diff = 0.0;
double min_diff = 1.0;
double total_diff = 0.0;
 
void initialize_graph(void)
{
  memset(adjacency_matrix, 0, sizeof(adjacency_matrix));
}

/**
 * @brief Calculates the pagerank of all vertices in the graph.
 * @param pagerank The array in which store the final pageranks.
 */
void calculate_pagerank(double pagerank[], int mpi_rank, int mpi_size)
{
    double initial_rank = 1.0 / GRAPH_ORDER;
 
    // Initialise all vertices to 1/n.
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        pagerank[i] = initial_rank;
    }
 
    double damping_value = (1.0 - DAMPING_FACTOR) / GRAPH_ORDER;
    double diff = 1.0;
    double global_diff;
    size_t iteration = 0;

    // MPI variables for gathering
    // Assume mpi_size > 1
    int recvcounts[mpi_size];
    int displs[mpi_size];
    int block_size = GRAPH_ORDER / (mpi_size - 1);

    displs[0] = 0;
    recvcounts[0] = block_size;
    for (int i = 1; i < mpi_size; i++) {
        displs[i] = displs[i-1] + block_size;
        recvcounts[i] = i == mpi_size-1 ? GRAPH_ORDER % (mpi_size - 1) : block_size;
    }

    // Local block size and block starting position
    int block_pos = displs[mpi_rank];
    if (mpi_rank == mpi_size - 1)
        block_size = recvcounts[mpi_size-1];

    double new_pagerank[block_size];

    // Pre-calculate the outdegree of all nodes
    int outdegrees[GRAPH_ORDER];
    memset(outdegrees, 0, sizeof(outdegrees));
    for (int i = 0; i < GRAPH_ORDER; i++)
        for (int k = 0; k < GRAPH_ORDER; k++)
            if (adjacency_matrix[i][k] == 1) outdegrees[i]++;

    // Time tracking
    double start, elapsed;
    double time_per_iteration = 0;
    if (mpi_rank == mpi_size - 1) {
        start = omp_get_wtime();
        elapsed = omp_get_wtime() - start;
    }
    // broadcast elapsed for synchronization
    MPI_Bcast(&elapsed, 1, MPI_DOUBLE, mpi_size-1, MPI_COMM_WORLD);

    // If we exceeded the MAX_TIME seconds, we stop. If we typically spend X seconds on an iteration, and we are less than X seconds away from MAX_TIME, we stop.
    while(elapsed < MAX_TIME && (elapsed + time_per_iteration) < MAX_TIME)
    {
        // double iteration_start = omp_get_wtime();
 
        memset(new_pagerank, 0, sizeof(new_pagerank));

        // Go through each destination and update it's page rank
        // using the incoming neighbour's page rank and outdegree.
        for (int i = 0; i < block_size; i++)
            for (int j = 0; j < GRAPH_ORDER; j++)
                if (adjacency_matrix[j][block_pos + i] == 1)
                    new_pagerank[i] += pagerank[j] / (double)outdegrees[j];     
 
        // Damping
        for(int i = 0; i < block_size; i++)
        {
            new_pagerank[i] = DAMPING_FACTOR * new_pagerank[i] + damping_value;
        }
 
        // Local diff calculation
        diff = 0.0;
        for(int i = 0; i < block_size; i++)
        {
            diff += fabs(new_pagerank[i] - pagerank[block_pos + i]);
        }
        // Global diff calculation on only 1 proc
        MPI_Reduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM, mpi_size-1, MPI_COMM_WORLD);
        if (mpi_rank == mpi_size - 1) {
            max_diff = (max_diff < diff) ? diff : max_diff;
            total_diff += diff;
            min_diff = (min_diff > diff) ? diff : min_diff;
        }

        // All procs gather page ranks from other procs
        MPI_Allgatherv(new_pagerank, block_size, MPI_DOUBLE,
                       pagerank, recvcounts, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);

        // Per iteration validation
        // NOTE: May replace with reduction instead if MPI communication is faster
        if (mpi_rank == mpi_size - 1) {
            double pagerank_total = 0.0;
            for(int i = 0; i < GRAPH_ORDER; i++)
            {
                pagerank_total += pagerank[i];
            }
            if(fabs(pagerank_total - 1.0) >= 1E-12)
            {
                printf("[ERROR] Iteration %zu: sum of all pageranks is not 1 but %.12f.\n", iteration, pagerank_total);
            }
        }
 
        // double iteration_end = omp_get_wtime();
        if (mpi_rank == mpi_size - 1)
            elapsed = omp_get_wtime() - start;
        // broadcast elapsed, hopefully also acts as a barrier
        MPI_Bcast(&elapsed, 1, MPI_DOUBLE, mpi_size-1, MPI_COMM_WORLD);

        iteration++;
        time_per_iteration = elapsed / iteration;
    }
    
    if (mpi_rank == mpi_size - 1)
        printf("%zu iterations achieved in %.2f seconds\n", iteration, elapsed);
}

/**
 * @brief Populates the edges in the graph for testing.
 **/
void generate_nice_graph(void)
{
    double start = omp_get_wtime();
    initialize_graph();
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        for(int j = 0; j < GRAPH_ORDER; j++)
        {
            int source = i;
            int destination = j;
            if(i != j)
            {
                adjacency_matrix[source][destination] = 1;
            }
        }
    }
    printf("%.2f seconds to generate the nice graph.\n", omp_get_wtime() - start);
}

/**
 * @brief Populates the edges in the graph for the challenge.
 **/
void generate_sneaky_graph(void)
{
    double start = omp_get_wtime();
    initialize_graph();
    for(int i = 0; i < GRAPH_ORDER; i++)
    {
        for(int j = 0; j < GRAPH_ORDER - i; j++)
        {
            int source = i;
            int destination = j;
            if(i != j)
            {
                adjacency_matrix[source][destination] = 1;
            }
        }
    }
    printf("%.2f seconds to generate the sneaky graph.\n", omp_get_wtime() - start);
}

int main(int argc, char* argv[])
{
    // We do not need argc, this line silences potential compilation warnings.
    (void) argc;
    // We do not need argv, this line silences potential compilation warnings.
    (void) argv;

    int size, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // This program has two graph generators: generate_nice_graph and generate_sneaky_graph. 
    // If you intend to submit, your code will be timed on the sneaky graph, remember to try both.

    // Get the time at the very start.
    double start = omp_get_wtime();
    
    // For now, let every MPI proc generate their own graph
    // Note: each proc can generate their own sub graph then all gather.
    generate_nice_graph();
 
    /// The array in which each vertex pagerank is stored.
    double pagerank[GRAPH_ORDER];
    calculate_pagerank(pagerank, rank, size);
 
    // Calculates the sum of all pageranks. It should be 1.0, so it can be used as a quick verification.
    if (rank == size - 1) {
        double sum_ranks = 0.0;
        for(int i = 0; i < GRAPH_ORDER; i++)
        {
            if(i % 100 == 0)
            {
                printf("PageRank of vertex %d: %.6f\n", i, pagerank[i]);
            }
            sum_ranks += pagerank[i];
        }
        printf("Sum of all pageranks = %.12f, total diff = %.12f, max diff = %.12f and min diff = %.12f.\n", sum_ranks, total_diff, max_diff, min_diff);
        double end = omp_get_wtime();
     
        printf("Total time taken: %.4f seconds.\n", end - start);
    }

    MPI_Finalize();
 
    return 0;
}
