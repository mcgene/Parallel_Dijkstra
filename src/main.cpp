#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<stack>
#include<mpi.h>
#include<omp.h>
#include "CycleTimer.h"
using namespace std;

/*global variable*/
const int infinite = 2139062143;

/*function prototype*/
void read_data(char *file_name, int node_num, int **arr);
void Usage(char* prog_name);
void Dijstra(int *arr, int num_nodes, int start, int end, int *precursor, int *length);
void Parallel_Dijstra(int *arr, int num_nodes, int start, int end, int *precursor, int *length);
void output(int start,int end, int *precursor, int length);

/*--------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    /*check parameter*/
    if (argc != 5) Usage(argv[0]);

    /*data*/
    int *data;
    int num_node = atoi(argv[1]);
    int start = atoi(argv[3]);//start node
    int end = atoi(argv[4]);// end note
    int *precursor = (int*)malloc(sizeof(int) * num_node );
    int *pprecursor = (int*)malloc(sizeof(int) * num_node );
    int pmin_length;//just omp
    int min_length;
    int  comm_sz; //communicator size
    int  my_rank;// my rank
    double startTime, endTime, pstartTime, pendTime;

    /*MPI Init*/
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0)
    {
        //double rstartTime = CycleTimer::currentSeconds();
        read_data(argv[2], num_node, &data);
        //double rendTime = CycleTimer::currentSeconds();
        //printf("read time: %f\n", rendTime - rstartTime );

        /*Serial version*/
        //startTime = CycleTimer::currentSeconds();
        //Dijstra(data, num_node, start, end, precursor, &min_length);
        //endTime = CycleTimer::currentSeconds();
        //printf("Serial time: %f\n", endTime - startTime );

        /*Parallel version*/
        //pstartTime = CycleTimer::currentSeconds();
        Parallel_Dijstra(data, num_node, start, end, pprecursor, &pmin_length);
        //pendTime = CycleTimer::currentSeconds();
        //printf("OMP Parallel time: %f\n", pendTime - pstartTime );

        /*check the result*/
        //if( pmin_length == min_length )
            //printf("Match\n");
        //else
            //printf("OMP unmatched\n");

        /*Speed up*/
        //printf("(Compared with Serial)OpenMP Speedup is %f\n", (endTime - startTime)*1.0/(pendTime - pstartTime));
        output(start,end,pprecursor,pmin_length);
    }
    MPI_Finalize();

    return 0;
}/* main */


/*--------------------------------------------------------------------
 * Function:    read_data
 * Purpose:     read graph from file, the data in file should be organize into "node node weight"
 * In arg:      file_name
 *              node_num: the node number of graph
 * Out arg:     arr
 */
void read_data(char *file_name, int node_num, int **arr)
{
    int x, y;
    FILE *f = fopen(file_name,"r");
    *arr = (int*)malloc(sizeof(int)*node_num * node_num);

    memset(*arr, 127, sizeof(int) * node_num * node_num);

    while(!feof(f))
    {
        fscanf(f, "%d", &x);
        fscanf(f, "%d", &y);
        fscanf(f, "%d", &(*arr)[x * node_num + y]);
    }

    fclose(f);
}/*read_data*/

/*--------------------------------------------------------------------
 * Function:    Usage
 * Purpose:     Print command line for function and terminate
 * In arg:      prog_name
 */
void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of nodes> <data path> <starting node> <ending node>\n", prog_name);
   exit(0);
}  /* Usage */

/*--------------------------------------------------------------------
 * Function:    Dijstra
 * Purpose:     find the shortest path from start to end.
 * In arg:      arr: adjecent matrix of graph
                num_nodes: total number of nodes
                start: start point
                end: end point     
 * Out arg:     precursor: record every vertex's parent, used to find the path
                length: length of shortest path
 */
void Dijstra(int *arr, int num_nodes, int start, int end, int *precursor, int *length) {
    int found[num_nodes];
    int distance[num_nodes];

    for( int i = 0; i < num_nodes; i++ )
    {
        found[i] = -1;
        distance[i] = arr[start*num_nodes + i] == 0?infinite:arr[start*num_nodes + i];
        precursor[i] = start;
    }

    found[start] = 1;
    distance[start] = 0;

    for( int i = 0; i < num_nodes; i++ )
    {
        int min_value = infinite;
        int min_node = 0;
        for( int j = 0; j < num_nodes; j++ )
        {
            if( found[j] == -1 && distance[j] < min_value)
            {
                min_node = j;
                min_value = distance[j];
            }
        }

        found[min_node] = 1;

        if( min_node == end )
            break;

        for( int j = 0; j < num_nodes; j++ )
        {
            if( found[j] == -1 && min_value + arr[min_node * num_nodes + j] < distance[j] && min_value + arr[min_node * num_nodes + j] > 0)
            {
                distance[j] = min_value + arr[min_node * num_nodes + j];
                precursor[j] = min_node;
            }
        }
    }
    *length = distance[end];
}  /* Dijstra */


/*--------------------------------------------------------------------
 * Function:    Parallel_Dijstra
 * Purpose:     find the shortest path from start to end in just OpenMP version.
 * In arg:      arr: adjecent matrix of graph
                num_nodes: total number of nodes
                start: start point
                end: end point     
 * Out arg:     precursor: record every vertex's parent, used to find the path
                length: length of shortest path
 */
void Parallel_Dijstra(int *arr, int num_nodes, int start, int end, int *precursor, int *length) {
    int found[num_nodes];
    int distance[num_nodes];
    int i, j;
    int min_value,min_node;

    #pragma omp parallel
    #pragma omp for private(i) 
    for( i = 0; i < num_nodes; i++ )
    {
        found[i] = -1;
        distance[i] = arr[start*num_nodes + i] == 0?1000000:arr[start*num_nodes + i];
        precursor[i] = start;
    }

    found[start] = 1;
    distance[start] = 0;

    #pragma omp parallel  private(i)
    for( i = 0; i < num_nodes; i++ )
    {
        #pragma single nowait
        {
            min_value = infinite;
            min_node = 0;
        }
        
        int temp_value = infinite;
        int temp_node = 0;

        #pragma omp for private(j) nowait
        for(  j = 0; j < num_nodes; j++ )
        {
            if( distance[j] != infinite && found[j] == -1 && distance[j] < temp_value)
            {
                temp_node = j;
                temp_value = distance[j];
            }
        }
        
        if( min_value > temp_value )
        {
            #pragma omp critical
            {
                if( min_value > temp_value )
                {
                    min_value = temp_value;
                    min_node = temp_node;
                }
            }
        }
        #pragma omp barrier

        found[min_node] = 1;

        
        if( min_node == end )
            break;

        #pragma omp for private(j)
        for( j = 0; j < num_nodes; j++ )
        {
            if( found[j] == -1 && min_value + arr[min_node * num_nodes + j] < distance[j] && min_value + arr[min_node * num_nodes + j] > 0)
            {
                distance[j] = min_value + arr[min_node * num_nodes + j];
                precursor[j] = min_node;
            }
        }
    }
    *length = distance[end];
}  /*Parallel Dijstra */

/*--------------------------------------------------------------------
 * Function:    output
 * Purpose:     print out the shortest path and length.
 * In arg:      start: start point
 *              end: end point
 *              precursor: parent of every node
 *              min_length: length of shortest path    
 * Out arg:     nothing
 */
void output(int start,int end, int *precursor, int length)
{
    stack<int> seq;
    int now = end;

    while( now != start )
    {
        seq.push(now);
        now = precursor[now];
    }
    seq.push(start);

    printf("shortest path is:\n");
    int cnt = seq.size();
    for( int i = 0; i < cnt; i++ )
    {
        int temp = seq.top();
        seq.pop();
        if( seq.empty())
        {
            printf("%d", temp);
        }else{
            printf("%d->", temp);
        }
    }
    printf("\nlength of shortest path is %d\n", length);
} /* output */

