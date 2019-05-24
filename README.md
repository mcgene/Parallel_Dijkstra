# Parallel version of Dijkstra

## Introduction

This code is a parallel version of Dijkstra algorithm combined with MPI and OpenMP.

## Prerequistes
* Linux
* mpicc
* OpenMP

## Run the code
```
mpirun -np [num of server] ./main [num node] [data path] [begin node] [end node]
```

If I have 3 server, my data is in the directory data, the graph contains 100001 nodes, and I want to know the shortest path from node 0 to node 8, I can run the code as

```
mpirun -np 3 ./main 10001 ./data/data.txt 0 8
```

But before running the code, you need to copy the code and data into all server, and run together.

## Build the project

If you want to build the project yourself, you can just type

```
mpi++ -fopenmp -o main src/main.cpp
```

If you want to show the speed up and run time, you can turn the comment into real code as you want.

## Data preparation

The data should be organized as

```
start1 end1 weight1
start2 end2 weight2
...
```