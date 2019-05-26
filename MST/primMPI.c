/**
 * Parallel implementation of Prim's MST algorithm
 *
 * Compile with: mpicc prim_mpi.c -o prim_mpi
 * Execute with: mpiexec -n 2 prim_mpi input_graph 0
 */


 /*
 *	Algorithm limitations
 *	- all vertices must be connected
 *	- the graph can't be too small(#vertices < #processors)
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>


void reduceMin(double *, double *, double *, MPI_Datatype *);
int getLocalVertex(int);

int number_of_processors; /* Number of processors */
int rank;                 /* MPI rank */
int root = 0;             /* Rank of root process (which holds the final result) */
MPI_Datatype edge;        /* Datatype containing tuple (u,v,d[u]) */
MPI_Op reduce_op;         /* Reduce operation for findimg min edge */
int number_of_vertices;   /* Number of vertices in the graph */
int vertices_per_process; /* Number of vertices that each process will get */
int extra_vertices;
int starting_vertex;      /* Initial vertice in MST */
double total_distance;       /* Total distance of MST (relevant only in root process) */

FILE * f;                 /* Input file containing graph. */

double ** weight;            /* Part of weighted adjacency matrix each process holds */
double * d;                  /* Array holding distances between MST and other vertices */
int * who;              /* who[i] holds the index of the vertex i would have to be
                             linked to in order to get a distance of d[i] */
unsigned char * in_tree;  /* 1 if the vertex is already in the MST; 0 otherwise*/

typedef struct edge_t {   /* struct representing one edge */
    double v;
    double u;
    double weight;
} edge_s;

edge_s * mst_edges;       /* Array of MST edges */
int edges_added_to_mst;   /* Number of edges currently added to MST */

double startTime;   /* Start time of computation */
double totalTime;     /* End time of computation */

//Initialise the MPI
void init_MPI(int argc, char** argv) {

	/* Parsing graph size */
	f = fopen(argv[1], "rb");
	fscanf(f, "%d", &number_of_vertices);
	//fread(&number_of_vertices, sizeof(number_of_vertices), 1, f);

	starting_vertex = atoi(argv[2]);

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Create 'edge' type (u, v, distance) */
	MPI_Type_contiguous(3, MPI_INT, &edge);
	MPI_Type_commit(&edge);

	/* Create reduce operation */
	MPI_Op_create((MPI_User_function *)reduceMin, 1, &reduce_op);
}

//Checks if the contraints to the problem are satisfacted
void check_conditions() {

	/* Check for error conditions */
	unsigned char isError = 0;
	if (number_of_vertices < number_of_processors) {
		printf("Number of vertices (%d) is smaller than number of processors (%d)\n", number_of_vertices, number_of_processors);
		isError = 1;
	}

	if (starting_vertex > number_of_vertices) {
		printf("Invalid start vertice. Must be < %d\n", number_of_vertices);
		isError = 1;
	}

	if (starting_vertex < 0) {
		printf("Invalid start vertice. Must be > -1\n");
		isError = 1;
	}

	if (isError == 1) {
		fclose(f);
		MPI_Type_free(&edge);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

}

//Read the file and initialise Data
void inputData() {

	/* Start and end index of vertex belonging to this process */
	int rank_range_start = rank * vertices_per_process;
	int rank_range_end  = (rank + 1) * vertices_per_process;
	if(rank == number_of_processors -1){
		rank_range_end += extra_vertices;
	}

	int total_edges, u, v;
	double edge_weight;
	fscanf(f, "%d", &total_edges);
	total_edges*=2;
	//fread(&total_edges, sizeof(int), 1, f);

	int i;
	double f_edge[3];
	for (i=0; i < total_edges; i++) {
		//fread(&f_edge, sizeof(double), 3, f);
		fscanf(f, "%d %d %lf", &v, &u, &edge_weight);
		//v = f_edge[0]; u = f_edge[1]; edge_weight = f_edge[2];

		if (v < 0 || u < 0) {
			printf("Invalid edge definition (%d, %d) = %d, index %d of %d edges\n", (int)v, (int)u, (int)edge_weight, i, total_edges);
			return;
		}

		if ((v >= rank_range_start) && (v < rank_range_end)) {
			weight[(int)u][getLocalVertex((int)v)] = edge_weight;
		}
		if ((u >= rank_range_start) && (u < rank_range_end)) {
			weight[(int)v][getLocalVertex((int)u)] = edge_weight;
		}
	}
	
	fclose(f);

	MPI_Barrier(MPI_COMM_WORLD);
}

//initialise Data
void initData() {

	check_conditions();

	//Variables initialisation

	/* Create data structures */
	vertices_per_process = number_of_vertices / number_of_processors;
	extra_vertices = number_of_vertices - vertices_per_process * number_of_processors;

	int i,j;

	/* Dynamic init of int weight[number_of_vertices][vertices_per_process] */
	weight = (double**)malloc(number_of_vertices * sizeof(double*));
	for (i = 0; i < number_of_vertices-1; i++) {
		weight[i] = (double*)malloc(vertices_per_process * sizeof(double));
		for (j = 0; j < vertices_per_process; j++) {
			weight[i][j] = 0;
		}
	}
	//for the extra vertices
	weight[number_of_vertices-1] = (double*)malloc((vertices_per_process + extra_vertices) * sizeof(double));
	for (j = 0; j < vertices_per_process; j++) {
		weight[number_of_vertices-1][j] = 0;
	}
	

	/* Dynamic init of d[vertices_per_process] */
	if(rank == number_of_processors-1){
		d = (double *) malloc((vertices_per_process+extra_vertices) * sizeof(double));
		who = (int *) malloc((vertices_per_process+extra_vertices) * sizeof(int));
		in_tree = (unsigned char *) malloc((vertices_per_process+extra_vertices)* sizeof(unsigned char));

		/* Initialize d with __DBL_MAX__ and inTree with 0 */
		for (i = 0; i < vertices_per_process+extra_vertices; i++) {
			d[i] = __DBL_MAX__;
			in_tree[i] = 0;
		}
	}
	
	else{
		d = (double *) malloc(vertices_per_process * sizeof(double));
		who = (int *) malloc(vertices_per_process * sizeof(int));
		in_tree = (unsigned char *) malloc(vertices_per_process * sizeof(unsigned char));

		/* Initialize d with __DBL_MAX__ and inTree with 0 */
		for (i = 0; i < vertices_per_process; i++) {
			d[i] = __DBL_MAX__;
			in_tree[i] = 0;
		}	
	}
	

	if (rank == root) {
		mst_edges = (edge_s *) malloc(number_of_vertices * sizeof(edge_s));
		edges_added_to_mst = 0;
	}

	inputData();
}

//they do the convertion between local and global indices
int getGlobalVertex(int localIndex) {
	/* Converts global vertex index to local depending on process */
	if (localIndex == -1 || rank == 0) {
		return localIndex;
	} else {
		return localIndex + (rank * vertices_per_process);
	}
}

int getLocalVertex(int globalIndex) {
	/* Converts local vertex index to global depending on process */
	int localIndex = -1;
	if ((globalIndex >= rank * vertices_per_process)) {
		if((globalIndex < (rank + 1) * vertices_per_process) ||(
            rank == number_of_processors-1 && 
            ((globalIndex < (rank + 1) * vertices_per_process)+ extra_vertices))){    
            localIndex = globalIndex - (rank * vertices_per_process);
        }
	}

	return localIndex;
}

//Refreshs d vectors
void refreshD(int addedVertex) {

	/* Updates d[] to contain minimum distances after new vertex has been added to MST */
	int i;
	for (i = 0; i < vertices_per_process; i++) {
		if ((weight[addedVertex][i] != 0) && (d[i] > weight[addedVertex][i])) {
			d[i] = weight[addedVertex][i];
			who[i] = addedVertex;
		}
	}

}

//Utilized in MPI_Reduce implementation
void reduceMin(double *invec, double *inoutvec, double *len, MPI_Datatype *type) {

	/* Find global min from MPI_Reduce buffer */
	if (invec[2] < inoutvec[2]) {
		inoutvec[0] = invec[0];
		inoutvec[1] = invec[1];
		inoutvec[2] = invec[2];
	}
}

//Implementation of prim using MPI
void find_mst() {

	/* Add starting vertex and update d[] for it */
	int localIndex = getLocalVertex(starting_vertex);
	if (localIndex != -1) {
		d[localIndex] = 0;
		in_tree[localIndex] = 1;
	}
	refreshD(starting_vertex);

	/* Buffers for MPI */
	double sendbuf[3];
	double recvbuf[3];

	/* Find MST */
	int vertices_added = 1;
	while (vertices_added < number_of_vertices) {

		/* Find the vertex with the smallest distance to the tree */
		int min = -1;
		int i;
		for (i = 0; i < vertices_per_process; i++) {
			if (in_tree[i] != 1) {
				if ((min == -1) || (d[min] > d[i])) {
					min = i;
				}
			}
		}
		if (min != -1)
		{
			sendbuf[0] = getGlobalVertex(min); // v
			sendbuf[1] = who[min];					   // u
			sendbuf[2] = d[min];					   //d[v]
		}
		else
		{
			sendbuf[0] = getGlobalVertex(min); // v
			sendbuf[1] = 0;							   // u
			sendbuf[2] = __DBL_MAX__;						//d[v]
		}

		/* Gather all results and find vertex with global minimum */
		MPI_Reduce(sendbuf, recvbuf, 3, MPI_DOUBLE, reduce_op, 0, MPI_COMM_WORLD);

		double u, v, globalMin;
		if (rank == root) {
			u = recvbuf[0];
			v = recvbuf[1];
			globalMin = recvbuf[2];

			edge_s e = { v, u, globalMin };
			mst_edges[edges_added_to_mst++] = e;
			total_distance += globalMin;
		}

		/* Broadcast vertex with global minimum */
		MPI_Bcast(&u, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

		/* Mark vertex as in tree for appropriate process */
		if (getLocalVertex((int)u) != -1) {
			in_tree[getLocalVertex((int)u)] = 1;
		}

		refreshD((int)u);

		vertices_added++;
	}
}

//Free the malloced Data and finalise the MPI
void GarbageCollector() {

	/* Free memory used by dynamic data structures */
	int i;
	for (i = 0; i < number_of_vertices; i++) {
		free(weight[i]);
	}

	free(weight);
	free(d);
	free(who);
	free(in_tree);

	/* Shut down MPI */
	MPI_Type_free(&edge);
	MPI_Finalize();

}

//Debug functions
void printWeight() {
	/* Helper function to print contents of weight[][] (used for debug) */
	int i,j;
	for (i = 0; i < number_of_vertices; i++) {
		for (j = 0; j < vertices_per_process; j++) {
			printf("Proc %d: weight[%d][%d]= %lf\n", rank, i, j, weight[i][j]);
		}
	}
}

void printDistances() {
	/* Helper function to print contents of d[] (used for debug) */
	int i;
	for (i = 0; i < vertices_per_process; i++) {
		printf("Proc %d: d[%d]= %lf\n", rank, i, d[i]);
	}
}

//Test functions

void startComputationTest() {
	startTime = MPI_Wtime();
}

void printMSTEdges() {
	printf("MST:\n");
	int i;
	for(i = 0; i < edges_added_to_mst; i++) {
		edge_s * e = &mst_edges[i];
		int v = (int)(*e).v;
		int u = (int)(*e).u;
		double w = (*e).weight;

		printf("(%d,%d)=%lf\n", v, u, w);
	}
}

void results(){
	/* Only root's has mst_edges array */
	if (rank == root) {
		totalTime = MPI_Wtime() - startTime;
		//printMSTEdges();
		printf("Total time : %lf\n", totalTime);
		printf("Total distance: %f\n", total_distance);
	}
}

//main 
int main(int argc, char** argv) {

	//Data check and inicialisation
	if (argc != 3) {
		printf("usage: %s filename starting_vertex\n", argv[0]);
		return 1;
	}
	init_MPI(argc, argv);
	initData();
	
	//Algorithm implementation
	startComputationTest();
	find_mst();
	
	//Results and return
	results();
	GarbageCollector();

	return 0;
}