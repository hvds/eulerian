
// simple script to count number of ways to select edges to satisfy given constraints on in/out degree
// the main loop can be trivially parallelized using OpenMP or multiple processes
// Compile: gcc SATCount.c -O3 -std=c11 -o SATCount
// Run: ./SATCount d > s{$d}.txt

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct t_SATData {
  int d, n, e;
  int* eVal; // size e
  int** edgeToCons; // size e x 2
  int* cNum; // size n
  int* cGoal; // size n
  int** consToEdges; // size n x d
} SATData;

// for debugging
void printSATData(SATData* data) {
  printf("----\n");
  printf("d: %d n: %d e: %d\n", data->d, data->n, data->e);
  printf("eVal: ");
  for (int i = 0; i < data->e; i++) {
    printf("%d ", data->eVal[i]);
  }
  printf("\nedgeToCons: ");
  for (int i = 0; i < data->e; i++) {
    printf("(%d, %d) ", data->edgeToCons[i][0], data->edgeToCons[i][1]);
  }
  printf("\ncNum: ");
  for (int i = 0; i < data->n; i++) {
    printf("%d ", data->cNum[i]);
  }
  printf("\ncGoal: ");
  for (int i = 0; i < data->n; i++) {
    printf("%d ", data->cGoal[i]);
  }
  printf("\nconsToEdges: ");
  for (int i = 0; i < data->n; i++) {
    printf("( ");
    for (int j = 0; j < data->d; j++) {
      printf("%d ", data->consToEdges[i][j]);
    }
    printf(") ");
  }
  printf("\n");
}

// recursive procedure to get number of edge orientations for given constraints on in/out degree
uint64_t countSAT(SATData* data) {
  //printSATData(data);
  // first check if any constraints are broken
  for (int i = 0; i < data->n; i++) {
    if (data->cNum[i] < abs(data->cGoal[i])) { return 0; }
  }
  // next check if all edges have been selected
  int done = 0;
  for (int i = 0; i < data->e; i++) {
    if (data->eVal[i] == 0) { done = 1; break; }
  }
  if (done == 0) { return 1; }
  // otherwise check to see if any edges are forced
  for (int i = 0; i < data->n; i++) {
    // a constraint forces unset edges when cNum = abs(cGoal)
    // that is to say the remaining number of edges must be just enough to satisfy the constraints
    if (data->cNum[i] != 0 && data->cNum[i] == abs(data->cGoal[i])) {
      // make a copy of eVal, cGoal, and cNum using memcpy
      int* eValCopy  = malloc(data->e * sizeof(int)); memcpy(eValCopy, data->eVal, data->e * sizeof(int));
      int* cGoalCopy = malloc(data->n * sizeof(int)); memcpy(cGoalCopy, data->cGoal, data->n * sizeof(int));
      int* cNumCopy  = malloc(data->n * sizeof(int)); memcpy(cNumCopy, data->cNum, data->n * sizeof(int));
      // if cGoal is positive, then all edges must point away from the constraint
      int s = (data->cGoal[i] > 0) ? -1 : 1;
      // loop over all edges in the constraint
      for (int j = 0; j < data->d; j++) {
        // get edge index
        int e = abs(data->consToEdges[i][j]);
        // if the edge is unset
        if (data->eVal[e] == 0) {
          data->cNum[data->edgeToCons[e][0]] -= 1;
          data->cNum[data->edgeToCons[e][1]] -= 1;
          // if the edge contributes negatively to the constraint, then it must point away from the constraint
          if (data->edgeToCons[e][0] == i) {
            data->eVal[e] = s;
            data->cGoal[data->edgeToCons[e][0]] += s;
            data->cGoal[data->edgeToCons[e][1]] -= s;
          } else { // otherwise it must point towards the constraint
            data->eVal[e] = -s;
            data->cGoal[data->edgeToCons[e][0]] += -s;
            data->cGoal[data->edgeToCons[e][1]] -= -s;
          }
        }
      }
      // recurse
      uint64_t sum = countSAT(data);
      // swap data->eVal, data->cGoal, and data->cNum using pointers
      int* temp = data->eVal; data->eVal = eValCopy; eValCopy = temp;
      temp = data->cGoal; data->cGoal = cGoalCopy; cGoalCopy = temp;
      temp = data->cNum; data->cNum = cNumCopy; cNumCopy = temp;
      // free memory
      free(eValCopy);
      free(cGoalCopy);
      free(cNumCopy);
      // return the sum
      return sum;
    }
  }
  // otherwise find the first unset edge (there are more clever ways to pick an edge here)
  int i = 0;
  while (data->eVal[i] != 0) { i++; }
  // set edge to 1
  data->eVal[i] = 1;
  data->cNum[data->edgeToCons[i][0]] -= 1;
  data->cNum[data->edgeToCons[i][1]] -= 1;
  data->cGoal[data->edgeToCons[i][0]] += 1;
  data->cGoal[data->edgeToCons[i][1]] -= 1;
  // recurse
  uint64_t sum = countSAT(data);
  // set edge to -1
  data->eVal[i] = -1;
  data->cGoal[data->edgeToCons[i][0]] -= 2;
  data->cGoal[data->edgeToCons[i][1]] += 2;
  // recurse
  sum += countSAT(data);
  // restore data->eVal, data->cNum, and data->cGoal
  data->cNum[data->edgeToCons[i][0]] += 1;
  data->cNum[data->edgeToCons[i][1]] += 1;
  data->cGoal[data->edgeToCons[i][0]] += 1;
  data->cGoal[data->edgeToCons[i][1]] -= 1;
  data->eVal[i] = 0;
  // return the sum
  return sum;
}

// generate a 2d array of size num_vertices x num_vertices for consistent edge numbering
int** getEdgeIds(int n) {
  int** eId = malloc(n * sizeof(int*));
  for (int i = 0; i < n; i++) {
    eId[i] = malloc(n * sizeof(int));
  }
  int id = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      // if hamming distance between i and j is 1
      if(__builtin_popcount(i ^ j) == 1) {
        eId[i][j] = id; eId[j][i] = id;
        id++;
      } else {
        eId[i][j] = -1; eId[j][i] = -1;
      }
    }
  }
  return eId;
}

// get table of constraints
SATData* getSATData(int d, int n, int e) {
  // get edge id array
  int** eId = getEdgeIds(n);
  // allocate memory for satData
  SATData* data = malloc(sizeof(SATData));
  // set constants
  data->d = d; data->n = n; data->e = e;
  // allocate memory for eVal and set values to 0
  data->eVal = malloc(e * sizeof(int));
  for (int i = 0; i < e; i++) { data->eVal[i] = 0; }
  // allocate memory for cNum and set values for cNum to d
  data->cNum = malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) { data->cNum[i] = d; }
  // allocate memory for cGoal (will be set later)
  data->cGoal = malloc(n * sizeof(int));
  // allocate memory for edgeToCons
  data->edgeToCons = malloc(e * sizeof(int*));
  for (int i = 0; i < e; i++) { data->edgeToCons[i] = malloc(2 * sizeof(int)); }
  // allocate memory for consToEdges
  data->consToEdges = malloc(n * sizeof(int*));
  for (int i = 0; i < n; i++) { data->consToEdges[i] = malloc(d * sizeof(int)); }
  // set edgeToCons and consToEdges
  // for each vertex (constraint)
  for (int i = 0; i < n; i++) {
    // iterator
    int c = 0;
    // for each other vertex
    for (int j = 0; j < n; j++) {
      if (eId[i][j] != -1) {
        // set constraint based on edge direction
        data->consToEdges[i][c] = ((i < j) ? 1 : -1) * eId[i][j];
        // gets set twice but that's fine
        data->edgeToCons[eId[i][j]][0] = i;
        data->edgeToCons[eId[i][j]][1] = j;
        // increment iterator
        c++;
      }
    }
  }
  // free edge id
  for (int i = 0; i < n; i++) { free(eId[i]); } free(eId);
  // return satData
  return data;
}

// the main function takes in 1 arg: the size of the hypercube
int main(int argc, char** argv) {
  // error out if no args provided
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <size of hypercube>\n", argv[0]);
    exit(1);
  }
  // get the size of the hypercube
  int d = atoi(argv[1]);
  // error out if d isn't in {0, 2, 4}
  if (!(d == 0 || d == 2 || d == 4)) {
    fprintf(stderr, "Error: d must be in [0, 4]\n");
    exit(1);
  }
  // get number of vertices in hypercube
  int n = 1 << d;
  // get number of edges in hypercube
  int e = n * d / 2;
  // build constraint table
  SATData* data = getSATData(d, n, e);
  // get total number of iterations
  uint64_t itr = ((uint64_t)1) << n;
  // for all 2^n possible vertex colorings
  // TO DO, optimize so last vertex is fixed to 0 for 2x speedup
  // TO DO, canonicalization of these d-cubes to reduce number of iterations by (d!)(2^d)
  for (uint64_t i = 0; i < itr; i++) {
    // for all dangling edges
    for (uint64_t j = 0; j < itr; j++) {
      // set constraints
      for (int k = 0; k < n; k++) {
        data->cGoal[k] = ((((i >> k) & 1) == 0) ? 1 : -1) + ((((j >> k) & 1) == 0) ? 1 : -1);
      }
      // count the number of ways to select edges to satisfy given constraints on in/out degree
      uint64_t count = countSAT(data);
      // if count > 0 then print out the vertex coloring, dangling edges, and count
      if (count > 0) { printf("%llu %llu %llu\n" , i, j, count); }
    }
  }
  // free memory
  free(data->eVal);
  free(data->cNum);
  free(data->cGoal);
  for (int i = 0; i < e; i++) { free(data->edgeToCons[i]); } free(data->edgeToCons);
  for (int i = 0; i < n; i++) { free(data->consToEdges[i]); } free(data->consToEdges);
  free(data);
  // return
  return 0;
}
