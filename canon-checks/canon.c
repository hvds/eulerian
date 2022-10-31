
// canon.c
// compile with gcc canon.c -o canon
// short script to count number of canonical vertex selections with respect to reflection symmetry for a hypercube

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// quicksort array of uint64_t in place (standard, forgot where I got this from)
void quicksort(uint64_t* a, int n) {
  if (n < 2) { return;}
  uint64_t p = a[n / 2];
  uint64_t* l = a;
  uint64_t* r = a + n - 1;
  while (l <= r) {
    if (*l < p) {
      l++;
    } else if (*r > p) {
      r--;
    } else {
      uint64_t t = *l;
      *l = *r;
      *r = t;
      l++;
      r--;
    }
  }
  quicksort(a, r - a + 1);
  quicksort(l, a + n - l);
}

// count distinct elements in array a of size n
int countDistinct(uint64_t *a, int n) {
  // sort array a 
  quicksort(a, n);
  // count distinct elements
  int c = 0;
  for (int i = 0; i < n; i++) {
    if (i == 0 || a[i] != a[i - 1]) {
      c++;
    }
  }
  return c;
}

// check if v is a canonical vertex selection against all (2^d)(d!) automophisms
// return -1 if it is not, otherwise return the number of distinct automorphisms
int isCanonical(uint64_t v, int d, int df, int n, int** r) {
  // malloc space for all (2^d)(d!) automorphisms 
  // to do, maybe move this out of the function so it doesn't have to be reallocated every time
  uint64_t* a = malloc(sizeof(uint64_t) * n * df);
  // for all possible bitstrings of length d
  // let the kth bit of i correspond to a reflection in the hypercube over that axis
  for (int i = 0; i < n; i++) {
    // build reflected version of v, bit by bit
    uint64_t w = 0;
    // for each bit in v
    for (int j = 0; j < n; j++) {
      int jp = j ^ i;
      // set the jth bit of w to the jpth bit of v
      w |= ((v >> jp) & 1) << j;
    }
    // check all d! permutations of the rotation array using recusive function
    for (int j = 0; j < df; j++) {
      // build permuted version of v, bit by bit
      uint64_t p = 0;
      // for each bit in v
      for (int k = 0; k < n; k++) {
        int kp = 0;
        for (int l = 0; l < d; l++) {
          // set the lth bit of kp to the r[j][l]th bit of k
          kp |= ((k >> r[j][l]) & 1) << l;
        }
        // set the kth bit of p to the kpth bit of w
        p |= ((w >> kp) & 1) << k;
      }
      // when v is the smallest out of all its reflections then it is canonical
      if (p < v) {
        free(a);
        return -1;
      }
      // add the permuted version to the array of automorphisms
      a[df * i + j] = p;
    }
  }
  // calculate number of distinct automorphisms
  int c = countDistinct(a, n * df);
  // free memory
  free(a);
  // return number of distinct automorphisms
  return c;
}

// heap algorithm to generate all permutations of an array (relatively standard)
void getPermutationTableRec(int** r, int* a, int i, int d) {
  // k is the index into the array of permutations
  static int k = 0;
  if (i == d) {
    // copy a to r
    for (int j = 0; j < d; j++) {
      r[k][j] = a[j];
    }
    k++;
    return;
  }
  for (int j = i; j < d; j++) {
    // swap a[i] and a[j]
    int t = a[i];
    a[i] = a[j];
    a[j] = t;
    // recurse
    getPermutationTableRec(r, a, i + 1, d);
    // swap a[i] and a[j]
    t = a[i];
    a[i] = a[j];
    a[j] = t;
  }
}

// setup to generate all permutations of the array [0, 1, ..., d - 1]
int** getPermutationTable(int d, int df) {
  // init memory for table
  int** r = malloc(sizeof(int*) * df);
  // init memory for rows
  for (int i = 0; i < df; i++) { r[i] = malloc(sizeof(int) * d); }
  // init initial array
  int* a = malloc(sizeof(int) * d);
  // init initial array with 0...d-1
  for (int i = 0; i < d; i++) { a[i] = i; }
  // call recursive function to fill table
  getPermutationTableRec(r, a, 0, d);
  // free a
  free(a);
  // return table
  return r;
}

// the main function takes in 1 arg, the size of the hypercube
int main(int argc, char** argv) {
    // error out if no arg provided
    if (argc != 2) {
      fprintf(stderr, "Usage: %s <size of hypercube>", argv[0]);
      exit(1);
    }
    // get the size of the hypercube
    int d = atoi(argv[1]);
    // error out if d isn't in [1, 5] since we are using 64 bit integers
    if (d < 1 || d > 5) {
      fprintf(stderr, "Error: d must be in [1, 5]");
      exit(1);
    }
    // get number of edges in hypercube
    int n = 1 << d;
    // get d!
    int df = 1; for (int i = 2; i <= d; i++) { df *= i; }
    // init store of all permutations
    int** r = getPermutationTable(d, df);
    // start counts (for error checking)
    uint64_t total = 0;
    uint64_t automorphSum = 0;
    // enumerate all possible vertex selections
    uint64_t v = (((uint64_t)1) << (n / 2)) - 1;
    while ((v >> n) == 0) {
      // check if the vertex selection is canonical
      int tmp = isCanonical(v, d, df, n, r);
      // if it is then
      if (tmp != -1) {
        // and add the number of distinct automorphisms to the sum
        automorphSum += tmp;
        // print the vertex selection and the number of distinct automorphisms
        printf("%llu %d\n", v, tmp);
      }
      // increment total count
      total += 1;
      // next permutation (bit twiddle hack)
      uint64_t t = (v | (v - 1)) + 1;
      v = t | ((((t & -t) / (v & -v)) >> 1) - 1);
    }
    // free the memory allocated for the rotations array r
    for (int i = 0; i < df; i++) { free(r[i]); }
    free(r);
    // error out if automorphSum != total
    if (automorphSum != total) {
      fprintf(stderr, "Error: automorphSum != total");
      exit(1);
    }
    // return 0
    return 0;
}