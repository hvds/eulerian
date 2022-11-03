
// given the output of canonCount (d - 1) and SATCount (d - 2), count the total number of eulerian orientations of a hypercube of size d
// quick and dirty main function to do everything, TO DO, clean this up
// Compile: gcc finalCount.c -O3 -std=c11 -lgmp -o finalCount
// Run: ./finalCount 4 s{$d-2}.txt c{$d-1}.txt

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

// FCPair, fringe, count, and next
typedef struct t_FCPair {
  uint16_t f, c;
  struct t_FCPair *next;
} FCPair;

int main(int argc, char** argv) {
  // error out if no arg provided
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <size of hypercube> <path to SATCount's results> <path to canonCount's results>\n", argv[0]);
    exit(1);
  }
  // get the size of the hypercube
  int d = atoi(argv[1]);
  // check that size is even otherwise error out
  if (d % 2 != 0) {
    fprintf(stderr, "Error: size of hypercube must be even\n");
    exit(1);
  }
  // try to open SATCount results
  FILE* fsat = fopen(argv[2], "r");
  if (fsat == NULL) {
    fprintf(stderr, "Error: could not open SATCount's results\n");
    exit(1);
  }
  // try to open canonCount's results
  FILE* fcanon = fopen(argv[3], "r");
  if (fcanon == NULL) {
    fprintf(stderr, "Error: could not open canonCount's results\n");
    exit(1);
  }
  // malloc array of arrayLists of FCNodes to store SATCount results (of size 2^(2^(d - 2)))
  FCPair** dcr = malloc((1 << (1 << (d - 2))) * sizeof(FCPair*));
  // get initial variables
  uint16_t ov = -1;
  FCPair *tail = NULL;
  // read SATCount output line by line
  while (1) {
    // read the line into a string
    char* line = NULL;
    size_t len = 0;
    ssize_t read = getline(&line, &len, fsat);
    // if the line is empty, break
    if (read == -1) { break; }
    // otherwise parse the line into three unsigned 16 bit integers
    uint16_t v, f, c;
    sscanf(line, "%hu %hu %hu", &v, &f, &c);
    // if v == ov, add a new node to the linked list
    if (v == ov) {
      tail->next = malloc(sizeof(FCPair));
      tail = tail->next;
      tail->f = f;
      tail->c = c;
      tail->next = NULL;
    } else {
       // otherwise, create a new linked list
      ov = v;
      dcr[v] = malloc(sizeof(FCPair));
      tail = dcr[v];
      tail->f = f;
      tail->c = c;
      tail->next = NULL;
    }
  }
  // create sum variable
  mpz_t sum; mpz_init(sum); mpz_set_ui(sum, 0);
  // for every line in canonCount's results
  while (1) {
    // read the line into a string
    char* line = NULL;
    size_t len = 0;
    ssize_t read = getline(&line, &len, fcanon);
    // if the line is empty, break
    if (read == -1) { break; }
    // otherwise parse the line into an unsigned 32 bit integer and a regular int
    uint32_t v; int c;
    sscanf(line, "%u %d", &v, &c);
    // get the lower 2^(d - 2) bits of v
    uint16_t lv = v & ((1 << (1 << (d - 2))) - 1);
    // get the upper 2^(d - 2) bits of v
    uint16_t uv = v >> (1 << (d - 2));
    // get temp local sum
    mpz_t lsum; mpz_init(lsum); mpz_set_ui(lsum, 0);
    // for every node in the linked list at dcr[uv]
    FCPair* n = dcr[uv];
    while (n != NULL) {
      // for every node in the linked list at dcr[lv]
      FCPair* m = dcr[lv];
      while (m != NULL) {
        // if the fringe of n xor the fringe of m is ((1 << (1 << (d - 2))) - 1)
        if ((n->f ^ m->f) == ((1 << (1 << (d - 2))) - 1)) {
          // add the product of the counts of n and m to the sum
          mpz_t tmp;
          mpz_init(tmp);
          mpz_set_ui(tmp, n->c);
          mpz_mul_ui(tmp, tmp, m->c);
          mpz_add(lsum, lsum, tmp);
        }
        // move to the next node in the linked list at dcr[lv]
        m = m->next;
      }
      // go to the next node
      n = n->next;
    }
    // add c times the square of lsum to the sum
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(tmp, c);
    mpz_mul(tmp, tmp, lsum);
    mpz_mul(tmp, tmp, lsum);
    mpz_add(sum, sum, tmp);
  }
  // print sum
  gmp_printf("%Zd\n", sum);
}
