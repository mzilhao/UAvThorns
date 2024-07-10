
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#define MAXBUF 1000  // max number of characters per line in input file. this
                     // should be safe

/* utility routine to check whether given array of chars is an empty line (or
   a line with only whitespace) */
static bool is_empty(const char *s)
{
  while (*s != '\0') {
    if (!isspace(*s))
      return false;
    s++;
  }
  return true;
}

int factor;

void UAv_ID_read_data(CCTK_INT *NF_p, CCTK_INT *NX_p, CCTK_REAL Xtmp[], CCTK_REAL thtmp[],
               CCTK_REAL F1[], CCTK_REAL F2[], CCTK_REAL F0[], CCTK_REAL phi0[], CCTK_REAL W[])
{
  DECLARE_CCTK_PARAMETERS;

  

  FILE *infile;
  /* open input file */
  infile = fopen(infilename, "rt");
  if (infile == NULL) {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
    "Unable to open file %s\n", infilename);
  } else {
    CCTK_VInfo(CCTK_THORNSTRING, "Reading data file %s", infilename);
  }

  // char flag[100]; // Adjust size as needed

  // // Read the flag
  // fscanf(infile, "%s\n", flag);

  // // Check the flag and act accordingly
  // if (strcmp(flag, "#4PIG=1") == 0) {
  //     factor = 0;
  // } else if (strcmp(flag, "#G=1") == 0) {
  //     factor = 1;
  // }


  /* read data from input file */
  char buf[MAXBUF];
  CCTK_INT NF = 0;     // NF will be the size of the full array
  CCTK_INT NX = 0;     // NX will be the number of X points
  bool first_block = true;
  while (fgets(buf, MAXBUF, infile) > 0) {

    // skip comments
    if (buf[0] == '#')
      continue;

    // empty lines mark a new theta coordinate. use that to count total number
    // of X points.
    if (is_empty(buf)) {
      if (first_block) {
        NX = NF;
        first_block = false;
      }
      continue;
    }

    /* printf("%s\n", buf); */
    sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf",
           &Xtmp[NF], &thtmp[NF], &F1[NF], &F2[NF], &F0[NF], &phi0[NF], &W[NF]);

    // take into account different normalization used for the stress-energy
    // tensor in the input files, which assumes 4 pi G = 1, whereas within ET it
    // is generally assumed that G = 1.
    if (norm == 1) {
    phi0[NF] *= 0.5/sqrt(M_PI);
  } else if (norm == 0) {
    phi0[NF] *= 1;
  }
    NF++;
  }

  // the following is a hack to set correctly the number of X points when we
  // only have one block (with the logic above and with just one block of data,
  // NX may never be assigned to NF after NF is properly incremented).
  if (NX == 0)
    NX = NF;

  *NF_p = NF;
  *NX_p = NX;

  fclose(infile);

  CCTK_INFO("Data file read.");

  return;
}
