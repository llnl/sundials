/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
                  Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for Butcher table structure
 * for the ARKODE infrastructure.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"

#define MAX_ORDER 10

/*---------------------------------------------------------------
  Routine to allocate an empty Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Alloc(int stages, sunbooleantype embedded)
{
  int i;
  ARKodeButcherTable B;

  /* Check for legal 'stages' value */
  if (stages < 1) { return (NULL); }

  /* Allocate Butcher table structure */
  B = NULL;
  B = (ARKodeButcherTable)malloc(sizeof(struct ARKodeButcherTableMem));
  if (B == NULL) { return (NULL); }

  /* initialize pointers in B structure to NULL */
  B->A = NULL;
  B->b = NULL;
  B->c = NULL;
  B->d = NULL;

  /* set stages into B structure */
  B->stages = stages;

  /*
   * Allocate fields within Butcher table structure
   */

  /* allocate rows of A */
  B->A = (sunrealtype**)calloc(stages, sizeof(sunrealtype*));
  if (B->A == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  /* initialize each row of A to NULL */
  for (i = 0; i < stages; i++) { B->A[i] = NULL; }

  /* allocate columns of A */
  for (i = 0; i < stages; i++)
  {
    B->A[i] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
    if (B->A[i] == NULL)
    {
      ARKodeButcherTable_Free(B);
      return (NULL);
    }
  }

  B->b = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (B->b == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  B->c = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (B->c == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  if (embedded)
  {
    B->d = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
    if (B->d == NULL)
    {
      ARKodeButcherTable_Free(B);
      return (NULL);
    }
  }

  /* initialize order parameters */
  B->q = 0;
  B->p = 0;

  return (B);
}

/*---------------------------------------------------------------
  Routine to allocate and fill a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p,
                                             sunrealtype* c, sunrealtype* A,
                                             sunrealtype* b, sunrealtype* d)
{
  int i, j;
  ARKodeButcherTable B;
  sunbooleantype embedded;

  /* Check for legal number of stages */
  if (s < 1) { return (NULL); }

  /* Does the table have an embedding? */
  embedded = (d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  B = ARKodeButcherTable_Alloc(s, embedded);
  if (B == NULL) { return (NULL); }

  /* set the relevant parameters */
  B->stages = s;
  B->q      = q;
  B->p      = p;

  for (i = 0; i < s; i++)
  {
    B->c[i] = c[i];
    B->b[i] = b[i];
    for (j = 0; j < s; j++) { B->A[i][j] = A[i * s + j]; }
  }

  if (embedded)
  {
    for (i = 0; i < s; i++) { B->d[i] = d[i]; }
  }

  return (B);
}

/*---------------------------------------------------------------
  Routine to copy a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B)
{
  int i, j, s;
  ARKodeButcherTable Bcopy;
  sunbooleantype embedded;

  /* Check for legal input */
  if (B == NULL) { return (NULL); }

  /* Get the number of stages */
  s = B->stages;

  /* Does the table have an embedding? */
  embedded = (B->d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  Bcopy = ARKodeButcherTable_Alloc(s, embedded);
  if (Bcopy == NULL) { return (NULL); }

  /* set the relevant parameters */
  Bcopy->stages = B->stages;
  Bcopy->q      = B->q;
  Bcopy->p      = B->p;

  /* Copy Butcher table */
  for (i = 0; i < s; i++)
  {
    Bcopy->c[i] = B->c[i];
    Bcopy->b[i] = B->b[i];
    for (j = 0; j < s; j++) { Bcopy->A[i][j] = B->A[i][j]; }
  }

  if (embedded)
  {
    for (i = 0; i < s; i++) { Bcopy->d[i] = B->d[i]; }
  }

  return (Bcopy);
}

/*---------------------------------------------------------------
  Routine to query the Butcher table structure workspace size
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Space(ARKodeButcherTable B, sunindextype* liw,
                              sunindextype* lrw)
{
  /* initialize outputs and return if B is not allocated */
  *liw = 0;
  *lrw = 0;
  if (B == NULL) { return; }

  /* fill outputs based on B */
  *liw = 3;
  if (B->d != NULL) { *lrw = B->stages * (B->stages + 3); }
  else { *lrw = B->stages * (B->stages + 2); }
}

/*---------------------------------------------------------------
  Routine to free a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Free(ARKodeButcherTable B)
{
  int i;

  /* Free each field within Butcher table structure, and then
     free structure itself */
  if (B != NULL)
  {
    if (B->d != NULL) { free(B->d); }
    if (B->c != NULL) { free(B->c); }
    if (B->b != NULL) { free(B->b); }
    if (B->A != NULL)
    {
      for (i = 0; i < B->stages; i++)
      {
        if (B->A[i] != NULL) { free(B->A[i]); }
      }
      free(B->A);
    }

    free(B);
  }
}

static sunbooleantype is_valid_table(ARKodeButcherTable table)
{
  if (table == NULL || table->stages < 1 || table->A == NULL ||
      table->b == NULL || table->c == NULL)
  {
    return SUNFALSE;
  }

  for (int i = 0; i < table->stages; i++)
  {
    if (table->A[i] == NULL) { return SUNFALSE; }
  }

  return SUNTRUE;
}

/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE* outfile)
{
  int i, j;

  if (!is_valid_table(B)) { return; }

  fprintf(outfile, "  A = \n");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, "      ");
    for (j = 0; j < B->stages; j++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", B->A[i][j]);
    }
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "  c = ");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", B->c[i]);
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "  b = ");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", B->b[i]);
  }
  fprintf(outfile, "\n");

  if (B->d != NULL)
  {
    fprintf(outfile, "  d = ");
    for (i = 0; i < B->stages; i++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", B->d[i]);
    }
    fprintf(outfile, "\n");
  }
}

sunbooleantype ARKodeButcherTable_IsStifflyAccurate(ARKodeButcherTable B)
{
  int i;
  for (i = 0; i < B->stages; i++)
  {
    if (SUNRabs(B->b[i] - B->A[B->stages - 1][i]) > SUN_UNIT_ROUNDOFF)
    {
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

/* Grafts a branch onto a base tree while maintaining a lexicographic ordering
 * of the children */
static void butcher_product(int* base, int branch, int* tree)
{
  int base_children = base[0];
  tree[0]           = base_children + 1;
  int i;
  for (i = 1; i <= base_children && base[i] < branch; i++)
  {
    tree[i] = base[i];
  }

  tree[i] = branch;

  for (; i <= base_children; i++) { tree[i + 1] = base[i]; }
}

/* Returns true if the trees (lexicographically ordered) are equal and false
 * otherwise. */
static sunbooleantype tree_equal(int* tree1, int* tree2)
{
  int children1 = tree1[0];
  int children2 = tree2[0];

  return children1 == children2 &&
         memcmp(&tree1[1], &tree2[1], children1 * sizeof(*tree1)) == 0;
}

typedef struct
{
  int* list;         /* A flattened array of all trees generated so far. A tree
                      * with n children is stored as n+1 integers in this list
                      * in the format {n, child_idx1, ..., child_idxn}. Children
                      * are represented by indices pointing to their position in
                      * this list. For example, the trees up to order 3 are
                      * {
                      * 1,
                      * 1, 0,
                      * 1, 1,
                      * 2, 0, 0
                      * } */
  int* current;      /* A memory buffer for constructing the next tree */
  int* order_offset; /* The indices into list at which each order starts */
  int length;        /* The number of ints used in list */
  int capacity;      /* The number of ints list can store */
  int order;         /* The current order */
  int root_order;    /* The current order of tree to use as a root */
  int root_offset;   /* The index offset for the current root tree for the
                      * root_order */
  int branch_offset; /* The index offset for the current branch tree to graft on
                      * the root */
} tree_generator;

/* A "constructor" for a tree_generator object which can produce rooted trees
 * one at a time */
static tree_generator tree_generator_create()
{
  return (tree_generator){.list          = NULL,
                          .current       = NULL,
                          .order_offset  = NULL,
                          .length        = 0,
                          .capacity      = 0,
                          .order         = 0,
                          .root_order    = 0,
                          .root_offset   = 0,
                          .branch_offset = 0};
}

/* Frees resources used by a tree_generator */
static void tree_generator_free(tree_generator* gen)
{
  free(gen->list);
  free(gen->current);
  free(gen->order_offset);
}

/* Adds gen->current to the list of trees if not already present. This may
 * require increasing the storage capacity of the generator list. */
static int tree_generator_push(tree_generator* gen)
{
  /* Loop over all trees of current order */
  for (int offset = gen->order_offset[gen->order - 1]; offset < gen->length;
       offset += gen->list[offset] + 1)
  {
    /* Check if current tree has already been generated */
    if (tree_equal(gen->current, &gen->list[offset])) { return ARK_WARNING; }
  }

  /* Check if additional capacity is needed */
  int tree_length = gen->current[0] + 1;
  int new_length  = gen->length + tree_length;
  if (new_length > gen->capacity)
  {
    gen->capacity = 2 * new_length;
    gen->list     = realloc(gen->list, gen->capacity * sizeof(*gen->list));
    if (gen->list == NULL) { return ARK_MEM_FAIL; }
  }

  /* Copy the tree from the current buffer to the list */
  memcpy(&gen->list[gen->length], gen->current, tree_length * sizeof(*gen->list));
  gen->length = new_length;
  return ARK_SUCCESS;
}

/* A utility function to write a tree in text for with t representing a leaf and
 * [...] representing joining of subtrees to a shared root */
static void tree_print(int* tree, tree_generator* gen, FILE* outfile)
{
  int children = tree[0];
  if (children == 0) { fprintf(outfile, "t"); }
  else
  {
    fprintf(outfile, "[");
    for (int i = 1; i <= children; i++)
    {
      tree_print(&gen->list[tree[i]], gen, outfile);
    }
    fprintf(outfile, "]");
  }
}

/* Generates the next rooted tree and places it into gen->current */
static int generate_tree(tree_generator* gen)
{
  /* Loop over orders */
  for (;;)
  {
    /* Loop over order of root tree */
    for (; gen->root_order < gen->order; gen->root_order++)
    {
      int root_min = gen->order_offset[gen->root_order - 1];
      int root_max = gen->order_offset[gen->root_order];

      /* Loop over trees of current root order */
      for (;;)
      {
        int root = root_min + gen->root_offset;
        if (root == root_max) { break; }

        int branch_order = gen->order - gen->root_order;
        int branch_min   = gen->order_offset[branch_order - 1];
        int branch_max   = gen->order_offset[branch_order];

        /* Loop over branches to graft to the root */
        for (;;)
        {
          int branch = branch_min + gen->branch_offset;
          if (branch == branch_max) { break; }

          butcher_product(&gen->list[root], branch, gen->current);
          gen->branch_offset += gen->list[branch] + 1;

          /* Add the tree if not already generated */
          int retval = tree_generator_push(gen);
          if (retval <= ARK_SUCCESS) { return retval; }
        }

        gen->root_offset += gen->list[root] + 1;
        gen->branch_offset = 0;
      }

      gen->root_offset = 0;
    }

    gen->root_order = 1;
    gen->order++;
    gen->current = realloc(gen->current, gen->order * sizeof(*gen->current));
    if (gen->current == NULL) { return ARK_MEM_FAIL; }
    gen->order_offset = realloc(gen->order_offset,
                                gen->order * sizeof(*gen->order_offset));
    if (gen->order_offset == NULL) { return ARK_MEM_FAIL; }
    gen->order_offset[gen->order - 1] = gen->length;

    if (gen->order == 1)
    {
      gen->current[0] = 0;
      return tree_generator_push(gen);
    }
  }
}

typedef struct
{
  sunrealtype* Phi;    /* Elementary weights for the stages */
  sunrealtype phi;     /* Elementary weight */
  sunrealtype phi_hat; /* Elementary weight for embedding */
  int gamma;           /* Density */
  int sigma;           /* Symmetry */
  int order;           /* Number of vertices */
} tree_props;

static void vec_set(sunrealtype* vec, sunrealtype value, int stages)
{
  for (int i = 0; i < stages; i++) { vec[i] = value; }
}

static void vec_times(sunrealtype* vec1, sunrealtype* vec2, int stages)
{
  for (int i = 0; i < stages; i++) { vec1[i] *= vec2[i]; }
}

static sunrealtype dot_prod(sunrealtype* vec1, sunrealtype* vec2, int stages)
{
  sunrealtype total = ZERO;
  sunrealtype err   = ZERO;
  for (int i = 0; i < stages; i++)
  {
    sunCompensatedSum(total, vec1[i] * vec2[i], &total, &err);
  }
  return total;
}

static sunrealtype* mat_vec(sunrealtype** mat, sunrealtype* vec,
                            sunrealtype* prod, int stages)
{
  for (int i = 0; i < stages; i++) { prod[i] = dot_prod(mat[i], vec, stages); }
  return prod;
}

static sunbooleantype rowsum(ARKodeButcherTable table, sunrealtype* inf_norm,
                             FILE* outfile)
{
  for (int i = 0; i < table->stages; i++)
  {
    sunrealtype row_sum     = ZERO;
    sunrealtype err         = ZERO;
    sunrealtype row_abs_sum = ZERO;
    for (int j = 0; j < table->stages; j++)
    {
      sunrealtype aij = table->A[i][j];
      sunCompensatedSum(row_sum, aij, &row_sum, &err);
      row_abs_sum += SUNRabs(aij);
    }

    if (row_abs_sum > *inf_norm) { *inf_norm = row_abs_sum; }

    /* Compensated summation has leading roundoff error of
     * 2 * epsilon * \sum_j |A_{i,j}|, so we use that to decide if row sum is
     * sufficiently close to c_i */
    sunrealtype residual = SUNRabs(row_sum - table->c[i]);
    if (residual > 2 * row_abs_sum * SUN_UNIT_ROUNDOFF)
    {
      if (outfile != NULL)
      {
        fprintf(outfile, "  row %i sum fails with residual %" RSYM "\n", i,
                residual);
      }
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

/* Recursively computes the properties of a tree with the color of each vertex
 * given by the bits of color */
static tree_props get_tree_props(int* tree, tree_generator* gen, int color,
                                 ARKodeButcherTable* tables, sunrealtype* buf,
                                 sunbooleantype root)
{
  tree_props props         = {.gamma = 1, .sigma = 1, .order = 1};
  ARKodeButcherTable table = tables[color & 1];
  int children             = tree[0];

  /* A leaf vertex corresponds to c coefficients */
  if (children == 0 && !root)
  {
    props.Phi = table->c;
    return props;
  }

  int s     = table->stages;
  props.Phi = buf;
  vec_set(props.Phi, ONE, s);

  /* Keep track of previous child color since sigma is color dependent */
  int prev_color    = -1;
  int num_duplicate = 1;
  for (int i = 1; i <= children; i++)
  {
    int* child             = &gen->list[tree[i]];
    int child_color        = color >> props.order;
    sunrealtype* child_buf = &buf[props.order * s];
    tree_props child_props = get_tree_props(child, gen, child_color, tables,
                                            child_buf, SUNFALSE);
    props.gamma *= child_props.gamma;
    props.sigma *= child_props.sigma;
    int masked_child_color = child_color & ((1 << child_props.order) - 1);
    /* This check relies on trees being in lexicographic order so duplicate
     * subtrees are consecutive */
    if (prev_color == masked_child_color && tree[i] == tree[i - 1])
    {
      num_duplicate++;
      props.sigma *= num_duplicate;
    }
    else { num_duplicate = 1; }
    props.order += child_props.order;
    prev_color = masked_child_color;
    vec_times(props.Phi, child_props.Phi, s);
  }

  props.gamma *= props.order;

  if (root)
  {
    props.phi = dot_prod(table->b, props.Phi, s);
    if (table->d != NULL) { props.phi_hat = dot_prod(table->d, props.Phi, s); }
  }
  else { props.Phi = mat_vec(table->A, props.Phi, &buf[s], s); }
  return props;
}

static int compare_orders(int given, int computed, int retval)
{
  if (given > computed || retval == 1) { return 1; }
  else if (given < computed || retval == -1) { return -1; }
  else { return 0; }
}

/* Iterates of trees and computes order condition residuals to determine the
 * order of a method and its embedding */
static int check_order(ARKodeButcherTable* tables, sunbooleantype ark, int* q,
                       int* p, sunrealtype inf_norm, FILE* outfile)
{
  tree_generator gen = tree_generator_create();
  sunrealtype* buf   = NULL;
  int retval         = ARK_SUCCESS;
  sunrealtype tol    = tables[0]->stages * inf_norm * SUN_UNIT_ROUNDOFF;

  while (*p < 0 || *q < 0)
  {
    retval = generate_tree(&gen);
    if (retval != ARK_SUCCESS) { break; }

    if (gen.order > MAX_ORDER)
    {
      if (outfile)
      {
        fprintf(outfile, "  reached maximum order of %d\n", MAX_ORDER);
      }

      if (*p < 0) { *p = MAX_ORDER; }
      if (*q < 0) { *q = MAX_ORDER; }
      break;
    }

    buf = realloc(buf, tables[0]->stages * gen.order * sizeof(*buf));
    if (buf == NULL)
    {
      retval = ARK_MEM_FAIL;
      break;
    }

    int max_color = ark ? (1 << gen.order) : 1;
    for (int color = 0; color < max_color; color++)
    {
      tree_props props = get_tree_props(gen.current, &gen, color, tables, buf,
                                        SUNTRUE);

      if (*q < 0)
      {
        sunrealtype residual = SUNRabs(props.phi - ONE / props.gamma) /
                               props.sigma;
        if (residual > tol)
        {
          *q = props.order - 1;
          if (outfile != NULL)
          {
            fprintf(outfile, "  method fails order %d condition for tree ",
                    props.order);
            tree_print(gen.current, &gen, outfile);
            fprintf(outfile, " with residual %" RSYM "\n", residual);
          }
        }
      }

      if (*p < 0)
      {
        sunrealtype embedded_residual =
          SUNRabs(props.phi_hat - ONE / props.gamma) / props.sigma;
        if (embedded_residual > tol)
        {
          *p = props.order - 1;
          if (outfile != NULL)
          {
            fprintf(outfile, "  embedding fails order %d condition for tree ",
                    props.order);
            tree_print(gen.current, &gen, outfile);
            fprintf(outfile, " with residual %" RSYM "\n", embedded_residual);
          }
        }
      }
    }
  }

  tree_generator_free(&gen);
  free(buf);
  return retval;
}

static int check_tables(ARKodeButcherTable* tables, sunbooleantype ark, int* q,
                        int* p, FILE* outfile)
{
  if (!is_valid_table(tables[0]) || (ark && !is_valid_table(tables[1])))
  {
    return -2;
  }

  if (outfile) { fprintf(outfile, "Order Conditions Check:\n"); }

  *q = *p              = -1;
  sunrealtype inf_norm = ZERO;
  if (rowsum(tables[0], &inf_norm, outfile) &&
      (!ark || rowsum(tables[1], &inf_norm, outfile)))
  {
    int retval = check_order(tables, ark, q, p, inf_norm, outfile);
    if (retval != ARK_SUCCESS) { return -2; }
  }

  int retval = 0;
  for (int i = 0; i < (ark ? 2 : 1); i++)
  {
    retval = compare_orders(tables[i]->q, *q, retval);
    retval = compare_orders(tables[i]->p, *p, retval);
  }

  return retval;
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified Butcher table.  We check the analytical [necessary]
  order conditions up through order 6.  After that, we revert to
  the [sufficient] Butcher simplifying assumptions.

  Inputs:
     B: Butcher table to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): internal {q,p} values match analytical order
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): internal p and q values are higher than analytical
         order
    -2 (failure): NULL-valued B (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, failure takes precedence over warning, which takes
  precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int* q, int* p,
                                  FILE* outfile)
{
  return check_tables(&B, SUNFALSE, q, p, outfile);
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified pair of Butcher tables in an ARK pair.  We check the
  analytical order conditions up through order 6.

  Inputs:
     B1, B2: Butcher tables to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): completed checks
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): NULL-valued B1, B2 (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, warning takes precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1, ARKodeButcherTable B2,
                                     int* q, int* p, FILE* outfile)
{
  ARKodeButcherTable tables[] = {B1, B2};
  return check_tables(tables, SUNTRUE, q, p, outfile);
}
