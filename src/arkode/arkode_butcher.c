/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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

#define ORDER_CONDITION_TOL (100 * SUN_UNIT_ROUNDOFF)

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

static sunbooleantype is_valid_table(const ARKodeButcherTable table)
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
    if (SUNRabs(B->b[i] - B->A[B->stages - 1][i]) > 100 * SUN_UNIT_ROUNDOFF)
    {
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

/* Grafts a branch onto a base tree while maintaining a lexicographic ordering
 * of the children */
static void butcher_product(const int* const base, const int branch,
                            int* const tree)
{
  const int base_children = base[0];
  tree[0]                 = base_children + 1;
  int i;
  for (i = 1; i <= base_children && base[i] < branch; i++)
  {
    tree[i] = base[i];
  }

  tree[i] = branch;

  for (; i <= base_children; i++) { tree[i + 1] = base[i]; }
}

/* Returns true if the trees are equal and false otherwise */
static sunbooleantype tree_equal(const int* const tree1, const int* const tree2)
{
  const int children1 = tree1[0];
  const int children2 = tree2[0];

  return children1 == children2 &&
         memcmp(&tree1[1], &tree2[1], children1 * sizeof(*tree1)) == 0;
}

typedef struct
{
  int* list;         /* A flattened array of all trees generated so far */
  int* current;      /* A memory buffer for constructing the next tree */
  int* order_offset; /* The indices into list at which each order starts */
  int length;        /* The number of ints used in list */
  int capacity;      /* The number of ints list can store */
  int order;         /* The current order */
  int root_order;    /* The current order of tree to use as a root */
  int root_offset; /* The index offset for the current root tree for the root_order */
  int branch_offset; /* The index offset for the current branch tree to graft on the root */
} tree_generator;

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

static void tree_generator_free(tree_generator* const gen)
{
  free(gen->list);
  free(gen->current);
  free(gen->order_offset);
}

static int tree_generator_push(tree_generator* const gen)
{
  for (int offset = gen->order_offset[gen->order - 1]; offset < gen->length;
       offset += gen->list[offset] + 1)
  {
    /* Check if current tree has already been generated */
    if (tree_equal(gen->current, &gen->list[offset])) { return ARK_WARNING; }
  }

  const int tree_length = gen->current[0] + 1;
  const int new_length  = gen->length + tree_length;
  if (new_length > gen->capacity)
  {
    gen->capacity = 2 * new_length;
    gen->list     = realloc(gen->list, gen->capacity * sizeof(*gen->list));
    if (gen->list == NULL) { return ARK_MEM_FAIL; }
  }

  memcpy(&gen->list[gen->length], gen->current, tree_length * sizeof(*gen->list));
  gen->length = new_length;
  return ARK_SUCCESS;
}

static void tree_print(const int* const tree, const tree_generator* const gen,
                       FILE* const outfile)
{
  const int children = tree[0];
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

static int generate_tree(tree_generator* const gen)
{
  /* Loop over orders */
  for (;;)
  {
    /* Loop over order of root tree */
    for (; gen->root_order < gen->order; gen->root_order++)
    {
      const int root_min = gen->order_offset[gen->root_order - 1];
      const int root_max = gen->order_offset[gen->root_order];

      /* Loop over trees of current root order */
      for (;;)
      {
        const int root = root_min + gen->root_offset;
        if (root == root_max) { break; }

        const int branch_order = gen->order - gen->root_order;
        const int branch_min   = gen->order_offset[branch_order - 1];
        const int branch_max   = gen->order_offset[branch_order];

        /* Loop over branches to graft to the root */
        for (;;)
        {
          const int branch = branch_min + gen->branch_offset;
          if (branch == branch_max) { break; }

          butcher_product(&gen->list[root], branch, gen->current);
          gen->branch_offset += gen->list[branch] + 1;

          const int retval = tree_generator_push(gen);
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
  sunrealtype* Phi;
  sunrealtype phi;
  sunrealtype phi_hat;
  int gamma;
  int sigma;
  int order;
} tree_props;

static void vec_set(sunrealtype* const vec, const sunrealtype value,
                    const int stages)
{
  for (int i = 0; i < stages; i++) { vec[i] = value; }
}

static void vec_times(sunrealtype* const vec1, const sunrealtype* const vec2,
                      int stages)
{
  for (int i = 0; i < stages; i++) { vec1[i] *= vec2[i]; }
}

static sunrealtype dot_prod(const sunrealtype* const vec1,
                            const sunrealtype* const vec2, const int stages)
{
  sunrealtype total = ZERO;
  for (int i = 0; i < stages; i++) { total += vec1[i] * vec2[i]; }
  return total;
}

static sunrealtype* mat_vec(sunrealtype* const* const mat,
                            const sunrealtype* const vec,
                            sunrealtype* const prod, const int stages)
{
  for (int i = 0; i < stages; i++) { prod[i] = dot_prod(mat[i], vec, stages); }
  return prod;
}

static sunbooleantype rowsum(const ARKodeButcherTable table)
{
  for (int i = 0; i < table->stages; i++)
  {
    sunrealtype rsum = SUN_RCONST(0.0);
    for (int j = 0; j < table->stages; j++) { rsum += table->A[i][j]; }
    if (SUNRabs(rsum - table->c[i]) > ORDER_CONDITION_TOL)
    {
      printf("ROWSUM ERR: %" RSYM "\n", SUNRabs(rsum - table->c[i]));
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

static tree_props get_tree_props(const int* const tree,
                                 const tree_generator* const gen, int color,
                                 const ARKodeButcherTable* const tables,
                                 sunrealtype* const buf, sunbooleantype root)
{
  tree_props props               = {.gamma = 1, .sigma = 1, .order = 1};
  const ARKodeButcherTable table = tables[color & 1];
  const int children             = tree[0];

  if (children == 0 && !root)
  {
    props.Phi = table->c;
    return props;
  }

  const int s = table->stages;
  props.Phi   = buf;
  vec_set(props.Phi, ONE, s);

  int prev_color    = -1;
  int num_duplicate = 1;
  for (int i = 1; i <= children; i++)
  {
    const int* const child       = &gen->list[tree[i]];
    const int child_color        = color >> props.order;
    sunrealtype* const child_buf = &buf[props.order * s];
    const tree_props child_props = get_tree_props(child, gen, child_color,
                                                  tables, child_buf, SUNFALSE);
    props.gamma *= child_props.gamma;
    props.sigma *= child_props.sigma;
    const int masked_child_color = child_color & ((1 << child_props.order) - 1);
    // This check relies on trees being in lexicographic order so duplicate subtrees are consecutive
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

static int compare_orders(const int given, const int computed, const int retval)
{
  if (given > computed || retval == 1) { return 1; }
  else if (given < computed || retval == -1) { return -1; }
  else { return 0; }
}

static int check_order(const ARKodeButcherTable* const tables,
                       const sunbooleantype ark, int* const q, int* const p,
                       FILE* const outfile)
{
  tree_generator gen = tree_generator_create();
  sunrealtype* buf   = NULL;
  int retval         = ARK_SUCCESS;

  while (SUNMIN(*p, *q) < 0)
  {
    retval = generate_tree(&gen);
    if (retval != ARK_SUCCESS) { break; }

    buf = realloc(buf, tables[0]->stages * gen.order * sizeof(*buf));
    if (buf == NULL)
    {
      retval = ARK_MEM_FAIL;
      break;
    }

    const int max_color = ark ? (1 << gen.order) : 1;
    for (int color = 0; color < max_color; color++)
    {
      const tree_props props = get_tree_props(gen.current, &gen, color, tables,
                                              buf, SUNTRUE);

      if (*q < 0)
      {
        const sunrealtype residual = SUNRabs(props.phi - ONE / props.gamma) /
                                     props.sigma;
        if (residual > ORDER_CONDITION_TOL)
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
        const sunrealtype embedded_residual =
          SUNRabs(props.phi_hat - ONE / props.gamma) / props.sigma;
        if (embedded_residual > ORDER_CONDITION_TOL)
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

static int check_tables(const ARKodeButcherTable* const tables,
                        const sunbooleantype ark, int* const q, int* const p,
                        FILE* const outfile)
{
  if (!is_valid_table(tables[0]) || (ark && !is_valid_table(tables[1])))
  {
    return -2;
  }

  if (outfile) { fprintf(outfile, "Order Conditions Check:\n"); }

  *q = *p = -1;
  if (rowsum(tables[0]) && (!ark || rowsum(tables[1])))
  {
    const int retval = check_order(tables, ark, q, p, outfile);
    if (retval != ARK_SUCCESS) { return -2; }
  }
  else if (outfile) { fprintf(outfile, "  method fails row sum condition\n"); }

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
  const ARKodeButcherTable tables[] = {B1, B2};
  return check_tables(tables, SUNTRUE, q, p, outfile);
}
