..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _ARKODE.ERKFullRHS:

ERKStep Full RHS
----------------

.. digraph:: erk_fullrhs

   node [shape=box]
   splines=ortho

   // -----------------
   // Before first step
   // -----------------

   init    [label="Initialize integrator\lCopy t0, y0 to tn, yn\l"]
   f0cur   [label="Set fn_is_current to False", style=filled, fillcolor=tomato1]
   h0      [label="Is h0 provided?", target="_top"]
   f0_q1   [label="Has f(tn, yn) been computed?"]
   eval1   [label="Evaluate f(tn, yn)\lStore in F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=lightskyblue]
   fcur1   [label="Set fn_is_current to True", style=filled, fillcolor=lightgreen]
   h0_comp [label="Compute h0"]
   start   [label="Start Step"]

   init    -> f0cur -> h0
   h0      -> start   [taillabel="Yes", labeldistance=2, labelangle=45]
   h0      -> f0_q1   [taillabel="No", labeldistance=2, labelangle=-45]
   f0_q1   -> h0_comp [taillabel="Yes", labeldistance=2, labelangle=-45]
   f0_q1   -> eval1   [taillabel="No", labeldistance=2, labelangle=-45]
   eval1   -> fcur1
   fcur1   -> h0_comp
   h0_comp -> start

   // ----------
   // First step
   // ----------

   f0q2     [label="Has f(tn, yn) been computed?"]
   step_q   [label="Is this the first step after initial setup?"]
   eval2    [label="Evaluate f(tn, yn)\lStore in F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=lightskyblue]
   fsal     [label="Is the method FSAL?", style=filled, fillcolor=slateblue1]
   eval3    [label="Evaluate f(tn, yn)\lStore in F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=slateblue1]
   copy1    [label="Copy F[s-1] to F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=slateblue1]
   fcur2    [label="Set fn_is_current to True", style=filled, fillcolor=lightgreen]
   stages1  [label="Compute Stages"]
   complete [label="Complete Step"]

   start   -> f0q2
   f0q2    -> stages1  [taillabel="Yes", labeldistance=2, labelangle=45]
   f0q2    -> step_q   [taillabel="No", labeldistance=2, labelangle=-45]
   step_q  -> eval2    [taillabel="Yes", labeldistance=2, labelangle=45]
   eval2   -> fcur2
   step_q  -> fsal     [taillabel="No", labeldistance=2, labelangle=-45]
   fsal    -> eval3    [taillabel="No", labeldistance=2, labelangle=-45]
   eval3   -> fcur2
   fsal    -> copy1    [taillabel="Yes", labeldistance=2, labelangle=45]
   copy1   -> fcur2
   fcur2   -> stages1 -> complete

   // -------------
   // Complete step
   // -------------

   interp_update [label="Interpolation enabled?"]
   interp_type1  [label="Using Hermite interpolation?"]
   f0q3          [label="Has f(tn, yn) been computed?"]
   update_yn     [label="Copy tcur, ycur to tn, yn"]
   fcur3         [label="Set fn_is_current to False", style=filled, fillcolor=tomato1]

   complete      -> interp_update
   interp_update -> update_yn     [taillabel="No", labeldistance=2, labelangle=45]
   interp_update -> interp_type1  [taillabel="Yes", labeldistance=2, labelangle=-45]
   interp_type1  -> update_yn     [taillabel="No", labeldistance=2, labelangle=45]
   interp_type1  -> f0q3          [taillabel="Yes", labeldistance=2, labelangle=-45]
   f0q3          -> update_yn     [taillabel="Yes", labeldistance=2, labelangle=45]
   update_yn     -> fcur3

   // ----------
   // After step
   // ----------

   interp_eval  [label="Evaluate interpolant?"]
   interp_type2 [label="Using Hermite interpolation?"]
   f1q1         [label="Has f(tn, yn) been computed?"]
   fsal1        [label="Is the method FSAL?", style=filled, fillcolor=slateblue1]
   eval4        [label="Evaluate f(tn, yn)\lStore in F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=slateblue1]
   copy2        [label="Copy F[s-1] to F[0]\lCopy F[0] to fn\l", style=filled, fillcolor=slateblue1]
   fcur4        [label="Set fn_is_current to True", style=filled, fillcolor=lightgreen]
   interp_yout  [label="Compute yout"]
   step2        [label="Start Next Step"]

   fcur3        -> interp_eval
   interp_eval  -> step2        [taillabel="No", labeldistance=2, labelangle=45]
   interp_eval  -> interp_type2 [taillabel="Yes", labeldistance=2, labelangle=-45]
   interp_type2 -> interp_yout  [taillabel="No", labeldistance=2, labelangle=45]
   interp_type2 -> f1q1         [taillabel="Yes", labeldistance=2, labelangle=-45]
   f1q1         -> fsal1        [taillabel="No", labeldistance=2, labelangle=45]
   fsal1        -> eval4        [taillabel="No", labeldistance=2, labelangle=-45]
   eval4        -> fcur4
   fsal1        -> copy2        [taillabel="Yes", labeldistance=2, labelangle=45]
   copy2        -> fcur4
   fcur4        -> interp_yout -> step2


ARKStep Full RHS
----------------

.. digraph:: ark_fullrhs_start
   :caption: ARKStep Full RHS Start

   node [shape=box, style=filled, fillcolor=white]
   splines=ortho
   bgcolor=lightskyblue

   // --------------
   // Full RHS Start
   // --------------

   fcur   [label="Has f(tn, yn) been computed?"]
   eval   [label="Evaluate fe(tn, yn), fi(tn,yn)\lStore in Fe[0], Fi[0]\l"]
   mass_a [label="Is there a mass matrix?"]
   mass_b [label="Is M time dependent?"]
   mass_c [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to f\l"]
   mass_d [label="Copy Fe[0] + Fi[0] to f\lSolve M x = f\lCopy x to f\l"]
   copy   [label="Copy Fe[0] + Fi[0] to f"]
   return [label="Set fn_is_current to True", style=filled, fillcolor=lightgreen]

   fcur   -> copy   [taillabel="Yes", labeldistance=2, labelangle=45]
   fcur   -> eval   [taillabel="No", labeldistance=2, labelangle=45]
   eval   -> mass_a
   mass_a -> copy   [taillabel="No", labeldistance=2, labelangle=45]
   copy   -> return
   mass_a -> mass_b [taillabel="Yes", labeldistance=2, labelangle=-45]
   mass_b -> mass_c [taillabel="Yes", labeldistance=2, labelangle=-45]
   mass_c -> return
   mass_b -> mass_d [taillabel="No", labeldistance=2, labelangle=45]
   mass_d -> return

.. digraph:: ark_fullrhs_start
   :caption: ARKStep Full RHS End

   node [shape=box, style=filled, fillcolor=white]
   splines=ortho
   bgcolor=slateblue1

   // ------------
   // Full RHS End
   // ------------

   fcur    [label="Has f(tn, yn) been computed?"]
   sa      [label="Is the method stiffly accurate?"]
   eval    [label="Evaluate fe(tn, yn), fi(tn, yn)\lStore in Fe[0], Fi[0]\l"]
   mass_a1 [label="Is there a mass matrix?"]
   mass_a2 [label="Is there a mass matrix?"]
   mass_b1 [label="Is M time dependent?"]
   mass_b2 [label="Is M time dependent?"]
   mass_c  [label="Solve M(t) u = Fe[0], M(t) v = Fi[0]\lStore u, v in Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to fn\l"]
   mass_d1 [label="Solve M x = fn\lCopy x to fn\l"]
   mass_d2 [label="Copy Fe[0] + Fi[0] to fn\lSolve M x = fn\lCopy x to fn\l"]
   copy_1  [label="Copy Fe[0] + Fi[0] to fn"]
   copy_2  [label="Copy Fe[s-1], Fi[s-1] to Fe[0], Fi[0]\lCopy Fe[0] + Fi[0] to fn\l"]
   return  [label="Set fn_is_current to True", style=filled, fillcolor=lightgreen]

   fcur         -> copy_1       [taillabel="Yes", labeldistance=2, labelangle=45]
   fcur         -> sa           [taillabel="No", labeldistance=2, labelangle=45]
   sa           -> copy_2       [taillabel="Yes", labeldistance=2, labelangle=45]
   copy_2       -> mass_a1
   mass_a1      -> return       [taillabel="No", labeldistance=2, labelangle=-45]
   mass_a1      -> mass_b1      [taillabel="Yes", labeldistance=2, labelangle=45]
   mass_b1      -> mass_d1      [taillabel="No", labeldistance=2, labelangle=-45]
   mass_d1      -> return
   mass_b1      -> return       [taillabel="Yes", labeldistance=2, labelangle=45]
   sa           -> eval         [taillabel="No", labeldistance=2, labelangle=-45]
   eval         -> mass_a2
   mass_a2      -> copy_1       [taillabel="No", labeldistance=2, labelangle=45]
   copy_1       -> return
   mass_a2      -> mass_b2      [taillabel="Yes", labeldistance=2, labelangle=-45]
   mass_b2      -> mass_c       [taillabel="Yes", labeldistance=2, labelangle=-45]
   mass_c       -> return
   mass_b2      -> mass_d2      [taillabel="No", labeldistance=2, labelangle=45]
   mass_d2      -> return


.. digraph:: ark_fullrhs

   node [shape=box]
   splines=ortho

   // -----------------
   // Before first step
   // -----------------

   init    [label="Initialize integrator\lCopy t0, y0 to tn, yn\l"]
   f0cur   [label="Set fn_is_current to False", style=filled, fillcolor=tomato1]
   h0      [label="Is h0 provided?"]
   f0_q    [label="Has f(tn, yn) been computed?"]
   rhs_1   [label="Call Full RHS Start", style=filled, fillcolor=lightskyblue]
   h0_comp [label="Compute h0"]
   start   [label="Start Step"]

   init -> f0cur -> h0
   h0      -> start   [taillabel="Yes", labeldistance=2, labelangle=45]
   h0      -> f0_q    [taillabel="No", labeldistance=2, labelangle=-45]
   f0_q    -> h0_comp [taillabel="Yes", labeldistance=2, labelangle=45]
   f0_q    -> rhs_1 [taillabel="No", labeldistance=2, labelangle=-45]
   rhs_1   -> h0_comp
   h0_comp -> start

   // ----------
   // Start step
   // ----------

   method_q [label="Is the first stage explicit?\nor\nIs the method stiffly accurate and Hermite interpolation is used?"]
   step_q   [label="Is this the first step after initial setup?"]
   fn_q     [label="Has f(tn, yn) been computed?"]
   rhs_2    [label="Call Full RHS Start", style=filled, fillcolor=lightskyblue]
   rhs_3    [label="Call Full RHS End", style=filled, fillcolor=slateblue1]
   stages   [label="Compute Stages"]
   complete [label="Complete Step"]

   start    -> method_q
   method_q -> stages   [taillabel="No", labeldistance=2, labelangle=45]
   method_q -> fn_q     [taillabel="Yes", labeldistance=2, labelangle=-45]
   fn_q -> stages       [taillabel="Yes", labeldistance=2, labelangle=-45]
   fn_q -> step_q       [taillabel="No", labeldistance=2, labelangle=-45]
   step_q -> rhs_2      [taillabel="Yes", labeldistance=2, labelangle=-45]
   step_q -> rhs_3      [taillabel="No", labeldistance=2, labelangle=-45]
   rhs_2 -> stages
   rhs_3 -> stages
   stages -> complete

   // -------------
   // Complete step
   // -------------

   interp_update [label="Interpolation enabled?"]
   interp_type1  [label="Using Hermite interpolation?"]
   f0q3          [label="Has f(tn, yn) been computed?"]
   rhs3          [label="Call Full RHS Start", style=filled, fillcolor=lightskyblue]
   update_yn     [label="Copy tcur, ycur to tn, yn"]
   fcur1         [label="Set fn_is_current to False", style=filled, fillcolor=tomato1]

   complete      -> interp_update
   interp_update -> update_yn     [taillabel="No", labeldistance=2, labelangle=45]
   interp_update -> interp_type1  [taillabel="Yes", labeldistance=2, labelangle=-45]
   interp_type1  -> update_yn     [taillabel="No", labeldistance=2, labelangle=45]
   interp_type1  -> f0q3          [taillabel="Yes", labeldistance=2, labelangle=-45]
   f0q3          -> update_yn     [taillabel="Yes", labeldistance=2, labelangle=45]
   f0q3          -> rhs3          [taillabel="No", labeldistance=2, labelangle=-45]
   rhs3          -> update_yn -> fcur1

   // ----------
   // After step
   // ----------

   interp_eval  [label="Evaluate interpolant?"]
   interp_type2 [label="Using Hermite interpolation?"]
   f1q1         [label="Has f(tn, yn) been computed?"]
   rhs4         [label="Call Full RHS End", style=filled, fillcolor=slateblue1]
   interp_yout  [label="Compute yout"]
   return       [label="Start Next Step"]

   fcur1        -> interp_eval
   interp_eval  -> return       [taillabel="No", labeldistance=2, labelangle=45]
   interp_eval  -> interp_type2 [taillabel="Yes", labeldistance=2, labelangle=-45]
   interp_type2 -> interp_yout  [taillabel="No", labeldistance=2, labelangle=45]
   interp_type2 -> f1q1         [taillabel="Yes", labeldistance=2, labelangle=-45]
   f1q1         -> interp_yout  [taillabel="Yes", labeldistance=2, labelangle=45]
   f1q1         -> rhs4         [taillabel="No", labeldistance=2, labelangle=45]
   rhs4         -> interp_yout -> return
