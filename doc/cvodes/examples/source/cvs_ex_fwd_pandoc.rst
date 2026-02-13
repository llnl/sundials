.. role:: math(raw)
   :format: html latex
..

.. role:: raw-latex(raw)
   :format: latex
..

.. _s:fwd_ex:

Forward sensitivity analysis example problems
=============================================

| For all the examples, any of three sensitivity method options
  (,
| , or ) can be used, and sensitivities may be
  included in the error test or not (error control set on or ,
  respectively).

The next three sections describe in detail two serial examples
( and ),
and a parallel one ().
For details on the other examples, the reader is directed to the comments in
their source files.

.. _ss:cvsAdvDiff_FSA_non:

A serial nonstiff example: cvsAdvDiff_FSA_non
---------------------------------------------

As a first example of using CVODES for forward sensitivity analysis,
we treat the simple advection-diffusion equation
for :math:`u=u(t,x)`

.. math::

   \label{e:cvsAdvDiff_FSA_non_PDE}
     \frac{\partial u}{\partial t}= q_1 \frac{\partial ^{2}u}{\partial x^{2}}
     + q_2 \frac{\partial u}{\partial x}

for :math:`0 \leq t \leq 5, ~~ 0\leq x \leq 2`, and subject to homogeneous
Dirichlet boundary conditions and initial values given by

.. math::

   \label{e:cvsAdvDiff_FSA_non_BC_IC}
     \begin{split}
       u(t,0) &= 0 \, , \quad u(t,2) = 0 \\
       u(0,x) &= x(2-x)e^{2x} \, .
     \end{split}

The nominal values of the problem parameters are :math:`q_1 = 1.0` and :math:`q_2 = 0.5`.
A system of ODEs is obtained by discretizing the :math:`x`-axis with +2
grid points and replacing the first and second order spatial derivatives
with their central difference approximations. Since the value of :math:`u` is
constant at the two endpoints, the semi-discrete equations for those points
can be eliminated. With :math:`u_{i}` as the approximation to :math:`u(t,x_{i})`,
:math:`x_{i} = i(\Delta x)`, and :math:`\Delta x=2/(\mathrm{MX}+1)`, the resulting system of
ODEs, :math:`{\dot u} = f(t,u)`, can now be written:

.. math::

   \label{e:cvsAdvDiff_FSA_non_ODE}
     {\dot u}_i= q_1 \frac{u_{i+1}-2u_{i}+u_{i-1}}{(\Delta x)^{2}}
     + q_2 \frac{u_{i+1}-u_{i-1}}{2(\Delta x)} \, .

This equation holds for :math:`i=1,2,\ldots ,` , with the understanding
that :math:`u_{0} = u_{MX+1}=0.`

The sensitivity systems for :math:`s^1 = \partial u / \partial q_1` and
:math:`s^2 = \partial u / \partial q_2` are simply

.. math::

   \label{e:cvsAdvDiff_FSA_non_S1}
     \begin{split}
       \frac{d s^1_i}{dt} 
       &= q_1 \frac{s^1_{i+1}-2s^1_{i}+s^1_{i-1}}{(\Delta x)^{2}}
       + q_2 \frac{s^1_{i+1}-s^1_{i-1}}{2(\Delta x)} 
       + \frac{u_{i+1}-2u_{i}+u_{i-1}}{(\Delta x)^{2}} \\
       s^1_i (0) &= 0.0 
     \end{split}

and

.. math::

   \label{e:cvsAdvDiff_FSA_non_S2}
     \begin{split}
       \frac{d s^2_i}{dt} 
       &= q_1 \frac{s^2_{i+1}-2s^2_{i}+s^2_{i-1}}{(\Delta x)^{2}}
       + q_2 \frac{s^2_{i+1}-s^2_{i-1}}{2(\Delta x)} 
       + \frac{u_{i+1}-u_{i-1}}{2(\Delta x)} \\
       s^1_i (0) &= 0.0  \, .
     \end{split}

This problem uses the Adams (non-stiff) integration formula and fixed-point iteration.
It is unrealistically simple [1]_,
but serves to illustrate use of the forward sensitivity capabilities in .

The file begins by including several header files, including
the main header file, the header file for the
definition of the type, and the header file
for the definitions of the serial type and operations on such vectors.
Following that are definitions of problem constants and a data block for communication
with the routine. That block includes the problem parameters and the mesh
dimension.

The program begins by processing and verifying the program arguments,
followed by allocation and initialization of the user-defined data structure. Next, the
vector of initial conditions is created (by calling ) and
initialized (in the function ). The next code block creates and allocates
memory for the object.

If sensitivity calculations were turned on through the command line arguments,
the main program continues with setting the scaling parameters
and the array of flags . In this example,
the scaling factors are used both for the finite difference approximation
to the right-hand sides of the sensitivity systems (`[e:cvsAdvDiff_FSA_non_S1] <#e:cvsAdvDiff_FSA_non_S1>`__)
and (`[e:cvsAdvDiff_FSA_non_S2] <#e:cvsAdvDiff_FSA_non_S2>`__) and in calculating the absolute tolerances for the
sensitivity variables.
The flags in are set to indicate that sensitivities with respect to both
problem parameters are desired.
The array of :math:`=2` vectors for the sensitivity variables is created
by calling and set to contain the initial values
(:math:`s^1_i(0) = 0.0`, :math:`s^2_i(0) = 0.0`).

The next three calls set optional inputs for sensitivity calculations: the sensitivity
variables are included or excluded from the error test (the boolean variable
is passed as a command line argument), the control variable is set to a value
:math:`=0` to indicate the use of second-order centered directional derivative
formulas for the approximations to the sensitivity right-hand sides, and the array of
scaling factors is passed to .
Memory for sensitivity calculations is allocated by calling
which also specifies the sensitivity solution method ( is passed
as a command line argument), and the initial conditions for the sensitivity variables.
The problem parameters and the arrays and are
passed to .

Next, in a loop over the output times, the program calls the integration
routine . On a successful return, the program prints the maximum norm
of the solution :math:`u` at the current time and, if sensitivities were also computed,
extracts and prints the maximum norms of :math:`s^1(t)` and :math:`s^2(t)`.
The program ends by printing some final integration statistics and freeing all
allocated memory.

The function is a straightforward implementation of Eqn.
(`[e:cvsAdvDiff_FSA_non_ODE] <#e:cvsAdvDiff_FSA_non_ODE>`__). The rest of the source file
contains definitions of private functions. The last
two, and , can be used with minor
modifications by any user code to print final
statistics and to check return flags from interface
functions, respectively.

Results generated by are shown in Fig. \ `[f:cvsAdvDiff_FSA_non] <#f:cvsAdvDiff_FSA_non>`__.

.. figure:: ../../../../doc/cvodes/cvsfwdnonx.pdf
   :alt: Results for the example problem.
   :width: 100%

   The time evolution of the squared solution norm, :math:`||u||^2`, is shown on the left.
   The figure on the right shows the evolution of the sensitivities of :math:`||u||^2`
   with respect to the two problem parameters.

[f:cvsAdvDiff_FSA_non]

The output generated by the example when computing
sensitivities with the method and full error
control () is as follows:

.. literalinclude:: ../../../../examples/cvodes/serial/cvsAdvDiff_FSA_non_-sensi_sim_t.out
   :language: none

.. _ss:cvsRoberts_FSA_dns:

A serial dense example: cvsRoberts_FSA_dns
------------------------------------------

This example is a modification of the chemical kinetics example
described in [cvode_ex]. It computes, in addition to the solution of the
IVP, sensitivities of the solution with respect to the three reaction rates
involved in the model. The ODEs are written as:

.. math::

   \label{e:cvsRoberts_FSA_dns_ode}
     \begin{split}
       {\dot y}_1 &= -p_1 y_1 + p_2 y_2 y_3   \\
       {\dot y}_2 &=  p_1 y_1 - p_2 y_2 y_3 - p_3 y_2^2 \\
       {\dot y}_3 &=  p_3 y_2^2 \, ,
     \end{split}

with initial conditions at :math:`t_0 = 0`, :math:`y_1 = 1` and :math:`y_2 = y_3 = 0`.
The nominal values of the reaction rate constants are
:math:`p_1 = 0.04`, :math:`p_2 = 10^4` and :math:`p_3 = 3\cdot 10^7`.
The sensitivity systems that are solved together with (`[e:cvsRoberts_FSA_dns_ode] <#e:cvsRoberts_FSA_dns_ode>`__) are

.. math::

   \label{e:cvsRoberts_FSA_dns_sens}
     \begin{split}
       & {\dot s}_i = 
       \begin{bmatrix}
         - p_1 &   p_2 y_3             &   p_2 y_2 \\
           p_1 & - p_2 y_3 - 2 p_3 y_2 & - p_2 y_2 \\
           0   &             2 p_3 y_2 &  0              
       \end{bmatrix}
       s_i + \frac{\partial f}{\partial p_i} ~,
       \quad s_i(t_0) = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}  ~,
       \quad i = 1,2,3 \\
       & \frac{\partial f}{\partial p_1} = \begin{bmatrix} -y_1 \\ y_1 \\ 0 \end{bmatrix}, \quad
       \frac{\partial f}{\partial p_2} = \begin{bmatrix} y_2 y_3 \\ -y_2 y_3 \\ 0 \end{bmatrix}, \quad
       \frac{\partial f}{\partial p_3} = \begin{bmatrix} 0 \\ - y_2^2 \\ y_2^2 \end{bmatrix} \, .
     \end{split}

The main program is described below with emphasis on the sensitivity related components.
These explanations, together with those given for the code
in [cvode_ex], will also provide the user with a template for instrumenting
an existing simulation code to perform forward sensitivity analysis.
As will be seen from this example, an existing simulation code can be modified to compute
sensitivity variables (in addition to state variables) by only inserting a few
calls into the main program.

First note that no new header files need be included. In addition to the constants already
defined in , we define the number of model parameters, (:math:`=3`),
the number of sensitivity parameters, (:math:`=3`), and a constant :math:`=0.0`.

As mentioned in , the user data structure
must provide access to the array of model parameters as the only way for
to communicate parameter values to the right-hand side function
. In the example this is done by defining
to be of type , i.e. a pointer to a structure which contains
an array of values.

Four user-supplied functions are defined. The function , passed to ,
computes the right-hand side of the ODE (`[e:cvsRoberts_FSA_dns_ode] <#e:cvsRoberts_FSA_dns_ode>`__), while
computes the dense Jacobian of the problem and is attached to the
dense linear solver module through a call to .
The function computes the right-hand side of each sensitivity system
(`[e:cvsRoberts_FSA_dns_sens] <#e:cvsRoberts_FSA_dns_sens>`__) for one parameter at a time and is therefore of type
. Finally, the function computes the error weights for the WRMS norm
estimations within .

The program prologue ends by defining six private helper functions. The first two,
and (which would not be present in a typical user code),
parse and verify the command line arguments to , respectively.
After each successful return from the main integrator, the functions
and print the state and sensitivity variables,
respectively. The function is called after completion of the
integration to print solver statistics.
The function is used to check the return flag from any of the
interface functions called by .

The program begins with definitions and type declarations.
Among these, it defines the vector of scaling factors for
the model parameters , and the array of vectors (of type )
which will contain the initial conditions and solutions for the sensitivity
variables. It also declares the variable of type
which will contain the user-defined data structure to be passed to
and used in the evaluation of the ODE right-hand sides.

The first code block in deals with reading and interpreting
the command line arguments. can be run with
or without sensitivity computations turned on and with different
selections for the sensitivity method and error control strategy.

The user’s data structure is then allocated and its field *p* is set to contain
the values of the three problem parameters.
The next block of code is identical to that in (see [cvode_ex])
and involves allocation and initialization of the state variables, and creation and
initialization of , the solver memory. It specifies that
a user-provided function () is to be used for computing the error weights.
It also attaches , with a non- Jacobian
function, as the linear solver to be used in the Newton nonlinear iteration.

If sensitivity analysis is enabled (through the command line arguments),
the main program will then set the scaling parameters
(:math:`_i` = :math:`_i`, which can typically be used for
nonzero model parameters).
Next, the program allocates memory for , by calling the function
, and initializes all sensitivity variables to :math:`0.0`.

The call to specifies the sensitivity solution
method through the argument (read from the command
line arguments) as one of , ,
or . It also specifies the user-defined routine, ,
for evaluation of the right-hand sides of sensitivity equations.

The next three calls specify optional inputs for forward sensitivity analysis:
specifying that sensitivity tolerances are to be based on ,
the error control strategy (read from the command line arguments), and
the information on the model parameters.
In this example, only is needed for the estimation of absolute sensitivity
variable tolerances; neither nor is required since the
sensitivity right-hand sides are computed in the user function . As a consequence,
we pass for the corresponding arguments in .

Note that this example uses the default estimates for the relative and absolute tolerances
and for sensitivity variables, based on the tolerances for state
variables and the scaling parameters (see for details).

Next, in a loop over the output times, the program calls the integration
routine which, if sensitivity analysis was initialized through the call
to , computes both state and sensitivity variables. However,
returns only the state solution at in the vector .
The program tests the return from for a value other than and
prints the state variables.
Sensitivity variables at are loaded into by calling .
The program tests the return from for a value other than
and then prints the sensitivity variables.

Finally, the program prints some statistics (function )
and deallocates memory through calls
to , ,
, and for the user data structure.

The user-supplied functions (for the right-hand side of the original ODEs) and
(for the system Jacobian) are identical to those in ,
with the notable exception that model parameters are extracted from the user-defined
data structure , which must first be cast to the type.
Similarly, the user-supplied function is identical to that in
. The user-supplied function computes the
sensitivity right-hand side for the -th sensitivity equation.

Results generated by are shown in
Fig. \ `[f:cvsRoberts_FSA_dns] <#f:cvsRoberts_FSA_dns>`__.

.. figure:: ../../../../doc/cvodes/cvsfwddenx.pdf
   :alt: Results for the example problem:
   :width: 100%

   Time evolution of :math:`y_1` and its sensitivities with respect to the
   three problem parameters. (Note the four different vertical scales.)

[f:cvsRoberts_FSA_dns]

The following output is generated by when computing
sensitivities with the method and full error
control ():

.. literalinclude:: ../../../../examples/cvodes/serial/cvsRoberts_FSA_dns_-sensi_sim_t.out
   :language: none

.. _ss:cvsDiurnal_FSA_kry_p:

A parallel example with user preconditioner: cvsDiurnal_FSA_kry_p
-----------------------------------------------------------------

As an example of using the forward sensitivity capabilities in
with the Krylov linear solver and the module, we describe
a test problem (derived from ) that solves the
semi-discrete form of a two-species diurnal kinetics advection-diffusion PDE
system in 2-D space, for which we also compute solution sensitivities with respect to
problem parameters (:math:`q_1` and :math:`q_2`) that appear in the kinetic rate terms.

The PDE system is

.. math::

   \label{e:cvsDiurnal_FSA_kry_p_PDE}
     \frac{\partial c^i}{\partial t} = K_h\frac{\partial^2 c^i}{\partial x^2}
     +V \frac{\partial c^i}{\partial x}
     + \frac{\partial} {\partial y} K_v(y) \frac{\partial c^i}{\partial y}
     + R^i(c^1,c^2,t) \quad (i=1,2) \, ,

where the superscripts :math:`i` are used to distinguish the two chemical
species, and where the reaction terms are given by

.. math::

   \label{e:cvsDiurnal_FSA_kry_p_R}
     \begin{split}
       R^1(c^1,c^2,t) & = -q_1c^1c^3-q_2c^1c^2+2q_3(t)c^3+q_4(t)c^2 ~, \\
       R^2(c^1,c^2,t) & = q_1c^1c^3-q_2c^1c^2-q_4(t)c^2 ~.
     \end{split}

The spatial domain is :math:`0 \leq x \leq 20,\;30 \leq y \leq 50` (in *km*).
The various constants and parameters are: :math:`K_h=4.0\cdot 10^{-6},
~ V=10^{-3},~ K_v=10^{-8}\exp (y/5),~ q_1=1.63\cdot 10^{-16},
~ q_2=4.66\cdot 10^{-16},~ c^3=3.7\cdot 10^{16},` and the diurnal
rate constants are defined as:

.. math::

   q_i(t) = 
     \left\{ \begin{array}{ll}
         \exp [-a_i/\sin \omega t], & \mbox{for } \sin \omega t>0 \\
         0, & \mbox{for } \sin \omega t\leq 0
       \end{array} \right\} ~~~(i=3,4) \, ,

where :math:`\omega =\pi /43200, ~ a_3=22.62,~ a_4=7.601.` The time interval of
integration is :math:`[0, 86400]`, representing 24 hours measured in seconds.

Homogeneous Neumann boundary conditions are imposed on each boundary, and the
initial conditions are

.. math::

   \label{e:cvsDiurnal_FSA_kry_p_IC}
     \begin{split}
     c^{1}(x,y,0) &= 10^{6}\alpha (x)\beta (y) ~,~~~ 
                       c^{2}(x,y,0)=10^{12}\alpha(x)\beta (y) ~, \\
     \alpha (x) &= 1-(0.1x-1)^{2}+(0.1x-1)^{4}/2 ~, \\
     \beta (y) &= 1-(0.1y-4)^{2}+(0.1y-4)^{4}/2 ~.
     \end{split}

We discretize the PDE system with central differencing, to
obtain an ODE system :math:`{\dot u} = f(t,u)` representing (`[e:cvsDiurnal_FSA_kry_p_PDE] <#e:cvsDiurnal_FSA_kry_p_PDE>`__).
In this case, the discrete solution vector is distributed across
many processes. Specifically, we may think of the processes as
being laid out in a rectangle, and each process being assigned a
subgrid of size :math:`\times` of the :math:`x-y` grid. If
there are processes in the :math:`x` direction and
processes in the :math:`y` direction, then the overall grid size is
:math:`\times` with :math:`=`\ :math:`\times` and
:math:`=`\ :math:`\times`, and the size of the ODE system is
:math:`2\cdot`\ :math:`\cdot`.

To compute :math:`f` in this setting, the processes pass and receive
information as follows. The solution components for the bottom row of
grid points assigned to the current process are passed to the process below
it, and the solution for the top row of grid points is received from
the process below the current process. The solution for the top
row of grid points for the current process is sent to the process
above the current process, while the solution for the bottom row of
grid points is received from that process by the current
process. Similarly, the solution for the first column of grid points
is sent from the current process to the process to its left, and
the last column of grid points is received from that process by the
current process. The communication for the solution at the right
edge of the process is similar. If this is the last process in a
particular direction, then message passing and receiving are bypassed
for that direction.

The overall structure of is very
similar to that of the code described above, with
differences arising from the use of the parallel module, .
On the other hand, the user-supplied routines in ,
for the right-hand side of the original system,
for the preconditioner setup, and for the
preconditioner solve, are identical to those defined in the example program
described in [cvode_ex]. The only difference is in the
routine , which operates on local data only and contains the actual
calculation of :math:`f(t,u)`, where the problem parameters are first extracted from the
user data structure . The program defines no
additional user-supplied routines, as it uses the internal difference quotient
routines to compute the sensitivity equation right-hand sides.

Sample results generated by are shown in
Fig. \ `[f:cvsDiurnal_FSA_kry_p] <#f:cvsDiurnal_FSA_kry_p>`__.
These results were generated on a :math:`(2\cdot40)\times(2\cdot40)` spatial grid.

.. figure:: ../../../../doc/cvodes/cvsfwdkryx_p.pdf
   :alt: Results for the example problem:
   :width: 100%

   Time evolution of :math:`c_1` and :math:`c_2` at the bottom-left and top-right corners
   (left) and of their sensitivities with respect to :math:`q_1`.

[f:cvsDiurnal_FSA_kry_p]

The following output is generated by when computing
sensitivities with the method and full error
control ():

.. literalinclude:: ../../../../examples/cvodes/parallel/cvsDiurnal_FSA_kry_p_-sensi_sim_t.out
   :language: none

.. [1]
   Increasing the number of grid points to better
   resolve the PDE spatially will lead to a stiffer ODE for which the Adams integration
   formula will not be suitable.
