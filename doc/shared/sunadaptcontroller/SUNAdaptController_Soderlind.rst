..
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
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
   ----------------------------------------------------------------

.. _SUNAdaptController.Soderlind:

The SUNAdaptController_Soderlind Module
=======================================

The Soderlind implementation of the SUNAdaptController class,
SUNAdaptController_Soderlind, implements a general structure for temporal
control proposed by G. Soderlind in :cite:p:`Sod:98`, :cite:p:`Sod:03`,
and :cite:p:`Sod:06`.  This controller has the form

.. math::
   h' = h_n \varepsilon_n^{-k_1/(p+1)} \varepsilon_{n-1}^{-k_2/(p+1)} \varepsilon_{n-2}^{-k_3/(p+1)} \left(\dfrac{h_n}{h_{n-1}}\right)^{k_4} \left(\dfrac{h_{n-1}}{h_{n-2}}\right)^{k_5}

with default parameter values :math:`k_1 = 1.25`, :math:`k_2 = 0.5`,
:math:`k_3 = -0.75`, :math:`k_4 = 0.25`, and :math:`k_5 = 0.75`, where
:math:`p` is the global order of the time integration method. If there is
insufficient history of past time steps and errors, i.e., on the first or second
time step, an I controller is used.

The SUNAdaptController_Soderlind controller is implemented as a derived
SUNAdaptController class, and defines its *content* field as:

.. code-block:: c

   struct _SUNAdaptControllerContent_Soderlind {
     sunrealtype k1;
     sunrealtype k2;
     sunrealtype k3;
     sunrealtype k4;
     sunrealtype k5;
     sunrealtype bias;
     sunrealtype ep;
     sunrealtype epp;
     sunrealtype hp;
     sunrealtype hpp;
     int firststeps;
     int historysize;
   };

These entries of the *content* field contain the following information:

* ``k1, k2, k3, k4, k5`` - controller parameters above.

* ``bias`` - error bias factor, that converts from an input temporal error
  estimate via :math:`\varepsilon = \text{bias}*\text{dsm}`.

* ``ep, epp`` - storage for the two previous error estimates,
  :math:`\varepsilon_{n-1}` and :math:`\varepsilon_{n-2}`.

* ``hp, hpp`` - storage for the previous two step sizes, :math:`h_{n-1}`
  and :math:`h_{n-2}`.

* ``firststeps`` - counter to handle first two steps (where previous
  step sizes and errors are unavailable).

* ``historysize`` - number of past step sizes or errors needed.

The header file to be included when using this module is
``sunadaptcontroller/sunadaptcontroller_soderlind.h``.

We note that through appropriate selection of the parameters :math:`k_*`,
this controller may create a wide range of proposed temporal adaptivity
controllers, including the :math:`H_{0}321`, I, PI, PID, as well as Gustafsson's
explicit and implicit controllers, :cite:p:`Gust:91` and :cite:p:`Gust:94`,
among others.  As a convenience, utility routines to create a variety of these
controllers and set their parameters (as special cases of SUNAdaptController_Soderlind) are provided.

The SUNAdaptController_Soderlind class provides implementations of all operations
relevant to a ``SUN_ADAPTCONTROLLER_H`` controller listed in
:numref:`SUNAdaptController.Description.operations`. This class
also provides the following additional user-callable routines:


.. c:function:: SUNAdaptController SUNAdaptController_Soderlind(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, and inserts its default parameters.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_Soderlind(SUNAdaptController C, sunrealtype k1, sunrealtype k2, sunrealtype k3, sunrealtype k4, sunrealtype k5)

   This user-callable function provides control over the relevant parameters
   above.  This should be called *before* the time integrator is called to evolve
   the problem.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1: parameter used within the controller time step estimate.
   :param k2: parameter used within the controller time step estimate.
   :param k3: parameter used within the controller time step estimate.
   :param k4: parameter used within the controller time step estimate.
   :param k5: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      /* Specify parameters for Soderlind's H_{0}312 controller */
      retval = SUNAdaptController_SetParams_Soderlind(C, 0.25, 0.5, 0.25, -0.75, -0.25);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_soderlind".


.. c:function:: SUNAdaptController SUNAdaptController_PID(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, set up to replicate a PID controller, and inserts its default parameters
   :math:`k_1=0.58`, :math:`k_2=-0.21`, :math:`k_3=0.1`, and :math:`k_4=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_PID(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_PID(SUNAdaptController C, sunrealtype k1, sunrealtype k2, sunrealtype k3)

   This user-callable function provides control over the relevant parameters
   above for a PID controller, setting :math:`k_4 = k_5 = 0`.  This should be
   called *before* the time integrator is called to evolve the problem.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1: parameter used within the controller time step estimate.
   :param k2: parameter used within the controller time step estimate.
   :param k3: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_PID(C, 0.58, -0.21, 0.1);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_pid".


.. c:function:: SUNAdaptController SUNAdaptController_PI(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, set up to replicate a PI controller, and inserts its default parameters
   :math:`k_1=0.8`, :math:`k_2=-0.31`, and :math:`k_3=k_4=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_PI(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_PI(SUNAdaptController C, sunrealtype k1, sunrealtype k2)

   This user-callable function provides control over the relevant parameters
   above for a PI controller, setting :math:`k_3 = k_4 = k_5 = 0`.  This should
   be called *before* the time integrator is called to evolve the problem.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1: parameter used within the controller time step estimate.
   :param k2: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_PI(C, 0.8, -0.31);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_pi".


.. c:function:: SUNAdaptController SUNAdaptController_I(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, set up to replicate an I controller, and inserts its default parameters
   :math:`k_1=1.0` and :math:`k_2=k_3=k_4=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_I(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1)

   This user-callable function provides control over the relevant parameters
   above for an I controller, setting :math:`k_2 = k_3 = k_4 = k_5 = 0`.  This
   should be called *before* the time integrator is called to evolve the problem.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1: parameter used within the controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_I(C, 1.0);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_i".


.. c:function:: SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, set up to replicate Gustafsson's explicit controller :cite:p:`Gust:91`, and
   inserts its default parameters :math:`k_1=0.635`, :math:`k_2=-0.268`, and
   :math:`k_3=k_4=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_ExpGus(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_ExpGus(SUNAdaptController C, sunrealtype k1_hat, sunrealtype k2_hat)

   This user-callable function provides control over the relevant parameters
   above for the explicit Gustafsson controller, setting :math:`k_3 = k_4 = k_5 = 0`.
   This should be called *before* the time integrator is called to evolve the problem.

   .. note::

      Gustafsson's explicit controller has the form

      .. math::
         h' = h_n \varepsilon_n^{-\hat{k}_1/(p+1)} \left(\frac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-\hat{k}_2/(p+1)}.

      The inputs to this function correspond to the values of :math:`\hat{k}_1` and :math:`\hat{k}_2`,
      which are internally transformed into the Soderlind coefficients :math:`k_1 = \hat{k}_1+\hat{k}_2`
      and :math:`k_2 = -\hat{k}_2`.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1_hat: parameter used within the explicit Gustafsson controller time step estimate.
   :param k2_hat: parameter used within the explicit Gustafsson controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_ExpGus(C, 0.367, 0.268);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_expgus".


.. c:function:: SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx)

   This constructor creates and allocates memory for a SUNAdaptController_Soderlind
   object, set up to replicate Gustafsson's implicit controller :cite:p:`Gust:94`, and
   inserts its default parameters :math:`k_1=1.93`, :math:`k_2=-0.95`, :math:`k_4=1`, and
   :math:`k_3=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :return: if successful, a usable :c:type:`SUNAdaptController` object;
            otherwise it will return ``NULL``.

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_ImpGus(sunctx);


.. c:function:: SUNErrCode SUNAdaptController_SetParams_ImpGus(SUNAdaptController C, sunrealtype k1_hat, sunrealtype k2_hat)

   This user-callable function provides control over the relevant parameters
   above for the implicit Gustafsson controller, setting :math:`k_4 = 1` and :math:`k_3 = k_5 = 0`.
   This should be called *before* the time integrator is called to evolve the problem.

   .. note::

      Gustafsson's implicit controller has the form

      .. math::
         h' = h_n \varepsilon_n^{-\hat{k}_1/(p+1)} \left(\frac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-\hat{k}_2/(p+1)} \left(\frac{h_n}{h_{n-1}}\right).

      The inputs to this function correspond to the values of :math:`\hat{k}_1` and :math:`\hat{k}_2`,
      which are internally transformed into the Soderlind coefficients :math:`k_1 = \hat{k}_1+\hat{k}_2`,
      :math:`k_2 = -\hat{k}_2`, and :math:`k_4=1`.

   :param C: the SUNAdaptController_Soderlind object.
   :param k1_hat: parameter used within the implicit Gustafsson controller time step estimate.
   :param k2_hat: parameter used within the implicit Gustafsson controller time step estimate.
   :return: :c:type:`SUNErrCode` indicating success or failure.

   Usage:

   .. code-block:: c

      retval = SUNAdaptController_SetParams_ImpGus(C, 0.98, 0.95);

   .. note::

      This routine will be called by :c:func:`SUNAdaptController_SetOptions`
      when using the key "Cid.params_impgus".

.. c:function:: SUNAdaptController SUNAdaptController_H0211(SUNContext sunctx)

   This constructor creates and allocates memory for a
   :c:type:`SUNAdaptController_Soderlind` object, set up to replicate the :math:`H_{0}211`
   controller from :cite:p:`Sod:03`, corresponding with the parameters :math:`k_1=0.5`,
   :math:`k_2=0.5`, :math:`k_4=-0.5`, and :math:`k_3=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :returns: if successful, a usable :c:type:`SUNAdaptController` object;
             otherwise it will return ``NULL``.

   .. versionadded:: 7.3.0

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_H0211(sunctx);

.. c:function:: SUNAdaptController SUNAdaptController_H0321(SUNContext sunctx)

   This constructor creates and allocates memory for a
   :c:type:`SUNAdaptController_Soderlind` object, set up to replicate the :math:`H_{0}321`
   controller from :cite:p:`Sod:03`, corresponding with the parameters :math:`k_1=1.25`,
   :math:`k_2=0.5`, :math:`k_3=-0.75`, :math:`k_4=0.25`, and :math:`k_5=0.75`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :returns: if successful, a usable :c:type:`SUNAdaptController` object;
             otherwise it will return ``NULL``.

   .. versionadded:: 7.3.0

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_H0321(sunctx);

.. c:function:: SUNAdaptController SUNAdaptController_H211(SUNContext sunctx)

   This constructor creates and allocates memory for a
   :c:type:`SUNAdaptController_Soderlind` object, set up to replicate the :math:`H211`
   controller from :cite:p:`Sod:03`, corresponding with the parameters :math:`k_1=0.25`,
   :math:`k_2=0.25`, :math:`k_4=-0.25`, and :math:`k_3=k_5=0`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :returns: if successful, a usable :c:type:`SUNAdaptController` object;
             otherwise it will return ``NULL``.

   .. versionadded:: 7.3.0

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_H211(sunctx);

.. c:function:: SUNAdaptController SUNAdaptController_H312(SUNContext sunctx)

   This constructor creates and allocates memory for a
   SUNAdaptController_Soderlind object, set up to replicate the :math:`H312`
   controller from :cite:p:`Sod:03`, corresponding with the parameters :math:`k_1=0.125`,
   :math:`k_2=0.25`, :math:`k_3=0.125`, :math:`k_4=-0.375`, and :math:`k_5=-0.125`.

   :param sunctx: the current :c:type:`SUNContext` object.
   :returns: if successful, a usable :c:type:`SUNAdaptController` object;
             otherwise it will return ``NULL``.

   .. versionadded:: 7.3.0

   Usage:

   .. code-block:: c

      SUNAdaptController C = SUNAdaptController_H312(sunctx);
