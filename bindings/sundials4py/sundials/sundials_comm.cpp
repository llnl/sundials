/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file contains the only use of the mpi4py C API. mpi4py.h defines
 * its API pointers as static per translation unit, so import_mpi4py()
 * and PyMPIComm_Get/PyMPIComm_New must stay in this file.
 * -----------------------------------------------------------------*/

#include "sundials_comm.hpp"

#include <cctype>
#include <cstring>
#include <stdexcept>
#include <string>

#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
// clang-format off
#include <mpi4py/mpi4py.h>
// clang-format on
#endif

namespace nb = nanobind;

namespace sundials4py {

namespace {

bool is_null_comm_object(nb::handle obj)
{
  if (obj.is_none()) { return true; }
  PyObject* p = obj.ptr();
  // Accept int-like objects (for example numpy integers) equal to zero.
  // bool is rejected explicitly.
  if (PyIndex_Check(p) && !PyBool_Check(p))
  {
    PyObject* idx = PyNumber_Index(p);
    if (!idx) { throw nb::python_error(); }
    nb::object idx_obj = nb::steal(idx);
    int overflow       = 0;
    long v             = PyLong_AsLongAndOverflow(idx, &overflow);
    if (v == -1 && PyErr_Occurred()) { throw nb::python_error(); }
    return v == 0 && overflow == 0;
  }
  return false;
}

std::string trim(std::string s)
{
  const char* ws = " \t\r\n";
  auto b         = s.find_first_not_of(ws);
  if (b == std::string::npos) { return std::string(); }
  auto e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

#if SUNDIALS_MPI_ENABLED

// Classifies an MPI library version string. Returns an empty string when the
// family cannot be determined (for example the MPI-ABI variants), in which
// case no hard mismatch is reported.
std::string mpi_family(const std::string& version)
{
  std::string v;
  for (char ch : version) { v.push_back(static_cast<char>(std::tolower(ch))); }
  if (v.find("open mpi") != std::string::npos ||
      v.find("open-mpi") != std::string::npos ||
      v.find("openmpi") != std::string::npos)
  {
    return "openmpi";
  }
  if (v.find("mpich") != std::string::npos ||
      v.find("mvapich") != std::string::npos ||
      v.find("intel(r) mpi") != std::string::npos)
  {
    return "mpich";
  }
  return std::string();
}

// Lazily import the mpi4py C API on first use so that importing sundials4py
// never imports mpi4py. Raises RuntimeError if mpi4py definitely uses a
// different MPI implementation family than sundials4py, and warns if the
// library version strings merely differ.
void ensure_mpi4py()
{
  static bool initialized = false;
  if (initialized) { return; }

  if (import_mpi4py() < 0) { throw nb::python_error(); }

  char buf[MPI_MAX_LIBRARY_VERSION_STRING];
  int len = 0;
  if (MPI_Get_library_version(buf, &len) == MPI_SUCCESS)
  {
    std::string ours = trim(std::string(buf, static_cast<size_t>(len)));
    nb::object theirs_obj =
      nb::module_::import_("mpi4py.MPI").attr("Get_library_version")();
    std::string theirs     = trim(nb::cast<std::string>(theirs_obj));
    std::string fam_ours   = mpi_family(ours);
    std::string fam_theirs = mpi_family(theirs);
    if (!fam_ours.empty() && !fam_theirs.empty() && fam_ours != fam_theirs)
    {
      // MPI_Comm has a different representation (int vs pointer) in the two
      // families, so using a handle across them would crash.
      throw std::runtime_error(
        "sundials4py was built against a different MPI implementation than "
        "mpi4py is using (sundials4py: '" +
        ours + "', mpi4py: '" + theirs +
        "'). Rebuild sundials4py and mpi4py against the same MPI.");
    }
    if (ours != theirs)
    {
      std::string msg =
        "sundials4py was built against a different MPI library than mpi4py "
        "is using (sundials4py: '" +
        ours + "', mpi4py: '" + theirs +
        "'). Rebuild sundials4py and mpi4py against the same MPI.";
      if (PyErr_WarnEx(PyExc_RuntimeWarning, msg.c_str(), 1) < 0)
      {
        throw nb::python_error();
      }
    }
  }

  initialized = true;
}

bool is_mpi4py_comm(nb::handle obj)
{
  nb::object comm_type = nb::module_::import_("mpi4py.MPI").attr("Comm");
  int r                = PyObject_IsInstance(obj.ptr(), comm_type.ptr());
  if (r < 0) { throw nb::python_error(); }
  return r == 1;
}

#endif

} // namespace

SUNComm to_suncomm(nb::handle obj)
{
  if (is_null_comm_object(obj)) { return SUN_COMM_NULL; }

#if SUNDIALS_MPI_ENABLED
  // If mpi4py.MPI was never imported, obj cannot be an mpi4py communicator.
  // Checking sys.modules avoids importing mpi4py for unrelated bad input.
  PyObject* mpi_mod = PyImport_GetModule(nb::str("mpi4py.MPI").ptr());
  if (mpi_mod)
  {
    Py_DECREF(mpi_mod);
    ensure_mpi4py();
    if (is_mpi4py_comm(obj))
    {
      MPI_Comm* c = PyMPIComm_Get(obj.ptr());
      if (!c) { throw nb::python_error(); }
      return *c;
    }
  }
  throw nb::type_error(
    "expected None, SUN_COMM_NULL, or an mpi4py.MPI.Comm as the communicator");
#else
  nb::object mod = nb::getattr(obj.type(), "__module__", nb::none());
  if (nb::isinstance<nb::str>(mod) &&
      std::strncmp(nb::cast<std::string>(mod).c_str(), "mpi4py", 6) == 0)
  {
    throw nb::type_error(
      "this sundials4py build was compiled without MPI support");
  }
  throw nb::type_error("expected None or SUN_COMM_NULL as the communicator");
#endif
}

nb::object from_suncomm(SUNComm c)
{
#if SUNDIALS_MPI_ENABLED
  if (c == SUN_COMM_NULL) { return nb::int_(0); }
  ensure_mpi4py();
  PyObject* o = PyMPIComm_New(c);
  if (!o) { throw nb::python_error(); }
  return nb::steal(o);
#else
  (void)c;
  return nb::int_(0);
#endif
}

} // namespace sundials4py
