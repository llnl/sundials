// #ifndef _SUNDIALS_TYPES_H
//
// #ifdef __cplusplus
// #endif
//
// #ifdef SWIG
//
// #else
//
// #endif
//

auto pyEnumSUNOutputFormat =
  nb::enum_<SUNOutputFormat>(m, "SUNOutputFormat", nb::is_arithmetic(), "")
    .value("SUN_OUTPUTFORMAT_TABLE", SUN_OUTPUTFORMAT_TABLE, "")
    .value("SUN_OUTPUTFORMAT_CSV", SUN_OUTPUTFORMAT_CSV, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//
// #ifndef SWIG
//
m.attr("SUN_COMM_NULL") = 0;
// #endif
//
// #ifdef __cplusplus
//
// #endif
//

auto pyEnumSUNDataIOMode =
  nb::enum_<SUNDataIOMode>(m, "SUNDataIOMode", nb::is_arithmetic(), "")
    .value("SUNDATAIOMODE_INMEM", SUNDATAIOMODE_INMEM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//
// #endif
