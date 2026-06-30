// #ifndef _SUNDIALS_CUDAMEMORY_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNMemoryHelper_Cuda",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
  {
    auto helper = SUNMemoryHelper_Cuda(sunctx);
    return our_make_shared<std::remove_pointer_t<SUNMemoryHelper>,
                           SUNMemoryHelperDeleter>(helper);
  },
  nb::arg("sunctx"), nb::keep_alive<0, 1>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
