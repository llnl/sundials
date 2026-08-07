// #ifndef _SUNDIALS_SYSMEMORY_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNMemoryHelper_Sys",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
  {
    auto SUNMemoryHelper_Sys_adapt_return_type_to_shared_ptr = [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
    {
      auto lambda_result = SUNMemoryHelper_Sys(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMemoryHelper>,
                             SUNMemoryHelperDeleter>(lambda_result);
    };

    return SUNMemoryHelper_Sys_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
