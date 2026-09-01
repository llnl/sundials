// #ifndef _SUNDIALS_CUDAMEMORY_H
//
// #ifdef __cplusplus
// #endif
//

m.def(
  "SUNMemoryHelper_Cuda",
  [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
  {
    auto SUNMemoryHelper_Cuda_adapt_return_type_to_shared_ptr =
      [](SUNContext sunctx)
      -> std::shared_ptr<std::remove_pointer_t<SUNMemoryHelper>>
    {
      auto lambda_result = SUNMemoryHelper_Cuda(sunctx);

      return our_make_shared<std::remove_pointer_t<SUNMemoryHelper>,
                             SUNMemoryHelperDeleter>(lambda_result);
    };

    return SUNMemoryHelper_Cuda_adapt_return_type_to_shared_ptr(sunctx);
  },
  nb::arg("sunctx"), "nb::keep_alive<0, 1>()", nb::keep_alive<0, 1>());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
