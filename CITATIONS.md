# Citing SUNDIALS

We ask users to cite one or more of the following papers in their publications
reporting work done with SUNDIALS.  High-level publications describing SUNDIALS
in general include:

```bibtex
@article{roberts2026sundials,
  title     = {New {Time} {Integrators} and {Capabilities} in {SUNDIALS} {Versions} 6.2.0-7.4.0},
  author    = {Roberts, Steven B and Aggul, Mustafa and Reynolds, Daniel R and Balos, Cody J and Gardner, David J and Woodward, Carol S},
  journal   = {ACM Transactions on Mathematical Software (TOMS)},
  publisher = {ACM},
  volume    = {52},
  number    = {1},
  pages     = {8:1--8:14},
  year      = {2026},
  doi       = {10.1145/3797888}
}
```

```bibtex
@article{gardner2022sundials,
  title     = {Enabling new flexibility in the {SUNDIALS} suite of nonlinear and differential/algebraic equation solvers},
  author    = {Gardner, David J and Reynolds, Daniel R and Woodward, Carol S and Balos, Cody J},
  journal   = {ACM Transactions on Mathematical Software (TOMS)},
  publisher = {ACM},
  volume    = {48},
  number    = {3},
  pages     = {1--24},
  year      = {2022},
  doi       = {10.1145/3539801}
}
```

```bibtex
@article{hindmarsh2005sundials,
  title     = {{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
  author    = {Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
  journal   = {ACM Transactions on Mathematical Software (TOMS)},
  publisher = {ACM},
  volume    = {31},
  number    = {3},
  pages     = {363--396},
  year      = {2005},
  doi       = {10.1145/1089014.1089020}
}
```

GPU features of SUNDIALS are described in:

```bibtex
@article{balos2021enabling,
  title     = {{Enabling GPU accelerated computing in the SUNDIALS time integration library}},
  author    = {Balos, Cody J and Gardner, David J and Woodward, Carol S and Reynolds, Daniel R},
  journal   = {Parallel Computing},
  publisher = {Elsevier},
  volume    = {108},
  pages     = {102836},
  year      = {2021},
  doi       = {10.1016/j.parco.2021.102836}
}
```

The ARKODE solver was introduced in:

```bibtex
@article{reynolds2023arkode,
  title   = {{ARKODE: A flexible IVP solver infrastructure for one-step methods}},
  author  = {Reynolds, Daniel R and Gardner, David J and Woodward, Carol S and Chinomona, Rujeko},
  journal = {ACM Transactions on Mathematical Software},
  volume  = {49},
  number  = {2},
  pages   = {1--26},
  year    = {2023},
  doi     = {10.1145/3594632}
}
```

Time adaptivity for ARKODE's multirate solvers is described in:

```bibtex
@article{reynolds2026efficient,
  title     = {Efficient and {Flexible} {Multirate} {Temporal} {Adaptivity}},
  author    = {Reynolds, Daniel R. and Amihere, Sylvia and Mitchell, Dashon and Luan, Vu Thai},
  journal   = {Journal of Computational and Applied Mathematics},
  publisher = {Elsevier},
  year      = {2026},
  pages     = {117773},
  doi       = {10.1016/j.cam.2026.117773}
}
```

If desired, users can instead cite the documentation for the package and version that
they are using rather than SUNDIALS as a whole:

```bibtex
@Misc{arkodeDocumentation,
  author = {Daniel R. Reynolds and David J. Gardner and Carol S. Woodward and Cody J. Balos},
  title  = {User Documentation for ARKODE},
  year   = {2026},
  note   = {v6.7.0}
}
```

```bibtex
@Misc{cvodeDocumentation,
  author = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
  title  = {User Documentation for CVODE},
  year   = {2026},
  note   = {v7.7.0}
}
```

```bibtex
@Misc{cvodesDocumentation,
  author = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
  title  = {User Documentation for CVODES},
  year   = {2026},
  note   = {v7.7.0}
}
```

```bibtex
@Misc{idaDocumentation,
  author = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
  title  = {User Documentation for IDA},
  year   = {2026},
  note   = {v7.7.0}
}
```

```bibtex
@Misc{idasDocumentation,
  author = {Radu Serban and Cosmin Petra and Alan C. Hindmarsh and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
  title  = {User Documentation for IDAS},
  year   = {2026},
  note   = {v6.7.0}
}
```

```bibtex
@Misc{kinsolDocumentation,
  author = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
  title  = {User Documentation for KINSOL},
  year   = {2026},
  note   = {v7.7.0}
}
```
