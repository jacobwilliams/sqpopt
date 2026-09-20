# sqpopt
Modern Fortran SQP OPTimizer

### Goals

A modern Fortran implementation of a Sequential Quadratic Programming (SQP) optimizer.

Features include:
- Modern Fortran implementation
- Modular architecture
- Open-source and actively maintained
- Sparse matrix support
- Easy integration with existing Fortran projects (uses the FPM build system)
- Selectable real kinds (single, double, quadruple)


### to Build

Use the `pixi` environment and FPM:

```
pixi shell
fpm build --profile release
fpm test --profile release
```


