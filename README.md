# Easy reading/writing of FITS files in Julia

[![Doc][doc-dev-img]][doc-dev-url]
[![License][license-img]][license-url]
[![Build Status][github-ci-img]][github-ci-url]
[![Build Status][appveyor-img]][appveyor-url]
[![Coverage][codecov-img]][codecov-url]

This package is to facilitate reading/writing
[FITS](https://fits.gsfc.nasa.gov/fits_standard.html) files (widely used in
Astronomy) in [Julia][julia-url]. `EasyFITS` uses [Clang.jl][clang-url] to call
the functions of the [CFITSIO][cfitsio-url] library.

See [documentation][doc-dev-url].

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/EasyFITS.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/EasyFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[github-ci-img]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/emmt/EasyFITS.jl/actions/workflows/CI.yml?query=branch%3Amaster

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/EasyFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/EasyFITS-jl/branch/master

[codecov-img]: http://codecov.io/github/emmt/EasyFITS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/EasyFITS.jl?branch=master

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/

[fitsbase-url]: https://github.com/emmt/FITSIO.jl
[fitsio-url]: https://github.com/JuliaAstro/FITSIO.jl
[cfitsio-url]: https://github.com/JuliaAstro/CFITSIO.jl
[julia-url]: http://julialang.org/
[libcfitsio-url]: http://heasarc.gsfc.nasa.gov/fitsio/
[clang-url]:https://github.com/JuliaInterop/Clang.jl
