# Convex1D

[![Build Status](https://travis-ci.com/PTW-Freiburg/Convex1D.jl.svg?branch=master)](https://travis-ci.com/PTW-Freiburg/Convex1D.jl)
[![Codecov](https://codecov.io/gh/PTW-Freiburg/Convex1D.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PTW-Freiburg/Convex1D.jl)
[![Coveralls](https://coveralls.io/repos/github/PTW-Freiburg/Convex1D.jl/badge.svg?branch=master)](https://coveralls.io/github/PTW-Freiburg/Convex1D.jl?branch=master)

# Usage
```julia
using Convex1D

sol = minimize(x -> x^2, (-1, 2), atol=1e-10)
@assert isapprox(sol.minimizer, 0, atol=1e-10)
```
