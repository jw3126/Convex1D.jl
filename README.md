# Convex1d

[![Build Status](https://travis-ci.com/jw3126/Convex1d.jl.svg?branch=master)](https://travis-ci.com/jw3126/Convex1d.jl)
[![Codecov](https://codecov.io/gh/jw3126/Convex1d.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jw3126/Convex1d.jl)
[![Coveralls](https://coveralls.io/repos/github/jw3126/Convex1d.jl/badge.svg?branch=master)](https://coveralls.io/github/jw3126/Convex1d.jl?branch=master)

# Usage
```julia
using Convex1d

sol = minimize(x -> x^2, (-1, 2), atol=1e-10)
@assert isapprox(sol.minimizer, 0, atol=1e-10)
```
