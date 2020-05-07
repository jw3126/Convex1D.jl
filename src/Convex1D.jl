module Convex1D
export minimize


let
    s = read(joinpath(dirname(@__DIR__), "README.md"), String)
    s = replace(s, "```julia" => "```jldoctest")
    @doc s Convex1D
end

function middle(x,y)
    # lets not depend on Statistics for this
    1//2 * (x + y)
end

"""
    minimize(f, (x_left, x_right); atol)

Minimize a convex function `f` on the interval `[x_left, x_right]`.

# Arguments
* f: Function that should be minimized, does not have to be differentiable.
* atol: Search tolereance, such that `isapprox(true_minimizer, sol.minimizer, atol=atol)`
will hold.
"""
function minimize(f, (x_left, x_right);
        atol=nothing
    )
    if atol === nothing
        atol = max(sqrt(eps(x_left)), sqrt(eps(x_right)))
    end
    # idea:
    # subdivide interval at five points
    # evaluate f at each point
    # recurse using x_left and x_right the neighbors of the minimizing point
    #
    # also we carefully reuse evals of f

    x1 = x_left
    x5 = x_right
    x3 = middle(x1, x5)
    x2 = middle(x1, x3)
    x4 = middle(x3, x5)
    xs = promote(x1,x2,x3,x4,x5)
    fs = map(f, xs)
    minimize_interval_subdivision_kernel(f, xs, fs, atol)
end

@noinline function minimize_interval_subdivision_kernel(f, xs, fs, atol)
    nevals = 5
    while true
        i = argmin(fs)
        at_left_bdry  = i == 1
        at_right_bdry = i == 5
        at_bdry = at_left_bdry | at_right_bdry
        if at_left_bdry
            ilo = 1; ihi = 2
        elseif at_right_bdry
            ilo = 4; ihi = 5
        else # interiour
            ilo = i-1
            ihi = i+1
        end
        @assert ilo < ihi
        x1 = xs[ilo]
        x5 = xs[ihi]
        f1 = fs[ilo]
        f5 = fs[ihi]
        if at_bdry
            x3 = middle(x1, x5)
            f3 = f(x3) ; nevals += 1
        else
            x3 = xs[i]
            f3 = fs[i]
        end

        converged = (x5 - x1) < atol
        if converged
            fs_opt = (f1,f3,f5)
            xs_opt = (x1,x3,x5)
            i_opt = argmin(fs_opt)
            return (
                minimizer = xs_opt[i_opt],
                minimum   = fs_opt[i_opt],
                nevals    = nevals,
            )
        end

        x2 = middle(x1, x3)
        f2 = f(x2); nevals += 1
        x4 = middle(x3, x5)
        f4 = f(x4); nevals += 1
        xs = (x1,x2,x3,x4,x5)
        fs = (f1,f2,f3,f4,f5)
    end
end

end # module
