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
            return minimize_interval_subdivision_finish((x1,x3,x5), (f1,f3,f5), nevals)
        end

        x2 = middle(x1, x3)
        f2 = f(x2); nevals += 1
        x4 = middle(x3, x5)
        f4 = f(x4); nevals += 1
        xs = (x1,x2,x3,x4,x5)
        fs = (f1,f2,f3,f4,f5)
    end
end

function minimize_interval_subdivision_finish(xs::NTuple{3}, fs::NTuple{3}, nevals)
    x1,x2,x3 = xs
    f1,f2,f3 = fs
    i_opt    = argmin(fs)
    # m_left   = (f2 - f1) / (x2 - x1)
    # m_right  = (f3 - f2) / (x3 - x2)
    # y_left   = f2 - m_right * (x2 - x1)
    # y_right  = f2 + m_left  * (x3 - x2)
    # @assert y_left <= f1
    # @assert y_right <= f3
    if i_opt == 1
        xbounds = (x1, x2)
        # ybounds = (y_left, f1)
    elseif i_opt == 2
        xbounds = (x1, x3)
        # ybounds = minmax(y_left, y_right)
    elseif i_opt == 3
        xbounds = (x2, x3)
        # ybounds = (y_right, f3)
    else
        error("Unreachable")
    end
    @assert issorted(xbounds)
    # @assert issorted(ybounds)
    return (
        minimizer = xs[i_opt],
        minimum   = fs[i_opt],
        nevals    = nevals,
        xbounds   = xbounds,
        # ybounds   = ybounds,
    )

end

end # module
