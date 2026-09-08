using AkaghefAlgebra
using AkaghefAlgebra.SymbolicTaft
import AkaghefAlgebra.TaftExact as TE

N=integer(:N)
a=integer(:a)
ctx=TaftContext(N)

# N stays symbolic during normalization. No basis states are enumerated.
f=finite_sum(:s,0,N-1) do s
    phase(ctx,a*s)
end
result=normalize(f,ctx)
println("Input:  ",f)
println("Output: ",result.expression)
println("Rules:  ",[step.rule for step in result.trace])
@assert result.status==:fixed_point
@assert evaluate_exact(result.expression,ctx,Dict(:N=>1000003,:a=>0);maxterms=0)==1000003
@assert evaluate_exact(result.expression,ctx,Dict(:N=>1000003,:a=>1);maxterms=0)==0
println("N=1000003: a=0 -> 1000003, a=1 -> 0; zero summands visited")

# A genuine product-coefficient contraction, with nilpotent support retained.
r=integer(:r)
contraction=finite_sum(:v,0,N-1) do v
    structure(ctx,:mul,(r,1),(1,0),(v,1))
end
contracted=normalize(contraction,ctx)
println("Contracted product coefficient: ",contracted.expression)
for rv in 0:3
    env=Dict(:N=>4,:r=>rv)
    @assert evaluate_exact(contraction,ctx,env)==evaluate_exact(contracted.expression,ctx,env)
end

# Concrete-N sparse noncommutative arithmetic in the same convention.
e,k=TE.E(4),TE.K(4)
@assert k*e==TE.zeta(4)*(e*k)
@assert iszero(e^4)
@assert e+2*e==3*e
@assert (1//2)*e+(1//2)*e==e
println("Exact N=4: KE=ζEK and E^4=0 verified")

# Deliberately unsupported simplification remains a finite sum.
gauss=finite_sum(:s,0,N-1) do s
    phase(ctx,s*s)
end
residual=normalize(gauss,ctx)
@assert residual.expression.op==:sum
println("Residual quadratic sum: ",residual.expression)
println("Explicit exact specialization at N=4: ",evaluate_exact(residual.expression,ctx,Dict(:N=>4)))
