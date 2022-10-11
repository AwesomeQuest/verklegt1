using Symbolics

"""
```julia
findErrorFromSym(symExpr:::Num; errorSuffix = "Err")
```
Find the error of an expression using the sum of the squares of it's derivatives,
"""
function findErrorFromSym(symExpr:::Num; errorSuffix = "Err")
	vars = Symbolics.get_variables(symExpr)
	varErrs = []
	for i in vars
		push!(varErrs, Symbolics.variable(string(i,errorSuffix)))
	end
	Dvars = [expand_derivatives(Differential(i)(symExpr)) for i in vars]
	symErr = sqrt(sum((Dvars[i]*varErrs[i])^2 for i in eachindex(vars)))
	return symErr
end