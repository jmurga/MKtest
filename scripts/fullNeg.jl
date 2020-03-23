using SpecialFunctions

function zeta(ppos::Float64,beta::Float64,al::Float64,x::Array)
	z(x,al=al,beta=beta,ppos=ppos) = (1.0-ppos)*(2.0^-al)*(beta^al)*(-SpecialFunctions.zeta(al,x+beta/2.0) + SpecialFunctions.zeta(al,(2+beta)/2.0))/((-1.0+x)*x)
	return x .|> z
end
