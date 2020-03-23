
function SfsPos(gamma::Int64,B::Float64,x::Array)
	gam = gamma*B
	sP(x,gam=gam) = 0.5*(exp(2*gam)*(1-exp(-2.0*gam*(1.0-x)))/((exp(2*gam)-1.0)*x*(1.0-x)))
	return x.|>sP
end


