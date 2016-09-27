# get beta mean from alpha and beta
betaMean = function(a,b)
	{
	a/(a+b)
	} 

# get beta variance from alpha and beta
betaVar = function(a,b)
	{
	(a*b)/(((a+b)^2)*(a+b+1))
	}

# get beta alpha from mean and variance
betaAlpha = function(m,SD)
	{
	(((1-m)/(SD^2))-(1/m))*m^2
	}

# get beta beta from mean and alpha
betaBeta = function(a,m)
	{
	a*((1/m)-1)
	}

