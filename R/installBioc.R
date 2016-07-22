install.bioc = function(package)
	{
	source("https://bioconductor.org/biocLite.R")
	biocLite(package,lib=.libPaths()[1])
	}
