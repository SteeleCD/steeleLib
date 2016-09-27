# install package from bioconductor
install.bioc = function(package)
	{
	source("https://bioconductor.org/biocLite.R")
	biocLite(package,lib=.libPaths()[1],lib.loc=.libPaths()[1])
	}
