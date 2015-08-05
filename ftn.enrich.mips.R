ftn.enrich.mips <- function() {

######################################################################
#
# > Purpose
#   To calculate the enrichment of MIPS funcat, using 3rd level 
#   categories ('foo.3rd')
#
# > Copyright(c)
#   Author		: H. F. Lee
#   E-mail		: lee@molgen.mpg.de
#   Date Created	: 2005.02.10
#   Last Modified	: 2006.05.21
#
######################################################################

# >>> To input the working directory <<<
cat("\nEnter the full path of the directory in which you would like to run this job:\n")
foo.path <- scan(stdin(), what="character", nlines=1, quiet=T)

setwd(foo.path)

# >>> To read the list of modules <<<
cat("\nEnter the file containing the list of modules:\n")
file.name <- scan(stdin(), what="character", nlines=1, quiet=T)
modules.list <- scan(file.name, what="character")

# >>> To input the directory of modules
cat("\nEnter the full path of the directory in which modules are:\n")
foo.path.mod <- scan(stdin(), what="character", nlines=1, quiet=T)

source("./ftn.orf2gene.R")

load("../data/funcat_14032005.RData")
funcat.data[grep("^Y[A-P][LR]", funcat.data[,1]),1] <- 
	toupper(funcat.data[grep("^Y[A-P][LR]", funcat.data[,1]),1])
n_all_orfs <- length(unique(funcat.by.level[,1]))
# !!! n_all_orfs, number of all orfs in MIPS funcat without 'UNCLASSIFIED' orfs
# !!! Alternatively,
# !!! (1) n_all_orfs <- length(grep('^Y[A-P]', unique(funcat.by.level[, 1])))
# !!! (2) n_all_orfs <- length(grep('^Y[A-P]', 
# !!! 		unique(funcat.by.level[funcat.by.level[,2] != '99', 1])))
# !!!		!!! When removing 'UNCLASSIFIED' orfs

for ( i in 1:length(modules.list) )
{
	setwd(foo.path.mod)
	foo.orf <- scan(file=modules.list[i], what="character", quiet=T)
	foo.gene <- ftn.orf2gene(foo.orf)
	foo.mips <- funcat.data[which(funcat.data[,1] %in% foo.orf), ]
	if ( length(which(funcat.data[,1] %in% foo.orf)) == 1 )
	{
		foo.mips <- t(foo.mips)
	}
	foo.mips <- cbind(foo.mips, funcat.scheme[match(foo.mips[,2], funcat.scheme[,1]), 2])
	foo.mtx <- cbind(foo.mips[,1], ftn.orf2gene(foo.mips[,1]), foo.mips[,2], foo.mips[,3])
	colnames(foo.mtx) <- c("ORF", "Gene", "MIPS_code", "MIPS_ftn")
	if ( nrow(foo.mtx) != 1 )
	{
		foo.mtx <- foo.mtx[order(foo.mtx[,3]), ]
	}

	setwd(foo.path)
	if ( !file.exists('enrich.MIPS') )
	{
		dir.create("enrich.MIPS")
	}
	setwd("./enrich.MIPS/")
	write.table(file=paste(modules.list[i], ".mips", sep=""), foo.mtx, quote=F, row.names=F, sep="\t")

	foo.n_mod_orfs <- length(unique(foo.mtx[,1]))

	foo.tmp <- lapply(strsplit(foo.mtx[, 3], ".", fixed=T), function(x) x[1:3])
	# !!! 3rd level category
	foo.3rd <- sub("[.]NA$", "", sapply(foo.tmp, function(x) 
		sub("[.]NA$", "", paste(x[1], x[2], x[3], sep="."))) )

	foo.mtx[, 3] <- foo.3rd
	foo.mtx[, 4] <- funcat.scheme[match(foo.3rd, funcat.scheme[,1]), 2]

	foo.by.funcat <- split(foo.mtx[,1], foo.mtx[,3])

	if ( nrow(foo.mtx) == 1 )
	{
		foo.pval.mips <- t(foo.mtx[, c(3,4,NA)])
	}
	else
	{
		foo.pval.mips <- unique(foo.mtx[, c(3,4,NA)])
	}
	foo.pval.mips <- cbind(foo.pval.mips, 
		as.numeric(sapply(foo.by.funcat, function(x) length(x))))
	colnames(foo.pval.mips)[3:4] <- c("P_val", "N_ORFs")

	# !!! Calculation of hypergeometric p-values for all funcats in each module
	for ( j in 1:length(foo.by.funcat) )
	{
		foo.funcat_size <- as.numeric(funcat.3rd.size[names(foo.by.funcat)[j]])
		# !!! foo.funcat_size, number of all orfs in each funcat
		foo.n_funcat_orfs <- length(foo.by.funcat[[j]])
		# !!! foo.n_funcat_orfs, number of orfs with the same funcat in each module
	
		foo.pval.mips[j,3] <- formatC(1 - sum(dhyper(0:(foo.n_funcat_orfs-1), 
			foo.funcat_size, n_all_orfs-foo.funcat_size, foo.n_mod_orfs)), 
			format="e", dig=3)
	}

	if ( nrow(foo.pval.mips) != 1 )
	{
		foo.pval.mips <- foo.pval.mips[order(as.numeric(foo.pval.mips[,3])), ]
	}

	write.table(file=paste(modules.list[i], ".mips.enrich", sep=""), 
		foo.pval.mips, quote=F, row.names=F, sep="\t")

}

setwd(foo.path)

cat('\n The output files have been saved in "enrich.MIPS" directory under\n\t', foo.path, '\n')

cat("\nDone!\n\n")

}
