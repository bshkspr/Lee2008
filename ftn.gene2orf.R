ftn.gene2orf <- function(genes) {

######################################################################
#
# > Purpose
#   To convert gene names to ORF Y-names based on MIPS info
#
# > Copyright(c)
#   Author		: H. F. Lee
#   E-mail		: lee@molgen.mpg.de
#   Date Created	: 2005.02.10
#   Last Modified	: 2006.05.21
#
######################################################################

	if ( length(grep('w32', sessionInfo()$R.version$os)) != 0 )
		load("C:/research/DATA/MIPS/0_lee/all.info.RData")
	else
		load("/project/Scer/DATA/MIPS/0_lee/all.info.RData")

	foo <- all.info[match(toupper(genes), all.info[,6]), 1]
	if ( sum(toupper(genes) %in% all.info[,1]) > 0 )
		foo[toupper(genes) %in% all.info[,1]] = 
			toupper(genes)[toupper(genes) %in% all.info[,1]]
	if ( sum(is.na(foo)) > 0 ) {
		foo.na <- which(is.na(foo))
		for ( i in 1:length(foo.na) ) {
			foo.ind = grep(genes[foo.na[i]], all.info[,7], ignore.case=T)
			foo.orf = all.info[foo.ind, 7]
			if ( length(foo.orf) != 0 ) {
				foo.orf = strsplit(foo.orf, ',')
				foo.orf = which(sapply(foo.orf, function(x) 
					sum(toupper(genes)[foo.na[i]] %in% toupper(x)) > 0 ) )
				if ( length(foo.orf) == 1 ) {
					foo[foo.na[i]] <- all.info[foo.ind[foo.orf], 1]
				} else {
					stop("\nThis gene, ", genes[foo.na[i]], 
						", has two ORF names! Check it manually.\n\n")
				}
			} else {
				foo[foo.na[i]] <- genes[foo.na[i]]
			}
		}
	}
	return(foo)
}

