ftn.orf2gene <- function(orfs) {

######################################################################
#
# > Purpose
#   To convert ORF Y-names to gene names based on MIPS info
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

	foo <- all.info[match(toupper(orfs), all.info[,1]), 6]
	if ( sum(is.na(foo)) > 0 ) {
		foo.na <- which(is.na(foo))
		for ( i in 1:length(foo.na) ) {
			foo[foo.na[i]] <- toupper(orfs[foo.na[i]])
		}
	}

	return(foo)
}
