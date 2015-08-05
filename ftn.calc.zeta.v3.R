ftn.calc.zeta.v3 <- function() {

# ! To unlimit memory usage (from Claudio Lottaz),
# % limit datasize unlimited
# % limit vmemoryuse unlimited
# % limit memoryuse unlimited

######################################################################
#
# > Purpose
#   To calculate the average(|PCC|) for mRNA expression values of 
#   all ORFs in each group ( ZETA = average(|PCC|) ) and 
#   P-values for ZETAs based on Pipel et al.'s method 
#   ( version 3 of 'ftn.calc.zeta.R' )
#
# > Algorithms
#   1. Calculate the ZETA for the common target genes of TFs, say, X and Y.
#      : Zeta_XY for the gene set, Genes_XY
#   2. Calculate the ZETA for the remaining target genes of X.
#      : Zeta_X|Y for the gene set, Genes_X - Genes_XY (== Genes_X|Y)
#   3. Calculate the ZETA for the remaining target genes of Y.
#      : Zeta_Y|X for the gene set, Genes_Y - Genes_XY (== Genes_Y|X)
#   4. Decide which TF gives the maximum Zeta. Assume that it is X.
#      : TF_max = X such that Zeta_X|Y >= Zeta_Y|X
#   5. Decide if X and Y are synergistic by requiring the following condition.
#      : Delta_Zeta = Zeta_XY - Zeta_X|Y > 0
#   6. Test if Delta_Zeta is statistically significant as follows.
#   7. Divide the gene set, Genes_X, into two random groups corresponding to 
#      the two sets, Genes_XY and Genes_X|Y, respectively.
#      : Genes_XY_rnd and Genes_X|Y_rnd
#   8. Calculate the difference between the Zeta's for the two random groups.
#      : Delta_Zeta_rnd = Zeta_XY_rnd - Zeta_X|Y_rnd
#   9. Make an emperical distribution of Delta_Zeta_rnd values by dividing Genes_X
#      many times repeatedly (e.g. 1000 times).
#  10. Calculate a P-value by counting the number of Delta_Zeta_rnd values which are 
#      greater than or equal to the observed difference, Delta_Zeta.
#      : P-value = sum(Delta_Zeta_rnd >= Delta_Zeta) / simul.times
#
# > Copyright(c)
#   Author		: H. F. Lee
#   E-mail		: lee@molgen.mpg.de
#   Date Created	: 2005.07.12
#   Last Modified	: 2006.07.31
#
######################################################################


# >>> To input the working directory <<<
cat("\nEnter the full path of the directory in which you would like to run this job :\n")
foo.path <- scan(stdin(), what="character", nlines=1, quiet=T)
setwd(foo.path)

# >>> To input the name of this job <<<
cat("\nGive a name of this job as you wish :\n")
foo.name = scan(stdin(), what="character", nlines=1, quiet=T)

# >>> To ask if q-values should be reported along with p-values <<<
input.qval = function() {
	cat("\nDo you want to have Q_values?\n")
	cat("Type 'yes' or 'no'.\n")
	return(scan(stdin(), what="character", nlines=1, sep="\n", quiet=T))
}
answer.qval = input.qval()

while ( answer.qval != "yes" & answer.qval != "no" ) {
	cat("Wrong answer. Please enter again.\n")
	answer.qval = input.qval()
}

# >>> Files to write outputs to <<<
FILENAME_1 = paste(foo.name, '.result', sep='')
FILENAME_2 = paste(foo.name, '.zero.sd.orfs.txt', sep='')
FILENAME_3 = paste(foo.name, '.missing.ORFs', sep='')
FILENAME_4 = paste(foo.name, '.no.calc.modules', sep='')


# >>> To read modules data and make a list object for them <<<
cat("\nEnter the file containing modules :\n")
filename_modules = scan(stdin(), what="character", nlines=1, quiet=T)
modules = scan(filename_modules, what="character", sep="\n", quiet=T)

mod.names = sub('^> ', '', grep('^> ', modules, value=T))
mod.list  = vector('list', length = length(mod.names))
names(mod.list) = mod.names

for ( ind_mod in names(mod.list) ) {
	mod.list[[ind_mod]]$tfs = sort(unlist(strsplit(
		modules[grep(paste('> ', ind_mod, '$', sep=''), modules) + 1], ' ')))
	mod.list[[ind_mod]]$tgs = sort(unlist(strsplit(
		modules[grep(paste('> ', ind_mod, '$', sep=''), modules) + 2], ' ')))
}


# >>> To ask if all modules are tested at one time or part of them separately <<<
input.allmod = function() {
	cat("\nAre they all modules to test?\n")
	cat("Type 'yes' or 'no'.\n")
	return(scan(stdin(), what="character", nlines=1, sep="\n", quiet=T))
}
ANSWER_ALL = input.allmod()

while ( ANSWER_ALL != "yes" & ANSWER_ALL != "no" ) {
	cat("Wrong answer. Please enter again.\n")
	ANSWER_ALL = input.allmod()
}


# >>> To input the ChIP-chip data file (R data format) from which the modules are derived <<<
cat("\nEnter the R data file of the ChIP-chip matrix from which the modules are derived :\n")
file.chip.mtx <- scan(stdin(), what="character", nlines=1, quiet=T)
foo.chip.mtx <- load(file.chip.mtx)
chip.data <- get(foo.chip.mtx)
rm(list=ls(pat=foo.chip.mtx), file.chip.mtx, foo.chip.mtx)


# >>> To input the ChIP-chip p-value threshold at which the modules are derived <<<
cat("\nEnter the ChIP-chip p-value threshold at which the modules are derived :\n")
chip.pval.thrs <- scan(stdin(), nlines=1, quiet=T)


# >>> To input the expression data file (R data format) <<<
cat("\nEnter the R data file of the expression matrix :\n")
file.expr.mtx <- scan(stdin(), what="character", nlines=1, quiet=T)
foo.expr.mtx <- load(file.expr.mtx)
expr.mtx <- get(foo.expr.mtx)
rm(list=ls(pat=foo.expr.mtx), file.expr.mtx, foo.expr.mtx)

# ============================================================== #
# To retain rows of a given expression matrix which have missing
# values up to a certain proportion only (30%)

ftn.retain.mtx <- function(foo.mtx, max.na.ratio) {
  foo.mtx[apply(foo.mtx, 1, function(x) 
	length(which(is.na(x)))) <= ncol(foo.mtx)*max.na.ratio, ]

}

expr.mtx = ftn.retain.mtx(expr.mtx, 0.3)

# ************************************************************** #

# =========================================== #
# To remove ORFs with zero standard deviation

zero.sd.orfs = which(apply(expr.mtx, 1, function(x) sd(x, na.rm=T) == 0 ))
if ( length(zero.sd.orfs) != 0 ) 
	expr.mtx = expr.mtx[-zero.sd.orfs, ]

# ******************************************* #

# >>> To load the 'combinat' package <<< #
os.win = length(grep('w32', sessionInfo()$R.version$os)) != 0
if ( !('combinat' %in% names(sessionInfo()$otherPkgs)) ) {
	if ( os.win ) library(combinat) else
		library(combinat, lib.loc='/home/lee/programming/R.BioC/lib/')
}

# >>> Preparations <<< #
N_mod = length(mod.list)
cat("\n\t The number of modules :", N_mod, "\n\n")

# !!! The number of random samples
cat("\n The number of random samples : ")
N_RANDOM = scan(stdin(), quiet=T, nlines=1)
cat("\n")

time.start <- date()
cat("Program Start Time :\n")
cat("\t", time.start, "\n\n")

cat("====================\n\n")

# ============================================= #
# To make a unique filename using date and time

foo.date <- unlist(strsplit(as.character(Sys.time()), " "))
foo.date <- paste(gsub("-", "", foo.date[1]), ".", 
	gsub(":", "", foo.date[2]), sep="")

LOG_FILE = paste(foo.name, foo.date, 'log', sep='.')

# ********************************************* #

cat(file=LOG_FILE, "\n> Program Goal :\n")
cat(file=LOG_FILE, "This program is to calculate the average of absolute Pearson coefficients for ", 
	"all pairs of mRNA expression profiles of all ORFs in each group ", 
	"( ZETA = average(|PCC|) ) and a P-value/Q-value for the difference of the ZETA ", 
	"from the maximum ZETA for the remaining target sets of individual TFs ", 
	"(based on the Pilpel et al. method).\n", sep='', append=T)
cat(file=LOG_FILE, "\n> The chosen directory to save output files :\n", 
	foo.path, "\n", sep='', append=T)
cat(file=LOG_FILE, "\n> The chosen name for this job :\n", 
	foo.name, "\n", sep='', append=T)
if ( answer.qval == 'yes' ) {
	cat(file=LOG_FILE, "\n> You have q values reported.\n", append=T)
}
if ( length(zero.sd.orfs) != 0 ) {
	writeLines(rownames(expr.mtx)[zero.sd.orfs], con=FILENAME_2)
	cat(file=LOG_FILE, "\n> In the expression matrix provided, there are ", length(zero.sd.orfs), 
		" ORFs with zero standard deviation over all conditions.\n", 
		"They were removed from the matrix, and the list of those ORFs ",
		"was saved in the following file :\n", FILENAME_2, "\n", sep='', append=T)
}


#####################################################################
##### 			>>>>> MAIN CODE <<<<< 			#####
#####################################################################

all.orfs <- rownames(expr.mtx)
pcc.mtx <- cor(t(expr.mtx), use='pairwise.complete.obs')

rm(expr.mtx)	### ! We don't need 'expr.mtx' any more. Save memory.

for ( k in 1:N_mod ) {
	cat("\t # Start of module", k, "...\n\n")

	orfs.each.mod = mod.list[[k]]$tgs
	N_orfs = length(orfs.each.mod)

	orfs.expr <- orfs.each.mod[orfs.each.mod %in% all.orfs]
	N_orfs.expr <- length(orfs.expr)

	orfs.mss <- setdiff(orfs.each.mod, orfs.expr)
	#! cf) orfs.mss <- orfs.each.mod[!orfs.each.mod %in% orfs.expr]
	if ( length(orfs.mss) != 0 ) {
		cat(file=FILENAME_3, paste(names(mod.list)[k], "\t", orfs.mss, "\n", sep=''), append=T)
	}

	if ( is.null(N_orfs.expr) ) {
		cat(file=FILENAME_4, paste(names(mod.list)[k], '\t0 ORF\n', sep=''), append=T)
		next
	} else if ( N_orfs.expr < 2 ) {
		cat(file=FILENAME_4, paste(names(mod.list)[k], '\t1 ORF\n', sep=''), append=T)
		next
	} else {

		# ================================
		#   1. Calculate the ZETA for the common target genes of TFs (assuming two TFs, X and Y).
		#      : Zeta_XY for the gene set, Genes_XY
		# ================================

		mod.pairs <- as.matrix(combn(orfs.expr, 2))
		mod.pairs.ind <- cbind(match(mod.pairs[1, ], all.orfs), 
			match(mod.pairs[2, ], colnames(pcc.mtx)))
		mod.pairs.pcc <- apply(mod.pairs.ind, 1, function(x) pcc.mtx[x[1], x[2]])
		ZETA.Mod.Real <- mean(abs(mod.pairs.pcc))


		# ================================
		#   2. Calculate the ZETA for the remaining target genes of X.
		#      : Zeta_X|Y for the gene set, Genes_X - Genes_XY (== Genes_X|Y)
		#   3. Calculate the ZETA for the remaining target genes of Y.
		#      : Zeta_Y|X for the gene set, Genes_Y - Genes_XY (== Genes_Y|X)
		# ================================

		tfs.each.mod = mod.list[[k]]$tfs
		ZETAs.for.TFs = vector('numeric', length = length(tfs.each.mod))
		names(ZETAs.for.TFs) = tfs.each.mod
		for ( ind_tf in 1:length(tfs.each.mod) ) {
			Genes_X = rownames(chip.data)[which(chip.data[, tfs.each.mod[ind_tf]] < chip.pval.thrs)]
			Genes_X = Genes_X[Genes_X %in% all.orfs]

			Genes_Xonly = Genes_X[!Genes_X %in% orfs.expr]

			if ( length(Genes_Xonly) >= 2 ) {
				mod.pairs.eachTF = as.matrix(combn(Genes_Xonly, 2))
				mod.pairs.ind.eachTF = cbind(match(mod.pairs.eachTF[1, ], all.orfs), 
					match(mod.pairs.eachTF[2, ], colnames(pcc.mtx)))
				mod.pairs.pcc.eachTF = apply(mod.pairs.ind.eachTF, 1, 
					function(x) pcc.mtx[x[1], x[2]])

				ZETA.Mod.eachTF = mean(abs(mod.pairs.pcc.eachTF))
				ZETAs.for.TFs[ind_tf] = ZETA.Mod.eachTF
			}
		}

		if ( any(ZETAs.for.TFs == 0) ) {
			Delta_Zeta = 0
			P_value = 1
			cat(file=FILENAME_1, names(mod.list)[k], "\t", paste(tfs.each.mod, collapse='_'), '\t', 
				N_orfs, "\t", N_orfs.expr, "\t", ZETA.Mod.Real, "\t", Delta_Zeta, "\t", 
				P_value, "\n", sep="", append=T)

			next
		}


		# ================================
		#   4. Decide which TF gives the maximum Zeta (most coherent). Assume that it is X.
		#      : TF_max = X such that Zeta_X|Y >= Zeta_Y|X
		# ================================

		Zeta_TFs_max = max(ZETAs.for.TFs)
		TF_max = names(which(ZETAs.for.TFs == Zeta_TFs_max))


		# ================================
		#   5. Decide if X and Y are synergistic by requiring the following condition.
		#      : Delta_Zeta = Zeta_XY - Zeta_X|Y > 0
		# ================================

		Delta_Zeta = ZETA.Mod.Real - Zeta_TFs_max

		if ( Delta_Zeta <= 0 ) {
			P_value = 1
			cat(file=FILENAME_1, names(mod.list)[k], "\t", 
				N_orfs, "\t", N_orfs.expr, "\t", ZETA.Mod.Real, "\t", Delta_Zeta, "\t", 
				P_value, "\n", sep="", append=T)

			next
		}


		# ================================
		#   6. Test if Delta_Zeta is statistically significant as follows.
		# ================================

		DELTAs_rnd = NULL
		for ( ind_rnd in 1:N_RANDOM ) {

			# >>> Target set of the TF_max <<< #
			Genes_X = rownames(chip.data)[which(chip.data[, TF_max] < chip.pval.thrs)]
			Genes_X = Genes_X[Genes_X %in% all.orfs]


			# ================================
			#   7. Divide the gene set, Genes_X, into two random groups corresponding to 
			#      the two sets, Genes_XY and Genes_X|Y, respectively.
			#      : Genes_XY_rnd and Genes_X|Y_rnd
			# ================================

			Genes_XY_rnd = sample(Genes_X, N_orfs.expr)
			Genes_Xonly_rnd = Genes_X[!Genes_X %in% Genes_XY_rnd]

			Genes_rnd = list(G_XY_rnd = Genes_XY_rnd, G_Xonly_rnd = Genes_Xonly_rnd)


			# ================================
			#   8. Calculate the difference between the Zeta's for the two random groups.
			#      : Delta_Zeta_rnd = Zeta_XY_rnd - Zeta_X|Y_rnd
			# ================================

			if ( length(Genes_Xonly_rnd) < 2 ) {
				DELTAs_rnd = rep(1, N_RANDOM)
				break
			} else {
				ZETAs_rnd = sapply(Genes_rnd, function(tmp.genes) {
					tmp.pairs <- as.matrix(combn(tmp.genes, 2))
					tmp.pairs.ind <- cbind(match(tmp.pairs[1, ], all.orfs), 
						match(tmp.pairs[2, ], colnames(pcc.mtx)))
					tmp.pairs.pcc <- apply(tmp.pairs.ind, 1, 
						function(x) pcc.mtx[x[1], x[2]])
					ZETA.tmp <- mean(abs(tmp.pairs.pcc))
				})

				Delta_Zeta_rnd = as.numeric(ZETAs_rnd['G_XY_rnd'] - ZETAs_rnd['G_Xonly_rnd'])
			}


			# ================================
			#   9. Make an emperical distribution of Delta_Zeta_rnd values by dividing Genes_X
			#      many times repeatedly (e.g. 1000 times).
			# ================================

			DELTAs_rnd = c(DELTAs_rnd, Delta_Zeta_rnd)
			cat(ind_rnd, ' ')
		}
		cat('\n')


		# ================================
		#  10. Calculate a P-value by counting the number of Delta_Zeta_rnd values which are 
		#      greater than or equal to the observed difference, Delta_Zeta.
		#      : P-value = sum(Delta_Zeta_rnd >= Delta_Zeta) / simul.times
		# ================================

		P_value = sum(DELTAs_rnd >= Delta_Zeta) / N_RANDOM

		cat(file=FILENAME_1, names(mod.list)[k], "\t", paste(tfs.each.mod, collapse='_'), '\t', 
			N_orfs, "\t", N_orfs.expr, "\t", ZETA.Mod.Real, "\t", Delta_Zeta, "\t", 
			P_value, "\n", sep="", append=T)


		# ================================
		#  11. To draw a histogram if necessary
		# ================================

#		if ( P_value == 1 ) next

#		tmp.range = c(DELTAs_rnd, Delta_Zeta)
#		jpeg(file=paste('./histograms/', names(mod.list)[k], '.jpg', sep=''), 
#		   quality=100, width=800, height=600)
#			hist(DELTAs_rnd, col='grey90', xlim=c(min(tmp.range), max(tmp.range)), 
#				main=paste(paste(tfs.each.mod, collapse='_'), 'vs.', TF_max))
#			abline(v=Delta_Zeta, lwd=5, col='red')
#		dev.off()

	}

	cat("\n\n\t\t ... Module", k, "Done!\n\n\n")

}

if ( ANSWER_ALL == 'yes' ) {
	# >>> To re-format the result table <<<
	foo.result <- as.matrix(read.table(file=FILENAME_1, sep="\t"))
	rownames(foo.result) <- NULL
	colnames(foo.result) <- c("Module", "TFs", "N_ORFs", "N_Expr", "ZETA", "DELTA", "P_value")

	foo.result <- foo.result[order(as.numeric(foo.result[, "P_value"])), ]
	rownames(foo.result) <- 1:nrow(foo.result)

	# >>> P-value correction by the FDR control <<<
	p.values <- as.numeric(foo.result[, "P_value"])
	foo.result <- cbind(foo.result, p.adjust(p.values, "fdr"))

	# >>> Q-value calculation if requested <<<
	if ( answer.qval == "yes" ) {
		library(qvalue)
		qobj <- qvalue(p.values)
		foo.result <- cbind(foo.result, qobj$qvalues)
	#	foo.result[, 5:9] <- formatC(as.numeric(foo.result[, 5:9]), format="e", digits=3)
		colnames(foo.result) <- c("Module", "TFs", "N_ORFs", "N_Expr", 
			"ZETA", "DELTA", "P_value", "P_FDR", "Q_value")
	} else if ( answer.qval == "no" ) {
	#	foo.result[, 5:8] <- formatC(as.numeric(foo.result[, 5:8]), format="e", digits=3)
		colnames(foo.result) <- c("Module", "TFs", "N_ORFs", "N_Expr", 
			"ZETA", "DELTA", "P_value", "P_FDR")
	}

	write.table(file=FILENAME_1, foo.result, sep="\t", quote=F)
}

time.end <- date()

cat("====================\n\n")
cat("The final result file, '", FILENAME_1, "', has been saved in the directory, '", 
	foo.path, "'.\n", sep="")
if ( file.exists(FILENAME_2) | file.exists(FILENAME_3) | file.exists(FILENAME_4) ) {
	cat("The following files also have been created in the same directory:\n")
	if ( file.exists(FILENAME_2) ) {
		cat(FILENAME_2, '\n')
	}
	if ( file.exists(FILENAME_3) ) {
		cat(FILENAME_3, '\n')
		cat("\n> There are ORFs with no expression values. They were saved in the following file :\n", 
			FILENAME_3, "\n", sep='', file = LOG_FILE, append=T)
	}
	if ( file.exists(FILENAME_4) ) {
		cat(FILENAME_4, '\n')
		cat("\n> There are modules with less than 2 ORFs. They were saved in the following file :\n", 
			FILENAME_4, "\n", sep='', file = LOG_FILE, append=T)
	}
}

cat("\nThe log file '", LOG_FILE, "' has been created in the same directory for the job summary.\n\n", sep='')

cat("\n!!!___Program Done___!!!\n\n")

cat(file=LOG_FILE, "\n> Program Start Time :\n", append=T)
cat(file=LOG_FILE, "\t", time.start, "\n", append=T)
cat(file=LOG_FILE, "\n> Program End Time :\n", append=T)
cat(file=LOG_FILE, "\t", time.end, "\n\n", append=T)

cat("Program Start Time:\n")
cat("\t", time.start, "\n")
cat("Program End Time:\n")
cat("\t", time.end, "\n\n")

}

