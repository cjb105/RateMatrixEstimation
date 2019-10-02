#
#	1) After running <EstimateQFromSiteFreq.R> for maximum likelihood estimator 
#	of the scaled rate matrix Q, plot log likelihood in vicinity of estimator 
#	to assess its accuracy
#
#	2) Calculate p-values to test significance in nested pairs of models 
#			Reversible Q < General Q  and  Strand-symmetric Q < General Q
#
#	3) Calculate heterozygosity for General Q, Reversible Q & Strand-symmetric Q
#
#	02/10/19 Conrad Burden
#
	if(.Platform$OS.type=="unix") {
		quartz()
		}
	if(.Platform$OS.type=="windows") {
		windows()
		}
#
	oldpar <- par(mfrow=c(3,4))
	xLabels <- c(expression(Q_AT), expression(Q_AG), expression(Q_AC), 
				expression(Q_TG), expression(Q_TC), expression(Q_GC), 
				expression(Q_TA), expression(Q_GA), expression(Q_CA), 
				expression(Q_GT), expression(Q_CT), expression(Q_CG)) 
#
#	arrays of counts, ordered (A, T, G, C)
#
	L <- c(L_A, L_T, L_G, L_C)
	LofY <- cbind(L_ATofY, L_AGofY, L_ACofY, 
				L_GTofY[(M - 1):1], L_CTofY[(M - 1):1], L_CGofY[(M - 1):1])
#
#	Function to evaluate log likelihood as a function of QQ
#
	logL <- function(QQ){
		pi_i <- eigen(t(QQ))$vectors[,4]
		pi_i <- pi_i/sum(pi_i)
		pi_iQij <- c(pi_i[1]*QQ[1,2], pi_i[1]*QQ[1,3], pi_i[1]*QQ[1,4], 
					pi_i[2]*QQ[2,3], pi_i[2]*QQ[2,4], pi_i[3]*QQ[3,4])
		pi_jQji <- c(pi_i[2]*QQ[2,1], pi_i[3]*QQ[3,1], pi_i[4]*QQ[4,1], 
					pi_i[3]*QQ[3,2], pi_i[4]*QQ[4,2], pi_i[4]*QQ[4,3])
		term1 <- sum(L*log(pi_i*(1 + HM*diag(QQ))))
		term2 <- sum(LofY*
			log(outer(1/((M - 1):1), pi_iQij) + outer(1/(1:(M - 1)), pi_jQji)))
		return(term1 + term2)
		}
#
	indexList <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3, 4), 
					c(2,1), c(3,1), c(4,1), c(3,2), c(4,2), c(4,3))
	for(iArg in 1:12){
		index1 <- indexList[[iArg]][1] 
		index2 <- indexList[[iArg]][2] 
		y <- function(x){
			QQ <- QHat
			QQ[index1, index2] <- x*QQ[index1, index2]
			diag(QQ) <- 0
			diag(QQ) <- -rowSums(QQ)
			return(logL(QQ))
			}
#		curve(Vectorize(y)(x), from=0.9995, to = 1.0005, 
		curve(Vectorize(y)(x), from=0.9, to = 1.1, 
				main=xLabels[iArg], 
				xlab=paste("x", QHat[index1, index2]), ylab="log-likelihood")
		abline(v=1)
		abline(h=y(1))
		}
	par(oldpar)
#
#	Likelihood ratio tests
#
	logLGeneral <- logL(QHat)
	logLRev <- logL(QHat_rev)
	logLSS <- logL(QHat_SS)
#
	pVal_gen_rev_log10 <- pchisq(2*(logLGeneral - logLRev), 
							df=3, 
							lower.tail=FALSE, 
							log.p=TRUE)/log(10)
	pVal_gen_SS_log10 <- pchisq(2*(logLGeneral - logLSS), 
							df=6, 
							lower.tail=FALSE, 
							log.p=TRUE)/log(10)
#
	cat("\n Pvalue Q general to Q rev: ", 10^pVal_gen_rev_log10, "      log Pvalue: ", pVal_gen_rev_log10)
	cat("\n Pvalue Q general to Q SS : ", 10^pVal_gen_SS_log10, "      log Pvalue: ", pVal_gen_SS_log10)
#
#	Calculate heterozygosity
#
	heterozygosity <- function(QQ){
		pi_i <- eigen(t(QQ))$vectors[,4]
		pi_i <- pi_i/sum(pi_i)
		print(pi_i)
		heterozygosity  <- -sum(pi_i*diag(QQ))
		return(heterozygosity)
		}
#
	cat("\n \n Heterozygosities: ", "\n    General Q:", heterozygosity(QHat),
								    "\n Reversible Q:", heterozygosity(QHat_rev),
								    "\n Strand-sym Q:", heterozygosity(QHat_SS))
#
