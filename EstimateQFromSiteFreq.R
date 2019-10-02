#
#	Estimate scaled rate matrix Q from site frequency data assuming
#
#
#	1) General 4x4 Q
#	2) Reversible 4x4 Q
#	3) Strand-symmetric 4x4 Q
#
#
#	02/10/19 Conrad Burden
#
#	setwd("")	# set working directory here
#
#
#	Read in site frequency Data
#
	siteFreqData <- read.table(file="BergmanData.txt")
#
	M <- dim(siteFreqData)[1] - 1		# number of individuals in population sampled
	HM <- sum(1/(1:(M - 1)))
	yArray <- 1:(M - 1)
	oneOverYSum <- 1/(M - yArray) + 1/yArray
	oneOverYDiff <- 1/(M - yArray) - 1/yArray  
#
	L_A <- siteFreqData$A.C[1]
	L_C <- siteFreqData$C.G[1]
	L_G <- siteFreqData$G.T[1]
	L_T <- siteFreqData$G.T[M + 1]
	L_i <- c(L_A, L_C, L_G, L_T)
#
	L_ij_ofY <- siteFreqData[2:M,]
	L_ij <- colSums(L_ij_ofY)
	L_ACofY <- L_ij_ofY$A.C
	L_AGofY <- L_ij_ofY$A.G
	L_ATofY <- L_ij_ofY$A.T
	L_CGofY <- L_ij_ofY$C.G
	L_CTofY <- L_ij_ofY$C.T
	L_GTofY <- L_ij_ofY$G.T
	L_AC <- L_CA <- sum(L_ACofY)
	L_AG <- L_GA <- sum(L_AGofY)
	L_AT <- L_TA <- sum(L_ATofY)
	L_CG <- L_GC <- sum(L_CGofY)
	L_CT <- L_TC <- sum(L_CTofY)
	L_GT <- L_TG <- sum(L_GTofY)
#
	L_m <- sum(L_i)
	L_p <- sum(L_ij)
	L_total <- L_m + L_p
#
###########################################################################
#
#	Estimate of Q for a reversible 4x4 Q 
#
	piHat_A_rev <- (L_A + (L_AC + L_AG + L_AT)/2)/L_total 
	piHat_C_rev <- (L_C + (L_AC + L_CG + L_CT)/2)/L_total 
	piHat_G_rev <- (L_G + (L_AG + L_CG + L_GT)/2)/L_total 
	piHat_T_rev <- (L_T + (L_AT + L_CT + L_GT)/2)/L_total 
	CHat_ij_rev <- L_ij/(2*L_total*HM)
#
#	Reconstruct reversible Q estimate (nucleotide ordering (A, T, G, C))
#
	piHat_rev <- c(piHat_A_rev, piHat_T_rev, piHat_G_rev, piHat_C_rev)
	CHat_rev <- t(matrix(c(	  0 , L_AT, L_AG, L_AC, 
							L_TA,   0 , L_TG, L_TC, 
							L_GA, L_GT,   0 , L_GC, 
							L_CA, L_CT, L_CG,   0  )/(2*L_total*HM), 
					  nrow=4, ncol=4))
	diag(CHat_rev) <- -rowSums(CHat_rev, na.rm=TRUE)	
#
	QHat_rev <- diag(1/piHat_rev) %*% CHat_rev
	dimnames(QHat_rev) <- 
			list(c("A", "T", "G", "C"), c("A", "T", "G", "C"))
	cat("\n QHat_rev =\n")
	options(digits=4)
	print(QHat_rev)
	options(digits=7)
#					
###########################################################################
#
#	Estimate of Q for a general 4x4 Q 
#
	CHat <- L_p/(2*L_total*HM)
	piPrimeHat_A <- L_A/L_m
	piPrimeHat_C <- L_C/L_m
	piPrimeHat_G <- L_G/L_m
	piPrimeHat_T <- L_T/L_m
#
#	Initial rough estimate of phi_ij using Eq.(48) of Burden&Tang(2017)
#
#	Function to evaluate log of the likelihood h(l_p ; cHat_rev_ij, phi_ij)
#	params[1:3] =  
#		(phi_AC, phi_AG, phi_CG)
#	and the remaining phi's are
#		phi_AT = - phi_AC - phi_AG = -params[1] - params[2]
#		phi_CT = + phi_AC - phi_CG = params[1] - params[3]
#		phi_GT = + phi_AG + phi_CG = params[2] + params[3]
#
	cHat_ij_rev <- CHat_ij_rev/CHat
	minus_log_h_init <- function(params){
		phi_ij <- c(params[1], params[2], -params[1] - params[2], 
						params[3], params[1] - params[3], params[2] + params[3])
		log_h_init <- sum(L_ij_ofY *
					log(outer(oneOverYSum, cHat_ij_rev) + outer(oneOverYDiff, phi_ij)))
		minus_log_h_init <- -log_h_init
		return(minus_log_h_init)
		}
#
#	Linear constriants for constrOptim(): init. rough estimate
#
	ui <- t(matrix(c(	 -1, 0, 0, 	#   - phi_AC > - c_AC 
					 	 +1, 0, 0, 	#   + phi_AC > - c_AC 
						  0,-1, 0, 	#   - phi_AG > - c_AG 
						  0,+1, 0, 	#   + phi_AG > - c_AG  
						 -1,-1, 0, 	#   - phi_AT > - c_AT 
						 +1,+1, 0, 	#   + phi_AT > - c_AT 
						  0, 0,-1, 	#   - phi_CG > - c_CG 
						  0, 0,+1, 	#   + phi_CG > - c_CG 
						 -1, 0,+1, 	#   - phi_CT > - c_CT  
						 +1, 0,-1, 	#   + phi_CT > - c_CT  
						 0,-1,-1, 	#   - phi_GT > - c_GT 
						 0,+1,+1) 	#   + phi_GT > - c_GT 
						, nrow=3, ncol=12))
	ci <- -rep(cHat_ij_rev, each=2)
#
#	Max of - log_h_init
#
	opt_init <- constrOptim(theta=c(0, 0, 0), 
					f= minus_log_h_init, 
					ui=ui, ci=ci, method="Nelder-Mead",
					control=list(reltol=1.e-14, maxit=100000))
#
#	Now estimate all c_ij & phi_ij using c_ij_rev and initial estimate 
#			of phi_ij as a starting point
#
#	Function to evaluate log of the likelihood h(l_p ; c_ij, phi_ij)
#		params are ordered
#	params[1:8] =  
#		(c_AC, c_AG, c_AT, c_CG, c_CT, phi_AC, phi_AG, phi_CG)
#	where the remaining c is 
#		c_GT = 1 - sum(the first 8 c's)
#	and the remaining phi's are
#		phi_AT = - phi_AC - phi_AG = -params[6] - params[7]
#		phi_CT = + phi_AC - phi_CG = params[6] - params[8]
#		phi_GT = + phi_AG + phi_CG = params[7] + params[8]
#
	minus_log_h <- function(params){
		c_ij <- c(params[1:5], 1 - sum(params[1:5]))
		phi_ij <- c(params[6], params[7], -params[6] - params[7], 
						params[8], params[6] - params[8], params[7] + params[8])
		log_h <- sum(L_ij_ofY *
					log(outer(oneOverYSum, c_ij) + outer(oneOverYDiff, phi_ij)))
		minus_log_h <- -log_h
		return(minus_log_h)
		}
#
#	Linear constriants for constrOptim()
#
	ui <- t(matrix(c(	 1, 0, 0, 0, 0,	-1, 0, 0, 	#  c_AC - phi_AC > 0 
					 	 1, 0, 0, 0, 0,	+1, 0, 0, 	#  c_AC + phi_AC > 0 
						 0, 1, 0, 0, 0,	 0,-1, 0, 	#  c_AG - phi_AG > 0 
						 0, 1, 0, 0, 0,	 0,+1, 0, 	#  c_AG + phi_AG > 0 
						 0, 0, 1, 0, 0,	-1,-1, 0, 	#  c_AT - phi_AT > 0 
						 0, 0, 1, 0, 0,	+1,+1, 0, 	#  c_AT + phi_AT > 0 
						 0, 0, 0, 1, 0,	 0, 0,-1, 	#  c_CG - phi_CG > 0 
						 0, 0, 0, 1, 0,	 0, 0,+1, 	#  c_CG + phi_CG > 0 
						 0, 0, 0, 0, 1,	-1, 0,+1, 	#  c_CT - phi_CT > 0 
						 0, 0, 0, 0, 1,	+1, 0,-1, 	#  c_CT + phi_CT > 0 
						-1,-1,-1,-1,-1,	 0,-1,-1, 	#  c_GT - phi_GT > 0 
						-1,-1,-1,-1,-1,	 0,+1,+1) 	#  c_GT + phi_GT > 0 
						, nrow=8, ncol=12))
	ci <- c(rep(0, 10), -1, -1)
#
#	Max of log_h
#
	theta0 <- c(cHat_ij_rev[1:5], opt_init$par) # starting parameters for opt
	opt <- constrOptim(theta=theta0, 
					f= minus_log_h, 
					ui=ui, ci=ci, method="Nelder-Mead",
					control=list(reltol=1.e-15, 
					parscale= theta0))
#
#	Reconstruct general Q estimate (nucleotide ordering (A, T, G, C))
#
	c_Hat_ij <- c(opt$par[1:5], 1 - sum(opt$par[1:5]))
#
	c_Hat_AC <- c_Hat_ij[1]; c_Hat_AG <- c_Hat_ij[2]; c_Hat_AT <- c_Hat_ij[3] 
	c_Hat_CG <- c_Hat_ij[4]; c_Hat_CT <- c_Hat_ij[5]; c_Hat_GT <- c_Hat_ij[6]
	C_Hat_matrix <- t(matrix(c(	  	0   , c_Hat_AT, c_Hat_AG, c_Hat_AC, 
							 c_Hat_AT,     0   , c_Hat_GT, c_Hat_CT, 
							 c_Hat_AG, c_Hat_GT,    0    , c_Hat_CG, 
							 c_Hat_AC, c_Hat_CT, c_Hat_CG,    0    )*CHat, 
					  			nrow=4, ncol=4))
	diag(C_Hat_matrix) <- -rowSums(C_Hat_matrix, na.rm=TRUE)	
#
	phi_Hat_AC <- opt$par[6]; phi_Hat_CA <- -phi_Hat_AC 
	phi_Hat_AG <- opt$par[7]; phi_Hat_GA <- -phi_Hat_AG
	phi_Hat_CG <- opt$par[8]; phi_Hat_GC <- -phi_Hat_CG		 
	phi_Hat_AT <- - phi_Hat_AC - phi_Hat_AG; phi_Hat_TA <- -phi_Hat_AT
	phi_Hat_CT <- + phi_Hat_AC - phi_Hat_CG; phi_Hat_TC <- -phi_Hat_CT 
	phi_Hat_GT <- + phi_Hat_AG + phi_Hat_CG; phi_Hat_TG <- -phi_Hat_GT 
#
	Phi_Hat_matrix <- t(matrix(c(	  	0    , phi_Hat_AT, phi_Hat_AG, phi_Hat_AC, 
							phi_Hat_TA,      0    , phi_Hat_TG, phi_Hat_TC, 
							phi_Hat_GA, phi_Hat_GT,     0     , phi_Hat_GC, 
							phi_Hat_CA, phi_Hat_CT, phi_Hat_CG,     0      )*CHat, 
					  			nrow=4, ncol=4))
#
	pi_Hat_A <- piPrimeHat_A*(1 - 2*HM*CHat) - HM*C_Hat_matrix[1, 1]
	pi_Hat_T <- piPrimeHat_T*(1 - 2*HM*CHat) - HM*C_Hat_matrix[2, 2]
	pi_Hat_G <- piPrimeHat_G*(1 - 2*HM*CHat) - HM*C_Hat_matrix[3, 3]
	pi_Hat_C <- piPrimeHat_C*(1 - 2*HM*CHat) - HM*C_Hat_matrix[4, 4]
	pi_Hat <- c(pi_Hat_A, pi_Hat_T, pi_Hat_G, pi_Hat_C)
#
	QHat <- diag(1/pi_Hat) %*% (C_Hat_matrix + Phi_Hat_matrix)
	dimnames(QHat) <- 
			list(c("A", "T", "G", "C"), c("A", "T", "G", "C"))
	cat("\n QHat =\n")
	options(digits=4)
	print(QHat)
	options(digits=7)
#
###########################################################################
#
#	Estimate of Q for a strand-symmetric Q 
#
#
#	Estimates of e' and b': first define log-likelihood, Eq.(382) blue notebook
#
	LOfEPrimeBPrime <- function(ePrimeBPrime){
#
		ePrime <- ePrimeBPrime[1]
		bPrime <- ePrimeBPrime[2]
		logL_eb <- 
			sum((L_ACofY + L_GTofY[M - yArray])*log(ePrime/(M - yArray) + (1 - bPrime)/yArray) +
			    (L_AGofY + L_CTofY[M - yArray])*log((1 - ePrime)/(M - yArray) + bPrime/yArray))
		return(-logL_eb)
		}
#
	uiEPrimeBPrime <- as.matrix(t(array(c(	 1, 0, 
											-1, 0, 
											 1, 0, 
											-1, 0), dim=c(2,4))))
#
	ciEPrimeBPrime <- as.matrix(c(	0, -1, 0, -1))
#
	EPrimeBPrimeHat <- constrOptim(theta=c(0.5, 0.5), 
				f=LOfEPrimeBPrime, 
				ui=uiEPrimeBPrime, ci=ciEPrimeBPrime, method="Nelder-Mead", 
				control=list(reltol=1.e-16))$par
#
	ePrimeHat <- EPrimeBPrimeHat[1]								
	bPrimeHat <- EPrimeBPrimeHat[2]								
#
#	Estimates of a, b, ..., f
#
	L_bracketsAT.CG <- L_AC + L_AG + L_TC + L_TG
	L_bracketsAT <- L_A + L_T + L_AT
	L_bracketsCG <- L_C + L_G + L_CG
#
	denominatorAT <- L_bracketsAT + 0.5*L_bracketsAT.CG
	denominatorCG <- L_bracketsCG + 0.5*L_bracketsAT.CG
#
	aHat <- Q_ATHat <- Q_TAHat <- L_AT/denominatorAT/HM
	bHat <- Q_CTHat <- Q_GAHat <- bPrimeHat*L_bracketsAT.CG/denominatorCG/2/HM
	cHat <- Q_AGHat <- Q_TCHat <- (1 - ePrimeHat)*L_bracketsAT.CG/denominatorAT/2/HM
	dHat <- Q_GTHat <- Q_CAHat <- (1 - bPrimeHat)*L_bracketsAT.CG/denominatorCG/2/HM
	eHat <- Q_ACHat <- Q_TGHat <- ePrimeHat*L_bracketsAT.CG/denominatorAT/2/HM
	fHat <- Q_CGHat <- Q_GCHat <- L_CG/denominatorCG/HM
#
#	Construct QHat_SS
#
	QHat_SS <- t(matrix(c(	   0 , aHat, cHat, eHat, 
							 aHat,   0 , eHat, cHat, 
							 bHat, dHat,   0 , fHat, 
							 dHat, bHat, fHat,   0 ), 
					  			nrow=4, ncol=4))
	diag(QHat_SS) <- -rowSums(QHat_SS, na.rm=TRUE)	
	dimnames(QHat_SS) <- 
			list(c("A", "T", "G", "C"), c("A", "T", "G", "C"))
	cat("\n QHat_SS =\n")
	options(digits=4)
	print(QHat_SS)
	options(digits=7)
#					
###########################################################################
	