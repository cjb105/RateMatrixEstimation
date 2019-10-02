# RateMatrixEstimation
Programs to accompany the paper "Maximum likelihood estimators for scaled mutation rates in an equilibrium mutation-drift model” by Claus Vogl, Lynette Mikula and Conrad Burden

R Programs:

EstimateQFromSiteFreq.R 

	Calculate maximum likelihood estimates of the scaled rate matrix Q from site frequency data assuming (1) General 4x4 Q; (2) Reversible 4x4 Q; and (3) Strand-symmetric 4x4 Q

CheckMinimumGeneralQ.R

	After running <EstimateQFromSiteFreq.R> for maximum likelihood estimator of the scaled rate matrix Q, (1) plot log likelihood in vicinity of estimator to assess its accuracy; 2) Calculate p-values to test significance in nested pairs of models (Reversible Q < General Q)  and  (Strand-symmetric Q < General Q); and 3) Calculate heterozygosity for General Q, Reversible Q & Strand-symmetric Q

Test Dataset: 

BergmanData.txt

	Site frequency counts from Bergman et al (2018), Genome Biology and Evolution 10, 269-275, Supplementary table S4.
