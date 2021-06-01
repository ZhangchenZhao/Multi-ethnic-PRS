##dat0 is the summary file
library(lassosum)
ref.bfile <- "/net/snowwhite/home/zczhao/PRS/1000G/EASID_1000G_plink" ###reference file, usually 1000 Genome
###test.bfile is the plink file you used for model training.
##LDblocks can be ASN, EUR, etc. Details can be found at https://github.com/tshmak/lassosum.
cor <- p2cor(p =dat0$p1, n =dat0$N1, sign=dat0$beta1) 
out1 <- lassosum.pipeline(cor=cor, chr=dat0$Chr, pos=dat0$Pos, snp=dat0$SNP,
			A1=dat0$A1,A2=dat0$A2,
                        ref.bfile=ref.bfile, test.bfile=test.bfile,  exclude.ambiguous=FALSE,
                        LDblocks = "ASN.hg19",destandardize=F) 

v1 <- validate(out1, pheno=Y,covar=Covar) ### all result can be found in v1.
print(v1$best.s); print(v1$best.lambda)  ## Best tuning parameters returned by lassosum.    
ss1=cbind(out1$sumstats,v1$best.beta) ##This ss1 contains the best beta to construct PRS by lassosum. 

