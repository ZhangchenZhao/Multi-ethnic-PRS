library(data.table)
library(lassosum)
library(Matrix)
library(parallel)


ped=read.table("/net/snowwhite/home/zczhao/PRS/Pipeline/Pheno/AFR_8traits_1113.ped",header=T);Y=ped[,c("FID","IID","LDL")]; colnames(Y)[3]="pheno"; Covar=ped[,c("FID","IID","Sex","BY","P1","P2","P3","P4")]
train_file="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink";test_file="/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/African/African_2ktrain_plink"
sum_stats=fread(paste0("/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/Method1/Training/out_","lassosum","_","AFR_5k","_","LDL","_","eur",".txt")); colnames(sum_stats)=c("SNP","A1","Beta"); LDblocks="EUR.hg19";Ytype="linear"
beta.file=paste0("/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/Method1/Training/out_","lassosum","_","AFR_5k","_","LDL","_","eur",".txt")


save.image(file="/net/snowwhite/home/zczhao/PRS/Pipeline/Rpackage/TestData.RData")

system.time({out.beta=PRS_TransferLearning(Y,Covar,train_file,test_file,sum_stats,LDblocks,Ytype="linear")})




ped=read.table("/net/snowwhite/home/zczhao/PRS/Pipeline/Pheno/AFR_8traits_1113.ped",header=T);Y=ped[,c("FID","IID","HDL")]; colnames(Y)[3]="pheno"; Covar=ped[,c("FID","IID","Sex","BY","P1","P2","P3","P4")]
train_file="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink";test_file="/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/African/African_2ktrain_plink"
sum_stats=fread(paste0("/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/Method1/Training/out_","lassosum","_","AFR_5k","_","HDL","_","eur",".txt")); colnames(sum_stats)=c("SNP","A1","Beta"); LDblocks="EUR.hg19";Ytype="linear"
beta.file=paste0("/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/Method1/Training/out_","lassosum","_","AFR_5k","_","HDL","_","eur",".txt")


save.image(file="/net/snowwhite/home/zczhao/PRS/Pipeline/Rpackage/TestData.RData")

system.time({out.beta=PRS_TransferLearning(Y,Covar,train_file,test_file,sum_stats,LDblocks,Ytype="linear",beta.file)})

library(snpStats)
bed="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink.bed"
bim="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink.bim"
fam="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink.fam"
snplist=read.table(bim,nrow=5)$V2

system.time({a=read.plink(bed, bim, fam,select.snps = snplist)})
