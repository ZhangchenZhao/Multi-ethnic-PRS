###PRSCS###
#####the PRScs package and the reference file can be downloaded at https://github.com/getian107/PRScs.
###Summary file for prscs: out.temp.file
###dat0 is the summary file.
###test.bfile is the plink file that you used for model fitting.
write.table(data.frame("SNP"=dat0$SNP,"A1"=dat0$A1,"A2"=dat0$A2,"BETA"=dat0$beta1,"P"=dat0$p1),file=out.temp.file,col.names=TRUE,row.names=FALSE,quote=FALSE)
cmd=paste0("python /net/snowwhite/home/zczhao/PRS/python/PRScs/PRScs.py --ref_dir=","/net/snowwhite/home/zczhao/PRS/Ref/ldblk_1kg_eas",
" --bim_prefix=",test.bfile," --sst_file=", out.temp.file," --n_gwas=",median(dat0$N1)," --out_dir=/net/csgspare2/spare1/snowwhite.archive/zczhao/PRS_Geno/scratch/out_prscs_",Cat,"_",kword,"_",ijk,"_eas")
system(cmd)

