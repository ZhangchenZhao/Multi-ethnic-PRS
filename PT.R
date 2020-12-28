###a function to run linear regression fitting with PRS#######
lm_result_Bayes<-function(PRS,Y,Covar ,name1){
	datatemp=merge(PRS,Y,by.x="FID",by.y="FID")
	datatemp2=merge(datatemp,Covar,by.x="FID",by.y="FID")
	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	modeltemp0=lm(pheno~Sex+BY+P1+P2+P3+P4, data=datatemp2)
	modeltemp=lm(pheno~PRS+Sex+BY+P1+P2+P3+P4, data=datatemp2)

	prs.coef <- summary(modeltemp)$coeff[2,]
	prs.beta <- as.numeric(prs.coef[1])
	prs.se <- as.numeric(prs.coef[2])
	prs.p <- as.numeric(prs.coef[4])
	model0.r2=summary(modeltemp0)$r.squared
	model1.r2=summary(modeltemp)$r.squared
	prs.r2 = summary(modeltemp)$r.squared-model0.r2
	list.out=cbind(class=name1, prs.beta,prs.se,prs.p, model0.r2, model1.r2, prs.r2)
	return(list.out)
}


###We run P+T through plink software.
##out.temp.file is the data.frame("SNP"=dat0$SNP,"A1"=dat0$A1, "Beta"=dat0$beta1)
##out.temp.file2 is the data.frame("SNP"=dat0$SNP,"P"=dat0$p1)
##range_list1 is the file you set a range, in my file, it is "S1  0  1"
ref.bfile <- "/net/snowwhite/home/zczhao/PRS/1000G/EASID_1000G_plink_no3alleles_forPT"

	out.all1=data.frame()
	for (RLD2 in c(0.1, 0.2, 0.5, 0.8)){
		for (PT in c(1.0, 0.8, 0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.05, 0.02, 0.01, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8)){
			
			cmd=paste0("plink-1.9 --bfile ",ref.bfile,"  --clump ",out.temp.file2, " --clump-p1 ",PT,"  --clump-r2 ",RLD2,
		" --out ", clump.file)
			system(cmd)
			cmd=paste0("plink-1.9 --bfile ",test.bfile,"  --score ",out.temp.file," sum --q-score-range /net/snowwhite/home/zczhao/PRS/range_list1 ",
		clump.file,".clumped 3 5 "," --out ",prs.file)
			system(cmd)
			resulttemp=read.table(paste0(prs.file, ".S1.profile"),header=T)
			parm=paste0(RLD2,"_",PT)
			if (parm=="0.1_1") {
					out.all1=c(lm_result_Bayes(resulttemp,Y,Covar,parm)[1:7],0);
			} else {
					out.all1=rbind(out.all1,c(lm_result_Bayes(resulttemp,Y,Covar,parm)[1:7],0))
			}

		}
	} 
####now out.all1 is the summary output from the grid of all parameters. You can find the best parameter based on r2.
