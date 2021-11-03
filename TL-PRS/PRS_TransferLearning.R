library(data.table)
library(lassosum)
library(Matrix)
library(parallel)

nagelkerke<-function (fit, null = NULL, restrictNobs = FALSE){
    TOGGLE = (class(fit)[1] == "lm" | class(fit)[1] == "gls" |
        class(fit)[1] == "lme" | class(fit)[1] == "glm" | class(fit)[1] ==
        "negbin" | class(fit)[1] == "zeroinfl" | class(fit)[1] ==
        "clm" | class(fit)[1] == "vglm" | class(fit)[1] == "betareg" |
        class(fit)[1] == "rq")
    BOGGLE = (class(fit)[1] == "nls" | class(fit)[1] == "lmerMod" |
        class(fit)[1] == "glmerMod" | class(fit)[1] == "merModLmerTest" |
        class(fit)[1] == "lmerModLmerTest" | class(fit)[1] ==
        "clmm")
    SMOGGLE = (class(fit)[1] == "lmerMod" | class(fit)[1] ==
        "glmerMod" | class(fit)[1] == "merModLmerTest" | class(fit)[1] ==
        "lmerModLmerTest" | class(fit)[1] == "vglm")
    ZOGGLE = (class(fit)[1] == "zeroinfl")
    ZOGGLE2 = (class(fit)[1] == "rq")
    NOGGLE = is.null(null)
    ERROR = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    if (!restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1)
    }
    if (restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1, data = fit$model)
    }
    if (restrictNobs & !NOGGLE) {
        null = update(null, data = fit$model)
    }
    if (NOGGLE & BOGGLE) {
        ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"
    }
    if ((!TOGGLE) & (!BOGGLE)) {
        ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"
    }
    SMOGGLE2 = (class(null)[1] == "lmerMod" | class(null)[1] ==
        "glmerMod" | class(null)[1] == "merModLmerTest" | class(null)[1] ==
        "lmerModLmerTest" | class(null)[1] == "vglm")
    Y = matrix(rep(NA, 2), ncol = 1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    Z = matrix(rep(NA, 3), ncol = 1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", "Nagelkerke (Cragg and Uhler)")
    X = matrix(rep(NA, 4), ncol = 4)
    colnames(X) = c("Df.diff", "LogLik.diff", "Chisq", "p.value")
    rownames(X) = ""
    U = matrix(rep(NA, 2), ncol = 1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    if (TOGGLE | BOGGLE) {
        if (!SMOGGLE) {
            Y[1] = toString(fit$call)
        }
        if (SMOGGLE) {
            Y[1] = toString(fit@call)
        }
    }
    if (TOGGLE | (BOGGLE & !NOGGLE)) {
        if (!SMOGGLE2) {
            Y[2] = toString(null$call)
        }
        if (SMOGGLE2) {
            Y[2] = toString(null@call)
        }
        if (!ZOGGLE & !ZOGGLE2) {
            N = nobs(fit)
            U[1, 1] = nobs(fit)
            U[2, 1] = nobs(null)
        }
        if (!ZOGGLE & ZOGGLE2) {
            N = length(fit$y)
            U[1, 1] = length(fit$y)
            U[2, 1] = length(null$y)
        }
        if (ZOGGLE) {
            N = fit$n
            U[1, 1] = fit$n
            U[2, 1] = null$n
        }
        if (U[1, 1] != U[2, 1]) {
            ERROR2 = "WARNING: Fitted and null models have different numbers of observations"
        }
        m = suppressWarnings(logLik(fit, REML = FALSE))[1]
        n = suppressWarnings(logLik(null, REML = FALSE))[1]
        mf = 1 - m/n
        Z[1, ] = signif(mf, digits = 6)
        cs = 1 - exp(-2/N * (m - n))
        Z[2, ] = signif(cs, digits = 6)
        nk = cs/(1 - exp(2/N * n))
        Z[3, ] = signif(nk, digits = 6)
        o = n - m
        dfm = attr(logLik(fit), "df")
        dfn = attr(logLik(null), "df")
        if (class(fit)[1] == "vglm") {
            dfm = df.residual(fit)
        }
        if (class(fit)[1] == "vglm") {
            dfn = df.residual(null)
        }
        dff = dfn - dfm
        CHI = 2 * (m - n)
        P = pchisq(CHI, abs(dff), lower.tail = FALSE)
        X[1, 1] = dff
        X[1, 2] = signif(o, digits = 5)
        X[1, 3] = signif(CHI, digits = 5)
        X[1, 4] = signif(P, digits = 5)
    }
    W = ERROR
    WW = ERROR2
    V = list(Y, Z, X, U, W, WW)
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null",
        "Likelihood.ratio.test", "Number.of.observations", "Messages",
        "Warnings")
    return(V)
}




##resulttemp,ped,Covar_name,Y_name
logistic_result_generator_old<-function(PRS,ped,Covar_name,Y_name){
        datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
  	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
        modeltemp0=glm(paste0(Y_name,"~",args), data=datatemp2, family = "binomial")
        modeltemp=glm(paste0(Y_name,"~ PRS +",args), data=datatemp2, family = "binomial")
	#casecontrol=paste0(sum(datatemp2$pheno==1,na.rm=T),":",sum(datatemp2$pheno==0,na.rm=T))
	r2=1-modeltemp$deviance/modeltemp0$deviance
	#ci=try(exp(cbind(coef(modeltemp), confint(modeltemp))["PRS",1:3] ),silent=TRUE)
	#if (class(ci)!="try-error"){
	#	or=ci[1]
	#	ci2=paste0("(",round(ci[2],2),",",round(ci[3],2),")")
	#} else {
	#	or=NA;ci2=NA
	#}
	#prs.p=coef(summary( modeltemp))["PRS",4]
	#
        #list.out=cbind(class=name1,casecontrol,r2,or,ci2, prs.p,r2)
        return(r2)
}

##resulttemp,ped,Covar_name,Y_name
logistic_result_generator<-function(PRS,ped,Covar_name,Y_name){
        datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
  	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
	if (args!=""){
        	modeltemp0=glm(paste0(Y_name,"~",args), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~ PRS +",args), data=datatemp2, family = "binomial")
	} else {
        	modeltemp0=glm(paste0(Y_name,"~1"), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~PRS"), data=datatemp2, family = "binomial")
	}
	r2=nagelkerke(modeltemp,null=modeltemp0)[[2]][3]
        return(r2)
}


##resulttemp,ped,Covar_name,Y_name
linear_result_generator<-function(PRS,ped,Covar_name,Y_name){
	datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
	if (args!=""){
		modeltemp0=lm(paste0(Y_name,"~",args), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~ PRS +",args), data=datatemp2)
	} else{
		modeltemp0=lm(paste0(Y_name,"~1"), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~PRS"), data=datatemp2)
	}
	#prs.coef <- summary(modeltemp)$coeff[2,]
	#prs.beta <- as.numeric(prs.coef[1])
	#prs.se <- as.numeric(prs.coef[2])
	#prs.p <- as.numeric(prs.coef[4])
	model0.r2=summary(modeltemp0)$r.squared
	model1.r2=summary(modeltemp)$r.squared
	prs.r2 = summary(modeltemp)$r.squared-model0.r2
	#list.out=cbind(class=name1, prs.beta,prs.se,prs.p, model0.r2, model1.r2, prs.r2)
	return(prs.r2)
}

result_generator<-function(PRS,ped,Covar_name,Y_name,Ytype){
	if (Ytype=="C"){
		return(linear_result_generator(PRS,ped,Covar_name,Y_name))
	} else {
		return(logistic_result_generator(PRS,ped,Covar_name,Y_name))
	}
}

##test.bfile=train_file;beta.file=sum_stats_file;
##covar=Covar_name;kword=Y_name;temp.file=paste0(tempfile,"_step0") 
calculate_betaPRS<-function(test.bfile,beta.file,ped,covar,kword,temp.file){
	cmd=paste0("plink-1.9 --bfile ",test.bfile,"  --score ",beta.file, " sum", " --out ",temp.file,".train.PRS")
	system(cmd)
	temp=read.table(paste0(temp.file,".train.PRS.profile"),header=T)[,c(1,2,6)]
	temp2=merge(temp, ped,by.x=c("FID","IID"),by.y=c("FID","IID"),sort=F)
	list1=which(is.na(temp2[,kword]))
	if (length(list1)>0){
		temp3=temp2[-which(is.na(temp2[,kword])),]
	} else {temp3=temp2}
	temp3[,3]=temp3[,3]-mean(temp3[,3],na.rm=T)
	for (i in 4:ncol(temp3)){
		temp3[,i]=scale(temp3[,i],center=TRUE,scale=TRUE) ##Covariates and the phenotype
	}
	fit1=lm(paste0(kword,"~."),data=temp3[,-c(1,2)])

	covtemp=as.matrix(cbind(1,temp3[,-which(colnames(temp3)%in% c("FID","IID",kword,"SCORESUM"))]))
	
	fit1$RES_withoutPRS=temp3[,kword]-covtemp%*%fit1$coef[which(names(fit1$coef)!="SCORESUM")]


	return(list("model"=fit1,"data"=temp3))
}




Calculate_PRS<-function(test.bfile,B.beta.info,B.beta.all){
	test.bim=fread(paste0(test.bfile,".bim"),header=F)
	test.bim$V1 <- as.character(sub("^chr", "", test.bim$V1, ignore.case = T))
	test.bim$.index.ref <- 1:nrow(test.bim)

	B.beta.info$.index.tomatch <- 1:nrow(B.beta.info)

	merged <- merge(test.bim,   B.beta.info,all=F, 
                  by.x="V2", by.y="SNP",sort=F)
	###check the sign##
	list1=which(merged$A2==merged$V5)
	if (length(list1)>0){
		B.beta.all[list1,]=-B.beta.all[list1,]
	}

	flag<-rep(FALSE, nrow(test.bim))
	flag[merged$.index.ref]<-TRUE

	BM=as.matrix(B.beta.all[merged$.index.tomatch,])
	BM0=list()
	BM0[[1]]=BM

	system.time({
	pgs <- lapply(BM0, function(x) pgs(bfile=test.bfile, weights = x, 
           extract=flag, keep=NULL, 
           cluster=NULL))
	})

	fam=read.table(paste0(test.bfile,".fam"))[,c(1,2)]
	colnames(fam)[1:2]=c("FID","IID")
	out.all=cbind(fam,pgs[[1]])
	return(out.all)
}


PRStr_main_check<-function(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks){
	out1=0
	if (!file.exists(ped_file)){out1="The ped file does not exist!"} else {
		temp=fread(ped_file,header=T,nrow=1)
		if (!Y_name %in% colnames(temp)){out1="Y does not exist in the ped file!"}
		if (sum(!c("FID","IID") %in% colnames(temp))>0){out1="FID and IID do not exist in the ped file!"}
		if (sum(Covar_name=="")<1){
			if (sum(!Covar_name %in% colnames(temp))>0){out1="The covariates do not exist in the ped file!"}	
		}
	}
	if (!Ytype %in% c("C", "B")) {out1="The Y type is wrong! Now we only support continuous(C) and binary(B) outcomes."}			
	if (file.exists(paste0(train_file,".bim")) & file.exists(paste0(train_file,".bed")) & file.exists(paste0(train_file,".fam"))){} else {out1="The train file doesn't exist!"}
	if (file.exists(paste0(test_file,".bim")) & file.exists(paste0(test_file,".bed")) & file.exists(paste0(test_file,".fam"))){} else {out1="The test file doesn't exist!"}
	if (!file.exists(sum_stats_file)){out1="The summary statistic file does not exist!"} else {
		temp=fread(sum_stats_file,nrow=1)
		if (ncol(temp)==3){
			if (sum(colnames(temp) %in% c("V1","V2","V3"))==3){} else{
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))==3){} else {
					out1="The structure of sum_stats_file is wrong!"
				}
			} 
		} else {
			if (ncol(temp)>3){
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))<3){ 
					out1="The structure of sum_stats_file is wrong!"
				}
			} else {
				out1="The structure of sum_stats_file is wrong!"
			}
		} 
	}
	if (!LDblocks %in% c("EUR.hg19", "AFR.hg19", "ASN.hg19")) {out1="The LDblocks name is wrong!"}
	return(out1)
}


splitgenome2<-function (CHR, POS, ref.CHR, ref.breaks, details = T, right = TRUE)
{
    CHR <- as.character(CHR)
    ref.CHR <- as.character(ref.CHR)
    POS <- as.integer(POS)
    ref.breaks <- as.integer(ref.breaks)
    stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
    stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
    stopifnot(length(CHR) == length(POS))
    stopifnot(length(ref.CHR) == length(ref.breaks))
    chr <- (unique(CHR))
    chr.ref <- (unique(ref.CHR))
    included.chr <- chr %in% chr.ref
    if (!all(included.chr))
        stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
    levels <- character(0)
    results <- character(length(POS))
    Details <- data.frame()
    for (C in chr.ref) {
        breaks <- sort(unique(ref.breaks[ref.CHR == C]))
        if (breaks[1] > 1)
            breaks <- c(1, breaks)
        if (breaks[length(breaks)] < Inf)
            breaks <- c(breaks, Inf)
        cut <- cut(POS[CHR == C], include.lowest = T, breaks = breaks,
            right = right)
        levels <- c(levels, paste0(C, "_", levels(cut)))
        cut <- paste0(C, "_", as.character(cut))
        results[CHR == C] <- cut
        if (details) {
            df <- data.frame(chr = C, start = breaks[-length(breaks)],
                end = breaks[-1])
            Details <- rbind(Details, df)
        }
    }
    results <- factor(results, levels = levels)
    if (details) {
        Details$counts <- as.integer(table(results))
        attr(results, "details") <- Details
    }
   out.item=list(results,Details)
   return(out.item)
}

split_SNPandA1<-function(x){ 
	xtemp=strsplit(x,"_")[[1]]
	if (length(xtemp)==2){
		return(xtemp[1:2])
	} else {
		xn=length(xtemp);
		return(c(paste0(xtemp[1:(xn-1)],collapse="_"),xtemp[xn]))
	}
} 

##cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; num=which(i==unique(LDblocks2[[1]]));nsnp=nrow(bim_sum_stats)
block_calculation<-function(cor,num,train_file,obj,nsnp,temp.file){
  temp_file=paste0(temp.file,"_block_",num)
  write.table(cor$V2,file=temp_file,col.names=F,row.names=F,quote=F)
  cmd = paste0("plink-1.9 --bfile ",train_file," --extract ",temp_file,   " --recodeA  --out ", temp_file,"_Geno.txt")
  system(cmd)
  
  Gtemp=try(as.data.frame(fread(paste0(temp_file,"_Geno.txt.raw"),header=T)),silent=T)
  if (file.exists(temp_file)) {file.remove(temp_file)}
  if (file.exists(paste0(temp_file,"_Geno.txt.nosex"))) {file.remove(paste0(temp_file,"_Geno.txt.nosex"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.log"))) {file.remove(paste0(temp_file,"_Geno.txt.log"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.raw"))) {file.remove(paste0(temp_file,"_Geno.txt.raw"))}

  if (class(Gtemp)=="try-error"){return(NULL)}

  GID=Gtemp[,c(1:2)];GID$order=1:nrow(GID)
  GID2=merge(obj$data[,c("FID","IID")],GID,by=c("FID","IID"),sort=F)
  
  geno=as.matrix(Gtemp[,c(7:ncol(Gtemp))])
  ##as.data.frame(t(sapply(colnames(Gtemp)[-(1:6)], split_SNPandA1  ) ))
  geno_info=as.data.frame(t(sapply(colnames(Gtemp)[-(1:6)],   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
  geno_mean=colMeans(geno,na.rm=T);geno_info$mean=geno_mean; geno_info$sd=NA
  for (j in 1:ncol(geno)){	
    flag=which(is.na(geno[,j]))
    if (length(flag)>0){ geno[flag,j]=geno_mean[j]}
    geno_info$sd[j]=sd(geno[,j])
    if (geno_info$sd[j]>0){
      geno[,j]=(geno[,j]-geno_info$mean[j])/geno_info$sd[j]
    }
  }
  list1=which(geno_info$sd==0)
  if (length(list1)>0){
    geno_info=geno_info[-list1,]
    geno=geno[,-list1]
  }
  if (nrow(geno_info)==0){return(NULL)}
  geno=as.matrix(geno)
  GG=t(geno)%*%geno/(nrow(geno)-1)
  geno_info$order=1:nrow(geno_info)
  
  geno_info2=merge(geno_info,cor[,c("V2","V5","Beta2")],by.x="SNP",by.y="V2",sort=F)
  
  
  if (nrow(geno_info2)>0){
    
    flag_nomatch=which(geno_info2$A1 != geno_info2$V5)
    if (length(flag_nomatch)>0){
      geno_info2$Beta2[flag_nomatch]=-geno_info2$Beta2[flag_nomatch]
    }
    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    
    geno2=as.matrix(geno[GID2$order,geno_info2$order])
    
    
    gy=t(geno2)%*%(obj$model$RES_withoutPRS)/(nrow(obj$data)-1)
    
    
    betatemp=geno_info2$Beta2*geno_info2$sd*obj$model$coef["SCORESUM"]
    u0=gy-GG2%*%betatemp
    beta.all=cbind(u0, betatemp)
    for (factor1 in c(1,10,100,1000)){
      k=1
      betatemp=beta.all[,2]
      u0=beta.all[,1]
      while (k<=15){
        ##betanew=c()
        learningrate=1/nsnp*factor1
        if (learningrate>1){learningrate=1}
        ##print(learningrate)
        for (j in 1:length(betatemp)){
          beta_old=betatemp[j]
          betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
          u0=u0-GG2[,j]*(betatemp[j]-beta_old)
        }
        beta.all=cbind(beta.all,betatemp)
        k=k+1
      } 
    }
    geno_info2=cbind(geno_info2,beta.all)
    return(geno_info2)
  } else {
    return(NULL)
  }##merge>0
}##function end

##cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; num=which(i==unique(LDblocks2[[1]]));nsnp=nrow(bim_sum_stats)
block_calculation2<-function(cor,num,train_file,nsnp,temp.file){
  temp_file=paste0(temp.file,"_block_",num)
  write.table(cor$V2,file=temp_file,col.names=F,row.names=F,quote=F)
  cmd = paste0("plink-1.9 --bfile ",train_file," --extract ",temp_file,   " --recodeA  --out ", temp_file,"_Geno.txt")
  system(cmd)
  
  Gtemp=try(as.data.frame(fread(paste0(temp_file,"_Geno.txt.raw"),header=T)),silent=T)
  if (file.exists(temp_file)) {file.remove(temp_file)}
  if (file.exists(paste0(temp_file,"_Geno.txt.nosex"))) {file.remove(paste0(temp_file,"_Geno.txt.nosex"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.log"))) {file.remove(paste0(temp_file,"_Geno.txt.log"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.raw"))) {file.remove(paste0(temp_file,"_Geno.txt.raw"))}

  if (class(Gtemp)=="try-error"){
    return(NULL)
    #GG=diag(nrow(cor));colnames(GG)=paste0(cor$V2,"_",cor$V5) 
    #geno_info=as.data.frame(t(sapply(colnames(GG),   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    #geno_info$mean=NA; geno_info$maf=NA; geno_info$sd=NA
  }else{
    GG=cor(as.matrix(Gtemp[,7:ncol(Gtemp)]))
    geno_info=as.data.frame(t(sapply(colnames(Gtemp)[7:ncol(Gtemp)],   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    geno_info$mean=colMeans(as.matrix(Gtemp[,7:ncol(Gtemp)]),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))

  }

  list1=which(geno_info$sd==0)
  if (length(list1)>0){
    geno_info=geno_info[-list1,]
    GG=GG[-list1,-list1]
  }
  if (nrow(geno_info)==0){
    return(NULL)
  } else {
    geno_info$order=1:nrow(geno_info)
    geno_info2=merge(cor[,c("V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)
    flag_nomatch=which(geno_info2$A1 != geno_info2$V5)
    if (length(flag_nomatch)>0){
      geno_info2$Beta2[flag_nomatch]=-geno_info2$Beta2[flag_nomatch]
      geno_info2$cor[flag_nomatch]=-geno_info2$cor[flag_nomatch]
    }
    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    gy=geno_info2$cor
    betatemp=geno_info2$Beta2*geno_info2$sd
    u0=gy-GG2%*%betatemp
    beta.all=cbind(u0, betatemp)
    for (factor1 in c(1,10,100,1000)){
      k=1
      betatemp=beta.all[,2]
      u0=beta.all[,1]
      while (k<=15){
        ##betanew=c()
        learningrate=1/nsnp*factor1
        if (learningrate>1){learningrate=1}
        ##print(learningrate)
        for (j in 1:length(betatemp)){
          beta_old=betatemp[j]
          betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
          u0=u0-GG2[,j]*(betatemp[j]-beta_old)
        }
        beta.all=cbind(beta.all,betatemp)
        k=k+1
      } 
    }
    geno_info2=cbind(geno_info2,beta.all)
    return(geno_info2)
  } 
}##function end


##temp.file=paste0(tempfile,"_step1")
PRStr_calculation<-function(obj, train_file, sum_stats, LDblocks, cluster=NULL,temp.file){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse=", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- fread(paste0(train_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  ref.bim$order=1:nrow(ref.bim)
  bim_sum_stats=merge(ref.bim, sum_stats,by.x="V2",by.y="SNP",order=F)
  bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  bim_sum_stats$Beta2=NA
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2]}
  
  bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2")]
  
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  
  if(!is.null(LDblocks)) {
      LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  if(is.null(cluster)) {
  	results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
    		block_calculation(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,obj=obj,nsnp=nrow(bim_sum_stats),temp.file)
  	})
  } else {
  	results.list <-  parallel::parLapplyLB(cluster,unique(LDblocks2[[1]]), function(i) {
    		block_calculation(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,obj=obj,nsnp=nrow(bim_sum_stats),temp.file)
  	})
  }
  
  results.list<-do.call("rbind", results.list)

  return(results.list)
}


##PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"))
##temp.file=paste0(tempfile,"_step1")
PRStr_calculation2<-function(sum_stats_target, train_file, sum_stats, LDblocks, cluster=NULL,temp.file){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse=", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- fread(paste0(train_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  ref.bim$order=1:nrow(ref.bim)
  bim_sum_stats=merge(ref.bim, sum_stats_target,by.x="V2",by.y="SNP",order=F)
  bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  bim_sum_stats$Beta2=NA
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2];  bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
  bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2","cor")]
  
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  
  if(!is.null(LDblocks)) {
      LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  if(is.null(cluster)) {
  	results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
    		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,nsnp=nrow(bim_sum_stats),temp.file)
  	})
  } else {
  	results.list <-  parallel::parLapplyLB(cluster,unique(LDblocks2[[1]]), function(i) {
    		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,nsnp=nrow(bim_sum_stats),temp.file)
  	})
  }
  
  results.list<-do.call("rbind", results.list)

  return(results.list)
}





PRStr_tuning<-function(Beta.all, ped, Covar_name,Y_name, Ytype, test_file){

    beta.info=Beta.all[,1:2]
    for (i in 9:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;##print("fdsf")
      }
    }
    
    beta.all=Beta.all[,9:ncol(Beta.all)]/Beta.all$sd
    PRS.all=Calculate_PRS(test_file,beta.info,beta.all)
    
 
    out.item=data.frame("order"=NA,"R2"=NA)
    for (flag in 3:ncol(PRS.all)){
      out.item[flag-2,1]=flag-2
      resulttemp=PRS.all[,c(1,2,flag)]
      colnames(resulttemp)[3]="SCORESUM"
      if (Ytype=="C"){
        out.item[flag-2,2]=as.numeric(linear_result_generator(resulttemp,ped,Covar_name,Y_name))
      } else {
        out.item[flag-2,2]=as.numeric(logistic_result_generator(resulttemp,ped,Covar_name,Y_name))
      }
    }
    flag=which(out.item$R2==max(out.item$R2))
    param_table=data.frame("lr"=c(0,rep(1/nrow(beta.all)*c(1,10,100,1000),each=15)),"iter"=c(0,rep(1:15,4)))   
    out.final=list()
    out.final$best.learning.rate=param_table$lr[flag]
    out.final$best.iteration=param_table$iter[flag]
    out.final$best.beta=cbind(beta.info,"beta"=beta.all[,flag]) 
    out.final$best.PRS=PRS.all[,c(1:3,flag+2)] ; colnames(out.final$best.PRS)[3:4]=c("PRS.NULL","PRS.TL")
    out.final$param_table=param_table 
    return(out.final)
}



######################The function used individual level data for training############### 
##ped_file=ped.file;Covar_name=cov_name;Y_name=kword;train_file=train.bfile;test_file=test.bfile;sum_stats_file=beta.file;LDblocks="EUR.hg19"
##PRS_TransferLearning(ped.file,"",kword, Ytype="C",train.bfile,test.bfile,beta.file,LDblocks="EUR.hg19",tempfile)}
PRS_TransferLearning_ind<-function(ped_file,Covar_name,Y_name, Ytype="C",train_file,test_file,sum_stats_file,LDblocks="EUR.hg19",tempfile,cluster=NULL){

	out1=PRStr_main_check(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks)
	if (out1!=0){stop(out1)}

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==3){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3"))==3){
			colnames(sum_stats)=c("SNP","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=F,row.names=F,quote=F)
	
	ped=data.frame(fread(ped_file,header=T))[,setdiff(c("FID","IID",Covar_name,Y_name),"")]

	obj=calculate_betaPRS(train_file,sum_stats_file,ped,Covar_name,Y_name,paste0(tempfile,"_step0") ) ##need to remove sum_stats_file and plink command later.

	
	beta_list=PRStr_calculation(obj, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"))
	write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)

	out1=PRStr_tuning(beta_list, ped,Covar_name, Y_name, Ytype,test_file)

  	if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.nosex"))) {file.remove(paste0(tempfile,"_step0.train.PRS.nosex"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.log"))) {file.remove(paste0(tempfile,"_step0.train.PRS.log"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.profile"))) {file.remove(paste0(tempfile,"_step0.train.PRS.profile"))}  
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.PRS,file=paste0(tempfile,"_best.PRS.txt"),row.names=F,quote=F,col.names=T)

	return(out1)
}

######################The function used summary statistics for training############### 
##ped_file=ped.file;Covar_name="";Y_name=kword;Ytype="C"; train_file=train.bfile;test_file=test.bfile;sum_stats_file=beta.file;LDblocks="EUR.hg19"
##ped.file,"",kword, Ytype="C",train.bfile,test.bfile,beta.file,target_sumstats_file,LDblocks="EUR.hg19",tempfile
PRS_TransferLearning<-function(ped_file,Covar_name,Y_name, Ytype="C",train_file,test_file,sum_stats_file,target_sumstats_file, LDblocks="EUR.hg19",tempfile,cluster=NULL){

	out1=PRStr_main_check(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks)
	if (out1!=0){stop(out1)}

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==3){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3"))==3){
			colnames(sum_stats)=c("SNP","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=F,row.names=F,quote=F)
	
	ped=data.frame(fread(ped_file,header=T))[,setdiff(c("FID","IID",Covar_name,Y_name),"")]

	##obj=calculate_betaPRS(train_file,sum_stats_file,ped,Covar_name,Y_name,paste0(tempfile,"_step0") ) ##need to remove sum_stats_file and plink command later.

	sum_stats_target=fread(target_sumstats_file)
	sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP",sort=F)
	if (sum(sum_stats_target$p<=1E-320)>0){ sum_stats_target$p[sum_stats_target$p<=1E-320]=1E-320}

	sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
	flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
	if (length(flag)>0){sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]}
	sum_stats_target=sum_stats_target[,c("SNP","A1.x","Beta","cor")];colnames(sum_stats_target)[2]="A1";
	gc()

	beta_list=PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"))
	beta_list=as.data.frame(beta_list[,-c("A1","order")])
	colnames(beta_list)[1:2]=c("SNP","A1")
	write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)

	out1=PRStr_tuning(beta_list, ped,Covar_name, Y_name, Ytype,test_file)

  	if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.nosex"))) {file.remove(paste0(tempfile,"_step0.train.PRS.nosex"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.log"))) {file.remove(paste0(tempfile,"_step0.train.PRS.log"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.profile"))) {file.remove(paste0(tempfile,"_step0.train.PRS.profile"))}  
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.PRS,file=paste0(tempfile,"_best.PRS.txt"),row.names=F,quote=F,col.names=T)

	return(out1)
}
