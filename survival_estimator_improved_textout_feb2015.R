library(data.table)
library(survival)

#### data files can be downloaded from cancerbrowser
####https://genome-cancer.ucsc.edu/download/public/HiSeqV2_PANCAN-2015-02-15.tgz
#### TCGA_pancan_normalized_RNAseq_data_ver_feb06_2015.txt = genomic
#### TCGA_pancan_clinical_data_ver_feb06_2015.txt = clinical

setwd("~/Dropbox/Programs/databases")
data2 <- fread("TCGA_pancan_normalized_RNAseq_data_ver_feb06_2015.txt", header=T, sep="\t", stringsAsFactors=F)
genomic<-as.data.frame(data2)

mat<-as.matrix(genomic[,2:length(genomic)])
row.names(mat)<-genomic$Sample
cln<-read.delim("TCGA_pancan_clinical_data_ver_feb06_2015.txt")
patient<-cln[,c(16,2,4,5,14,15,17,19,18,22)]
row.names(patient)<-cln[,1]

story<-"ribo_splice_new2015"

setwd("../sig2survival")

imported_list<-read.delim("import.txt",header=FALSE)
filter<-match(imported_list[,1],row.names(mat))
input_mat<-as.matrix(mat[filter,])

if(dim(imported_list)[1]>1)
{
	# determine high correlates
	tmat<-t(input_mat)
	correl<-as.data.frame(as.table(cor(tmat,method="spearman")))
	correl<-na.omit(correl)
	list<-correl[(correl$Freq>0.4 & correl$Freq!=1),1:2]
	geneset<-unique(list[,1])
	filtered_mat<-na.omit(input_mat[geneset,])
	scale_mat<-scale(t(filtered_mat))
	#scale_mult<-imported_list[match(colnames(scale_mat),imported_list$gene),2]
	#scale_mult[scale_mult>0]<-1
	#scale_mult[scale_mult<0]<-(-1)
	#scale_mat<-sweep(scale_mat,2,scale_mult,"*")
	sample_mean<-data.frame(apply(scale_mat,1,median),scale_mat) 
}

if (dim(imported_list)[1]==1)
{
	geneset<-as.list(imported_list)
	sample_mean<-scale(input_mat)
}

final<-merge(sample_mean,patient,by="row.names")
names(final)[1:2]<-c("sample_name","median")
OS<-na.omit(final[,c(1:(length(final)-10),(length(final)-7),(length(final)-6),(length(final)-3):(length(final)))])
OS<-OS[OS$sample_type!="Solid Tissue Normal",]
OS$X_OS<-as.numeric(OS$X_OS)
OS$X_OS_IND<-as.numeric(OS$X_OS_IND)
RFS<-na.omit(final[,c(1:(length(final)-10),(length(final)-5),(length(final)-4),(length(final)-3):(length(final)))])
RFS<-RFS[RFS$sample_type!="Solid Tissue Normal",]
RFS$X_RFS<-as.numeric(RFS$X_RFS)
RFS$X_RFS_IND<-as.numeric(RFS$X_RFS_IND)

OS_list<-unique(OS$X_cohort)
RFS_list<-unique(RFS$X_cohort)

stats_total<-{}
stats2<-{}
sig<-gsub("[.]","-",names(final)[3:(length(final)-10)])
for (i in 1:length(OS_list))
{	
	OS_loop<-OS[OS$X_cohort==OS_list[i],]	
	for (k in 2:(length(OS)-6))
	{
		novel<-OS_loop
		merged_set3<-{}
		novel$group<- 2
		low<-novel[,k]<=quantile(novel[,k],0.33)
		high<-novel[,k]>=quantile(novel[,k],0.66)
		
		if(quantile(novel[,k],0.33)!=quantile(novel[,k],0.66))
		{
			novel$group[low]<-1
			novel$group[high]<-3
			merged_set3<-novel[novel$group!=2,]

			#generate km data
			kmsurvival<-survfit(Surv(as.numeric(merged_set3$X_OS),merged_set3$X_OS_IND)~merged_set3$group)
			HR<-summary(coxph(Surv(as.numeric(merged_set3$X_OS),merged_set3$X_OS_IND)~merged_set3$group))
			stats2[[k]]<-cbind("OS",names(novel)[k],OS_list[i],as.data.frame(t(HR$logtest)),HR$conf.int,paste(round(quantile(novel[,k],0.3),3),round(quantile(novel[,k],0.66),3),sep="|"),paste(table(novel$group),collapse="|"),sum(table(novel$group)),paste(sig,collapse="|"))
		}
	}
	stats_total[[i]]<-Reduce(rbind,stats2)
}
OS_final<-na.omit(Reduce(rbind,stats_total))
names(OS_final)[c(1,2,3,11,12,13,14)]<-c("data_type","factor","cancer_type","expression(low|high)","counts(low|mid|high)","total","correlated_signature")

stats_total<-{}
stats2<-{}

for (i in 1:length(RFS_list))
{	
	RFS_loop<-RFS[RFS$X_cohort==RFS_list[i],]	
	for (k in 2:(length(RFS)-6))
	{
		novel<-RFS_loop
		merged_set3<-{}
		novel$group<- 2
		low<-novel[,k]<=quantile(novel[,k],0.33)
		high<-novel[,k]>=quantile(novel[,k],0.66)
		
		if(quantile(novel[,k],0.33)!=quantile(novel[,k],0.66))
		{
			novel$group[low]<-1
			novel$group[high]<-3
			merged_set3<-novel[novel$group!=2,]

			#generate km data
			kmsurvival<-survfit(Surv(as.numeric(merged_set3$X_RFS),merged_set3$X_RFS_IND)~merged_set3$group)
			HR<-summary(coxph(Surv(as.numeric(merged_set3$X_RFS),merged_set3$X_RFS_IND)~merged_set3$group))
			stats2[[k]]<-cbind("RFS",names(novel)[k],RFS_list[i],as.data.frame(t(HR$logtest)),HR$conf.int,paste(round(quantile(novel[,k],0.3),3),round(quantile(novel[,k],0.66),3),sep="|"),paste(table(novel$group),collapse="|"),sum(table(novel$group)),paste(sig,collapse="|"))
		}
	}
	stats_total[[i]]<-Reduce(rbind,stats2)
}
RFS_final<-na.omit(Reduce(rbind,stats_total))
names(RFS_final)[c(1,2,3,11,12,13,14)]<-c("data_type","factor","cancer_type","expression(low|high)","counts(low|mid|high)","total","correlated_signature")

overall<-rbind(OS_final,RFS_final)

matrix<-reshape(overall,idvar=c("cancer_type","data_type"),timevar="factor",v.names="exp(coef)",direction="wide")
names(matrix)<-gsub("exp\\(coef\\).","",names(matrix))

pmatrix<-reshape(overall,idvar=c("cancer_type","data_type"),timevar="factor",v.names="pvalue",direction="wide")
names(pmatrix)<-gsub("pvalue.","",names(pmatrix))


ordered_proper<-12+rev(order(apply(matrix[,13:length(matrix)],2,function(x) cor(matrix[,13],x,method="spearman",use="pairwise"))))
fmatrix<-matrix[rev(order(matrix$median)),c(1,2,5,10,11,ordered_proper)]

pordered_proper<-12+rev(order(apply(pmatrix[,13:length(pmatrix)],2,function(x) cor(pmatrix[,13],x,method="spearman",use="pairwise"))))
fpmatrix<-pmatrix[order(pmatrix$median),c(1,2,5,10,11,pordered_proper)]

dir.create("tested_signatures")
setwd("tested_signatures/")

dir.create(story)
setwd(paste(story,"/",sep=""))
write.table(imported_list,paste(story,sep="","_input_list.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(sig,paste(story,sep="","_correl_gene_list.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.csv(final,paste(story,sep="","_survival_table.csv"))
write.csv(overall,paste(story,sep="","_output.csv"))
write.csv(fmatrix,paste(story,sep="","_matrix_output.csv"))
write.csv(fpmatrix,paste(story,sep="","_pvaluematrix_output.csv"))
write.table(data.frame("+",sig),paste(story,sep="","_CBinput.txt"),sep="",quote=F,row.names=F,col.names=F)
setwd("../")
setwd("../")
