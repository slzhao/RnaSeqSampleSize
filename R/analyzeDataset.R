##' analyze_dataset
##' 
##' A function analyze data set
##' 
##'
##' 
##' @param expObj RangedSummarizedExperiment object.
##' @param expObjGroups sample groups. Should be a vector of 0 and 1. 0 as control samples.
##' @param fdrCut FDR cutoff to select differential genes.
##' @param subset RangedSummarizedExperiment object.
##' @param repN Number of replications.
##' @param useAllSamplesAsNegativeControl Logic. If true, will Use all samples in the obj as negative control
##' @return Figures and a list of result data.
##' @export
##' @examples 1
analyze_dataset=function(expObj,expObjGroups=NULL,
												 fdrCut=0.05,subset=0,repN=2,useAllSamplesAsNegativeControl=FALSE) {
	if (is.null(expObjGroups)) {
		message()
		expObjGroups=sample(c(0,1),ncol(expObj),replace=TRUE)
	}
	
	
	if (useAllSamplesAsNegativeControl) {
		expObjNegativeControl=expObj
		expObjGroupsNegativeControl=expObjGroups
	} else {
		expObjNegativeControl=expObj[,which(expObjGroups==0)]
		expObjGroupsNegativeControl=sample(c(0,1),ncol(expObjNegativeControl),replace = TRUE)
	}
	if (ncol(expObjNegativeControl)<=3) {
		stop("At least 4 control samples needed for dataset analysis")
	}
	
	############################
	#running analysis for Negative Control
	############################
	print("Running analysis for Negative Control")
	reSampleDiffNumAll=NULL
	reSampleDiffFcAll=NULL
	for (i in 1:repN) {
		reSampleGroups=sample(seq_along(expObjGroupsNegativeControl),length(expObjGroupsNegativeControl)
													,replace = FALSE)
		reSampleGroups=expObjGroupsNegativeControl[reSampleGroups]
		#print(reSampleGroups)
		
		y <- estimateDisp(expObjNegativeControl,design= model.matrix(~reSampleGroups))
		y$samples$group=reSampleGroups
		et <- exactTest(y)
		results=topTags(et,n=Inf)
		#print(head(results))
		
		#FC change in resampling
		#hist(results$table$logFC)
		#print(quantile(results$table$logFC,c(0.005,0.025,0.975,0.995)))
		reSampleDiffFcAll=c(reSampleDiffFcAll,results$table$logFC)
		
		reSampleDiffNum=length(which(results$table$FDR<=fdrCut))
		reSampleDiffNumAll=c(reSampleDiffNumAll,reSampleDiffNum)
	}
	print("Number of differential genes in test of negative contol:")
	print(summary(reSampleDiffNumAll))
	print("Quantiles of FC in test of negative contol:")
	print(quantile(reSampleDiffFcAll,c(0.005,0.025,0.05,0.95,0.975,0.995)))
	pFCNegative=ggplot(results$table,aes(x=logFC))+geom_histogram()
	pFCNegative=pFCNegative+theme(text = element_text(size=20))
	

	############################
	#running analysis for Positive Control
	############################
	print("Running analysis for Positive Control")
	y <- estimateDisp(expObj,design= model.matrix(~expObjGroups))
	y$samples$group=expObjGroups
	et <- exactTest(y)
	results=topTags(et,n=Inf)
	
	print("Number of differential genes in test of Positive Control:")
	print(length(which(results$table$FDR<=fdrCut)))
	print("Quantiles of FC in test of Positive Control:")
	print(quantile(results$table$logFC,c(0.005,0.025,0.05,0.95,0.975,0.995)))
	pFCPositive=ggplot(results$table,aes(x=logFC))+geom_histogram()
	pFCPositive=pFCPositive+theme(text = element_text(size=20))
	
	#estimate dispersion for control samples only and compare with dispersion from all samples
	controlSampleInd=which(expObjGroups==0)
	yControl <- estimateDisp(SummarizedExperiment::assay(expObj)[,controlSampleInd],group= expObjGroups[controlSampleInd])

	temp=data.frame(DispersionInTreatmentVsControl=y$tagwise.dispersion,DispersionInControlOnly=yControl$tagwise.dispersion)
	p=ggplot(temp,aes(x=DispersionInTreatmentVsControl,y=DispersionInControlOnly))+geom_point()
	pDispersion=p+theme(text = element_text(size=20))
	
	plot(ggpubr::ggarrange(
		pFCNegative,
		pFCPositive,
		pDispersion,
		ncol = 2, nrow = 2
	))
	
	return(invisible(list(
		NegativeDiffNum=reSampleDiffNumAll,
		NegativeFC=reSampleDiffFcAll,
		PositiveResults=results,
		PositiveEdgeR=y,
		ControlSamplesEdgeR=yControl)))
}


