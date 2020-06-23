#'Make predictions about the chemoresistance for PDAC samples based on RNA-seq gene expression levels
#'@param Sample_Expression_Frame A dataframe of normalized gene expression values for each PDAC sample. Rows should represent genes with row names containing the gene names. Columns should represent samples with column names providing the sample names. It is important that the gene expression are normalized for read depth.
#'@export
#'@examples
#'data("CCLE_Share")
#'findSensitivity(CCLE_Share)

findSensitivity<-function(Sample_Expression_Frame){
  #Load Screen data
  data("Panc1_ScreenData")
  data("Bxpc3_ScreenData")

  #Find common gene names
  commonGeneNames<-P1_SAM_Result_Frame_Share$Gene_Name[P1_SAM_Result_Frame_Share$Gene_Name%in%rownames(Sample_Expression_Frame)]
  Sample_Expression_Frame_Filt<-Sample_Expression_Frame[commonGeneNames,]
  message(paste0(length(commonGeneNames), " out of 18998 total genes found from screen"))

  #Compute screen weights for each drug
  B3_Sums<-list(B3_SAM_Result_Frame_Share[,1]+B3_SAM_Result_Frame_Share[,2],
                B3_SAM_Result_Frame_Share[,4]+B3_SAM_Result_Frame_Share[,5],
                B3_SAM_Result_Frame_Share[,7]+B3_SAM_Result_Frame_Share[,8],
                B3_SAM_Result_Frame_Share[,10]+B3_SAM_Result_Frame_Share[,11])

  P1_Sums<-list(P1_SAM_Result_Frame_Share[,1]+P1_SAM_Result_Frame_Share[,2],
                P1_SAM_Result_Frame_Share[,4]+P1_SAM_Result_Frame_Share[,5],
                P1_SAM_Result_Frame_Share[,7]+P1_SAM_Result_Frame_Share[,8],
                P1_SAM_Result_Frame_Share[,10]+P1_SAM_Result_Frame_Share[,11])

  Combined_Min<-do.call(cbind,lapply(1:4,function(x) apply(cbind(B3_Sums[[x]],P1_Sums[[x]]),1,min)))
  Combined_Max<-do.call(cbind,lapply(1:4,function(x) apply(cbind(B3_Sums[[x]],P1_Sums[[x]]),1,max)))
  row.names(Combined_Min)<- P1_SAM_Result_Frame_Share[,13]
  row.names(Combined_Max)<- P1_SAM_Result_Frame_Share[,13]
  Combined_Min[which(Combined_Min<0)]<-0
  Combined_Min<-Combined_Min^2
  Combined_Max[which(Combined_Max>0)]<-0
  Combined_Max<-Combined_Max^2

  #Scale expression data
  Sample_Expression_Frame_Filt<-t(scale(t(Sample_Expression_Frame_Filt)))

  #Compute final drug sensitivity
  Resistance_Scores<-lapply(1:4, function(x) scale(colSums(Sample_Expression_Frame_Filt*Combined_Min[commonGeneNames,x],na.rm=T)-colSums(Sample_Expression_Frame_Filt*Combined_Max[commonGeneNames,x],na.rm=T))[,1])
  names(Resistance_Scores)<-c("Gemcitabine Resistance", "Oxaliplatin Resistance", "Irinotecan resistance", "Five-FU Resistance")
  return(Resistance_Scores)
}

