#### Read the GeneListData
ReadGeneList <- function(ArrayData) {
  readFile <- read.csv2(ArrayData, header = FALSE, sep=",")
  Unique_Samples_Name <- colnames(readFile)
  readFile <- as.matrix(readFile)
  LengthOfFirstRow <- length(readFile[1,])
  LengthOfFirstCOlumn <- length(readFile[,1])
  GeneList_L <- readFile[2:LengthOfFirstCOlumn,1]
  Categories <- unique(readFile[1,2:LengthOfFirstRow])
  Category_Label <- vector(length=LengthOfFirstRow-1)
  
  Index = 2
  for(each in readFile[1,2:LengthOfFirstRow]){
    if(each==Categories[1]) {Category_Label[Index-1] <- 0
    Unique_Samples_Name[Index] <- paste(readFile[1,Index],"_",Unique_Samples_Name[Index],sep="")
    Index = Index + 1
    }
    else {Category_Label[Index-1] <- 1
    Unique_Samples_Name[Index] <- paste(readFile[1,Index],"_",Unique_Samples_Name[Index],sep="")
    Index = Index + 1
    }
    
  }
  Ordered_Column_Index = order(Category_Label, decreasing=FALSE)
  Category_Label <- Category_Label[Ordered_Column_Index]
  Ordered_Column_Index<- Ordered_Column_Index+1
  Unique_Samples_Name <- Unique_Samples_Name[Ordered_Column_Index]
  readFile[1,2:LengthOfFirstRow] <- Unique_Samples_Name
  for (i in 2:LengthOfFirstCOlumn) {
    readFile[i,2:LengthOfFirstRow] <- readFile[i,Ordered_Column_Index]
    
  }
  data_ret <- matrix(as.numeric(0),nrow=LengthOfFirstCOlumn,ncol=LengthOfFirstRow)
  for(i in 1:LengthOfFirstCOlumn)
  {
    if(i!=1){data_ret[i,2:LengthOfFirstRow] <- as.numeric(c(readFile[i,2:LengthOfFirstRow]))
    data_ret[i,1]<-readFile[i,1]}
    else{data_ret[i,]<-as.character(c(readFile[i,]))}
  }
  return(list(Categories=Categories,Category_Label=Category_Label, Unique_Samples_Name = Unique_Samples_Name,GeneList_L=GeneList_L, data=data_ret ))
}

RankedList_TwoStates <- function(SortedGeneArrayData, Category_Label, Unique_Samples_Name,GeneList_L)
{
  Numerical_SortedGeneArrayData <- SortedGeneArrayData[2:length(SortedGeneArrayData[,1]),2:length(SortedGeneArrayData[1,])]
  Numerical_SortedGeneArrayData <- apply(Numerical_SortedGeneArrayData,2,as.numeric)
  jpeg(file="Before_Transformation.jpeg")
  heatmap(Numerical_SortedGeneArrayData,Rowv = NA, Colv = NA, col=rainbow(150),scale="row", labRow = FALSE)
  dev.off()
  Split_Index <- split(Category_Label,Category_Label)
  Group_1_Data <- Numerical_SortedGeneArrayData[,1:length(Split_Index[[1]])]
  Group_2_Data <- Numerical_SortedGeneArrayData[,(length(Split_Index[[1]])+1):length(Category_Label)]
  ###Trying Simple Average First
  Group_1_Averaged_Data <- apply(Group_1_Data,1,mean)
  Group_2_Averaged_Data <- apply(Group_2_Data,1,mean)
  Difference_Data <- Group_1_Averaged_Data - Group_2_Averaged_Data
  RankedList <- order(Difference_Data, decreasing = TRUE)
  Ranked_Gene_List <- GeneList_L[RankedList]
  SortedData_BasedOn_RankedGeneList <- Numerical_SortedGeneArrayData[RankedList,]
  jpeg(file="After_AverageDifference_Reranked.jpeg")
  heatmap(SortedData_BasedOn_RankedGeneList, Rowv = NA, Colv = NA, col=rainbow(150),scale="row", labRow = FALSE)
  dev.off()
  return(list(Ranked_Gene_List=Ranked_Gene_List,SortedData_BasedOn_RankedGeneList=SortedData_BasedOn_RankedGeneList))
  
  
}

GSEA_SignalToNoise <- function(GeneArrayData, Category_Label, Unique_Samples_Name)
{
  NumGenes_ <- length(Unique_Samples_Name)
  NumExperiments_ <- length(GeneArrayData[1,])-1
  NumericalData_ <- data.matrix(GeneArrayData[2:length(GeneArrayData[,1]),2:length(GeneArrayData[1,])])
  NumericalData_ <- apply(NumericalData_,2,as.numeric)
  ColMean <- apply(NumericalData_, MARGIN=2, mean)
  ColSD <- apply(NumericalData_, MARGIN=2, sd)
  print(1)
  for (i in 1:length(NumericalData_[1,])){
    if (ColSD[i]==0) {
      NumericalData_[i,]=0
    }
    else{
      NumericalData_[,i] <- (NumericalData_[,i]-ColMean[i])/ColSD[i]
    }
  }  
  subset_mask <- matrix(0,nrow=NumExperiments_,ncol=100)
  randomized_Category_C1 <- matrix(0,nrow=NumExperiments_,ncol=100)
  randomized_Category_C2 <- matrix(0,nrow=NumExperiments_,ncol=100)
  Category_C1 <- matrix(0,nrow=NumExperiments_,ncol=100)
  Category_C2 <- matrix(0,nrow=NumExperiments_,ncol=100)
  C <- split(Category_Label,Category_Label)
  Cat1_ <- C[[1]]
  Cat2_ <- C[[2]]
  Index_Cat1_ <- seq(1,length(Cat1_),1)
  Index_Cat2_ <- seq(length(Cat1_)+1,length(Cat1_)+length(Cat2_),1)
  ###Permutation Step 
  for (perm in 1:100)
  {
    subset_C1_ <- sample(Index_Cat1_,size=length(Cat1_), replace=TRUE)
    subset_C2_ <- sample(Index_Cat2_,size=length(Cat2_), replace=TRUE)
    C1_ <- rep(0,length(Cat1_))
    for(i in 1:length(Cat1_))
    {
      if(is.element(Index_Cat1_[i],subset_C1_)){C1_[i] <- 1}
    }
    C2_ <- rep(0,length(Cat2_))
    for(j in 1:length(Cat2_))
    {
      if(is.element(Index_Cat2_[j],subset_C2_)){C2_[j] <- 1}
    }
    subset_mask[,perm] <- as.numeric(c(C1_,C2_))
    full_subset <- c(subset_C1_,subset_C2_)
    subset_Cat1_ <- sample(full_subset,size=length(Cat1_))
    for(numExp in 1:NumExperiments_)
    {
      x <- sum(!is.na(match(subset_Cat1_,numExp)))
      y <- sum(!is.na(match(full_subset,numExp)))
      randomized_Category_C1[numExp,perm] <- x
      randomized_Category_C2[numExp,perm] <- y-x
      if(numExp <= length(Cat1_))
      {
        Category_C1[numExp,perm] <- y
        Category_C2[numExp,perm] <- 0
      }
      else
      {
        Category_C1[numExp,perm] <- 0
        Category_C2[numExp,perm] <- y        
      }
    }
  }
  ##Now compute the Signal to Noise for random permutation matrix
  Product_1 <- randomized_Category_C1*subset_mask
  sum_1 <- sum(Product_1[,1])
  Matrix_Mul_1 <- NumericalData_%*%Product_1
  Matrix_Mul_1 <- Matrix_Mul_1/sum_1
  Matrix_Mul_1_square <- Matrix_Mul_1*Matrix_Mul_1
  Data_Square_1 <- NumericalData_*NumericalData_
  Sub_1_ <- Data_Square_1 %*% Product_1
  Sub_1_ <- Sub_1_/sum_1 - Matrix_Mul_1_square
  Sub_1_ <- sqrt(abs(sum_1/(sum_1 - 1)*Sub_1_))
  
  Product_2 <- randomized_Category_C2*subset_mask
  sum_2 <- sum(Product_2[,1])
  Matrix_Mul_2 <- NumericalData_%*%Product_2
  Matrix_Mul_2 <- Matrix_Mul_2/sum_2
  Matrix_Mul_2_square <- Matrix_Mul_2*Matrix_Mul_2
  Data_Square_2 <- NumericalData_*NumericalData_
  Sub_2_ <- Data_Square_2 %*% Product_2
  Sub_2_ <- Sub_2_/sum_2 - Matrix_Mul_2_square
  Sub_2_ <- sqrt(abs(sum_2/(sum_2 - 1)*Sub_2_))
  
  Signal <- Matrix_Mul_1 - Matrix_Mul_2
  Noise <- Sub_1_ + Sub_2_
  SignalToNoiSe_Matrix_ <- matrix(0,nrow=NumGenes_,ncol=100)
  SignalToNoiSe_Matrix_ <- Signal/Noise
  
  ###Now Compute Signal to Noise for the GeneList As seen in the experiment
  Product_1 <- Category_C1*subset_mask
  sum_1 <- sum(Product_1[,1])
  Matrix_Mul_1 <- NumericalData_%*%Product_1
  Matrix_Mul_1 <- Matrix_Mul_1/sum_1
  Matrix_Mul_1_square <- Matrix_Mul_1*Matrix_Mul_1
  Data_Square_1 <- NumericalData_*NumericalData_
  Sub_1_ <- Data_Square_1 %*% Product_1
  Sub_1_ <- Sub_1_/sum_1 - Matrix_Mul_1_square
  Sub_1_ <- sqrt(abs(sum_1/(sum_1 - 1)*Sub_1_))
  
  Product_2 <- Category_C2*subset_mask
  sum_2 <- sum(Product_2[,1])
  Matrix_Mul_2 <- NumericalData_%*%Product_2
  Matrix_Mul_2 <- Matrix_Mul_2/sum_2
  Matrix_Mul_2_square <- Matrix_Mul_2*Matrix_Mul_2
  Data_Square_2 <- NumericalData_*NumericalData_
  Sub_2_ <- Data_Square_2 %*% Product_2
  Sub_2_ <- Sub_2_/sum_2 - Matrix_Mul_2_square
  Sub_2_ <- sqrt(abs(sum_2/(sum_2 - 1)*Sub_2_))
  
  Signal <- Matrix_Mul_1 - Matrix_Mul_2
  Noise <- Sub_1_ + Sub_2_
  Experimental_SignalToNoiSe_Matrix_ <- matrix(0,nrow=NumGenes_,ncol=100)
  Experimental_SignalToNoiSe_Matrix_ <- Signal/Noise
  
  File <- "./Histogram/Histogram Comparison.jpg"
  dir.create(dirname(File), showWarnings = FALSE)
  
  jpeg(File)
  par(mfrow=c(1,2))
  
  hist(GSEA_S2N$Random_SignalToNoiseRatio,main = "Histogram of the analysis with 100 Random permutations",xlab = "Random Ratio",xlim = c(-1,1))
  hist(GSEA_S2N$Experimental_SignalToNoiSe_Matrix,main = "Histogram of the analysis of the given sample",xlab = "Experimental Ratio",xlim = c(-1,1))
  dev.off()
  return(list(Random_SignalToNoiseRatio = SignalToNoiSe_Matrix_, Experimental_SignalToNoiSe_Matrix = Experimental_SignalToNoiSe_Matrix_))
  
  }

####Read The GeneSet Data
ReadGeneSet <- function(GeneSetData, RankedGeneList){
  readFile <- read.csv2(GeneSetData, header = FALSE, sep=",")
  readFile <- as.matrix(readFile)
  NumOfGeneSets <- nrow(readFile)
  max_norm_res <- vector()
  GeneSet_S_Num <- 1
  RunningIndex <- matrix((0),nrow=1,ncol=2)
  while(GeneSet_S_Num <= NumOfGeneSets)
  {
    GeneSet_S <- readFile[GeneSet_S_Num,-(1:2)]
    Present_Tag <- sign(match(RankedGeneList,GeneSet_S, nomatch=0))
    Absent_Tag <- 1 -Present_Tag
    Correl_ <- rep(1,length(RankedGeneList))
    Correl_ <- abs(Correl_**1)
    Sum_Correl_Tag <- sum(Correl_[Present_Tag==1])
    Norm_Tag_Present <- 1.0/Sum_Correl_Tag
    Norm_Tag_Absent <- 1.0/(length(RankedGeneList)-length(GeneSet_S))
    RES <- cumsum(Present_Tag*Correl_*Norm_Tag_Present - Absent_Tag*Norm_Tag_Absent)
    RES[is.na(RES)] <- 0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > - min.ES) {
      #      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      #      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
    }
    print(paste("Printing Graph No",GeneSet_S_Num,sep=""))
    file=paste("./RES plots/",trimws(readFile[GeneSet_S_Num,1]),".jpeg",sep="")
    dir.create(dirname(file), showWarnings = FALSE)
    
    jpeg(file)
    
    matplot(x=1:length(RankedGeneList),y=RES,type="l",xlab="RankedGeneList")
    abline(v=arg.ES,col="red")
    dev.off()
   
    Norm_Res <- RES/mean(RES)
    print(paste("MAx_Norm_Res:",max(Norm_Res),sep=""))
    max_norm_res[GeneSet_S_Num] <-max(Norm_Res)
    GeneSet_S_Num = GeneSet_S_Num + 1
  }
return(max_norm_res)  
}

################################################################################################
################################################################################################
ReadGeneListInfo <- ReadGeneList("./Leukemia.csv")
RankedList <- RankedList_TwoStates(ReadGeneListInfo$data, ReadGeneListInfo$Category_Label,ReadGeneListInfo$Unique_Samples_Name, ReadGeneListInfo$GeneList_L)
GSEA_S2N <- GSEA_SignalToNoise(ReadGeneListInfo$data, ReadGeneListInfo$Category_Label,ReadGeneListInfo$Unique_Samples_Name) 
ReadGeneSetInfo <- ReadGeneSet("./pathways.csv", RankedList$Ranked_Gene_List)