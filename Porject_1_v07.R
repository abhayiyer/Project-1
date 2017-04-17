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



####Read The GeneSet Data
ReadGeneSet <- function(GeneSetData, RankedGeneList){
  readFile <- read.csv2(GeneSetData, header = FALSE, sep=",")
  readFile <- as.matrix(readFile)
  NumOfGeneSets <- nrow(readFile)
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
    jpeg(file=paste(trimws(readFile[GeneSet_S_Num,1]),".jpeg",sep=""))
    matplot(x=1:length(RankedGeneList),y=RES,type="l",xlab="RankedGeneList")
    abline(v=arg.ES,col="red")
    dev.off()
    GeneSet_S_Num = GeneSet_S_Num + 1
  }
  
}

################################################################################################
################################################################################################
ReadGeneListInfo <- ReadGeneList("./Leukemia.csv")
RankedList <- RankedList_TwoStates(ReadGeneListInfo$data, ReadGeneListInfo$Category_Label,ReadGeneListInfo$Unique_Samples_Name, ReadGeneListInfo$GeneList_L)
ReadGeneSetInfo <- ReadGeneSet("./pathways.csv", RankedList$Ranked_Gene_List)