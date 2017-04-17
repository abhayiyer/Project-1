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
  Ordered_Column_Index<- Ordered_Column_Index
  Category_Label <- Category_Label[Ordered_Column_Index]
  Unique_Samples_Name <- Unique_Samples_Name[Ordered_Column_Index]
  readFile[1,2:LengthOfFirstRow] <- Unique_Samples_Name
  for (i in 2:LengthOfFirstCOlumn) {
    readFile[i,2:LengthOfFirstRow] <- readFile[i,Ordered_Column_Index]
    
  }
  data <- readFile
  return(list(Categories=Categories,Category_Label=Category_Label, Unique_Samples_Name = Unique_Samples_Name,GeneList_L=GeneList_L, data=data ))
}

####Read The GeneSet Data
ReadGeneSet <- function(GeneSetData){
  readFile <- read.csv2(GeneSetData, header = FALSE, sep=",")
  NumOfGeneSets <- nrow(readFile)
  NumOfGenes_EachPathway <- vector(length=NumOfGeneSets)
  for(i in 1:NumOfGeneSets){
    NumGene=1
    for(each in readFile[i,]){
      if(each !=""){
        NumGene = NumGene + 1
      }
      
    }
    NumOfGenes_EachPathway[i] <- NumGene
  }
  max_X_dim_ <- max(NumOfGenes_EachPathway)+1
  max_Y_dim_ <- NumOfGeneSets
  Matrix <- matrix(rep(0,max_Y_dim_*max_X_dim_),nrow=max_Y_dim_,ncol = max_X_dim_)
  Matrix[,1]<- as.vector(readFile[,1])
  
}

### This function normalizes the COlumns in the GeneList Data
NormalizeColumnsInArrayData <- function(ArrayData) {
  Colmean <- apply(ArrayData, MARGIN=2, fun=mean)
  ColSD <- apply(ArrayData, MARGIN=2, fun=sd)
  for (i in 1:length(ArrayData[1,])){
    if (ColSD[i]==0) {
      ArrayData[i,]=0
    }
    else{
      ArrayData[,i] <- (ArrayData[,i]-COlMean[i])/ColSD[i]
    }
  }
  return(ArrayData)
}


### GSEA function
GSEA_main <- function(GeneArrayData, Category_Label, Unique_Samples_Name)
{
  NumGenes_ <- length(Unique_Samples_Name)
  NumExperiments_ <- length(GeneArrayData[1,])-1
  print(NumExperiments_)
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
    # randomized_Category_C1[,perm] <- rep(0,NumExperiments_)
    # randomized_Category_C2[,perm] <- rep(0,NumExperiments_)
    # Category_C1[,perm] <- rep(0,NumExperiments_)
    # Category_C2[,perm] <- rep(0,NumExperiments_)
    for(numExp in NumExperiments_)
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
  print(randomized_Category_C1)
  ##Now compute the Signal to Noise for random permutation matrix
}
ReadGeneListInfo <- ReadGeneList("./Leukemia.csv")
# print(ReadGeneListInfo)
# invisible(readline(prompt="Press [enter] to continue"))

# ReadGeneSetInfo <- ReadGeneSet("./pathways.csv")
GSEA_results <- GSEA_main(ReadGeneListInfo$data, ReadGeneListInfo$Category_Label,ReadGeneListInfo$Unique_Samples_Name)