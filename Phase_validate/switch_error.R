getOptimalPath<- function(dataPhased,dataSolution,nInd,nMark,ploidy)
{
  Allpossible <- permutations(ploidy,ploidy, v=1:ploidy, set=TRUE, repeats.allowed=FALSE)
  totcombi <- dim(Allpossible)[1]
  hapCount<-1
  path_table <- data.frame()
  
  #get the best possible association of true vs inferred haplotypes
  for(indCount in 1:nInd)
  {
    xTempData <- dataPhased[hapCount:(hapCount + ploidy - 1),]
    xTempREAL <- solutionPhase[hapCount:(hapCount + ploidy - 1),]
    heteroMarks <- apply(xTempREAL,2,sum)   
    heteroMarks <- which((heteroMarks !=ploidy)&(heteroMarks !=0))
    nhetero <- length(heteroMarks) 
    
    #get all the possible combinations
    state_table_Data <- data.frame() 
    cost_vec_Data <- c() 
    curMark_Data <- c()
    prevMark_Data <- c()    
    state_vec_Data <- c(1)
    counter_Data <- 1
    
    #collect the states at each marker in the table state_table
    for(markCount in 1:nMark)
    {
      for(infHap in 1:totcombi)
      {
        if(sum(xTempData[Allpossible[infHap,],markCount]== xTempREAL[,markCount])==ploidy)
        {
          state_table_Data <- rbind(state_table_Data,infHap)
          counter_Data <- counter_Data +1
        }       
      }
      state_table_Data <- rbind(state_table_Data,0)
      counter_Data <- counter_Data +1
      state_vec_Data <- c(state_vec_Data,counter_Data)
    }
    #cost at first marker is zero for all states
    state_table_Data$cost <- c(0)
    state_table_Data$prevState <- c("-1")   
    prevMark_Data <- state_table_Data[state_vec_Data[1]:(state_vec_Data[2]-2),1]    
    counter_Data<-length(prevMark_Data) +1
    prevMark_counter_Data <-0     
    cost_val<-0
    
    #START the viterbi algo.  
    for(stcount in 2:nMark)
    {
      curMark_Data <- state_table_Data[state_vec_Data[stcount]:(state_vec_Data[stcount+1]-2),1] 
      counter_Data <- counter_Data+1
      
      #get the state from prev marker with min cost (num of switches) for each state at this marker
      for(i in 1: length(curMark_Data))
      {   
        mincost <- c()
        for(j in 1: length(prevMark_Data))
        {
          mincost <- c(mincost,sum(Allpossible[curMark_Data[i],] != Allpossible[prevMark_Data[j],]))
        }
        state_table_Data$prevState[counter_Data] <- paste(prevMark_Data[which(mincost == min(mincost))],collapse=" ")
        cost_val <- state_table_Data$cost[(prevMark_counter_Data+which.min(mincost))]+ min(mincost)
        state_table_Data$cost[counter_Data] <- cost_val
        counter_Data <- counter_Data+1
      }
      
      prevMark_counter_Data <- prevMark_counter_Data+length(prevMark_Data)+1       
      prevMark_Data <- curMark_Data
    }
    
    #Now traverse backwards from the last marker to get the least-cost path
    last_best_state_Data <- as.numeric(unlist((strsplit(state_table_Data[(which.min(state_table_Data$cost[(prevMark_counter_Data+1):(prevMark_counter_Data+length(prevMark_Data))]))+prevMark_counter_Data,3]," "))))[1]
    best_state_Data <- as.numeric(state_table_Data[(which.min(state_table_Data$cost[(prevMark_counter_Data+1):(prevMark_counter_Data+length(prevMark_Data))]))+prevMark_counter_Data,1])
    path_Data <- c(best_state_Data)
    
    for(i in (nMark-1):1)
    {
      index<-state_vec_Data[i]-1
      cur_states <- as.numeric(state_table_Data[state_vec_Data[i]:(state_vec_Data[i+1]-2),1])
      best_state_Data <- as.numeric(unlist(strsplit(state_table_Data[index+which(cur_states==last_best_state_Data),3], " ")))[1]     
      path_Data <- c(path_Data,last_best_state_Data)
      last_best_state_Data <- best_state_Data
    }
    path_Data <- rev(path_Data)
    path_table <- rbind(path_table,c(path_Data, nhetero))    
    hapCount <- hapCount+ploidy
  }
  
  return(path_table)
}

compute_switchError <- function(path_table=path_table,nInd,nMark)
{
  Allpossible <- permutations(ploidy,ploidy, v=1:ploidy, set=TRUE, repeats.allowed=FALSE)
  totcombi <- dim(Allpossible)[1]
  switch_Vec <- c()
  
  for(i in 1:nInd)
  {
    Switch_Error <- 0
    
    cur_ind <- unlist(path_table[i,])
    
    for(j in 2:nMark)
    {
      if(cur_ind[j-1] != cur_ind[j])
      {
        diff <- sum(Allpossible[cur_ind[j-1],] != Allpossible[cur_ind[j],])
        Switch_Error = Switch_Error + 1
      }
    }
    Switch_Error <- Switch_Error/(path_table[i,nMark+1]-1)
    switch_Vec <- c(switch_Vec,Switch_Error)
  }
  
  colnames(path_table) <- c(1:nMark,"Num_hetero")
  path_table$Switch_error <- switch_Vec 
  
  return(path_table)
}
process_MyprogOutput <- function (fNamemyProg,ploidy)
{
  #read output from my program
  phasedTxtMyprog <- readLines(fNamemyProg,encoding="UTF-8")
  nInd <- as.numeric(phasedTxtMyprog[1])
  nMark <- as.numeric(phasedTxtMyprog[2])
  phasedTxtMyprog <- phasedTxtMyprog[!grepl(pattern = "(^[#])", x = phasedTxtMyprog)]
  phasedTxtMyprog <- phasedTxtMyprog[3:length(phasedTxtMyprog)]
  data_phasedMyprog <- data.frame(phasedTxtMyprog,stringsAsFactors=F)
  data_phasedMyprog <- data.frame(apply(data_phasedMyprog,1,function(x)as.numeric(unlist(strsplit(x," ")))))
  data_phasedMyprog <-t(data_phasedMyprog)
  
  return(data_phasedMyprog)
}

library("gtools", lib.loc="/home/soumya/R/x86_64-suse-linux-gnu-library/3.1")
args <- commandArgs(trailingOnly = TRUE)
fNamemyProg <- args[1]
fNameSol <- args[2]
ploidy <- as.numeric(args[3])
nInd <- as.numeric(args[4])
nMark <- as.numeric(args[5])

numHaplo <- nInd*ploidy
#read phased out from phasing method
data_phasedMyprog <- process_MyprogOutput(fNamemyProg=fNamemyProg,ploidy=ploidy)
#read solution
sol_rawData<-read.table(fNameSol,header=F)
solutionPhase<-t(sol_rawData[1:nMark,3:(2+numHaplo)])
path_table_Myprog <- getOptimalPath(dataPhased=data_phasedMyprog,
                                    dataSolution=solutionPhase,nInd=nInd,nMark=nMark,ploidy=ploidy)
write.csv(path_table_Myprog,paste("switch_error ",nInd,"_m_",nMark,".csv",sep=""))


