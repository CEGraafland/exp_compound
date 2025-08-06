PropagationCompSimple <- function(baysnet, data.c, nodeComp, level.t2m,level.tp, nodesEvidence, level.ev, perm){
 envlist <- list(baysnet = baysnet,
                 data.c = data.c,
                 nodeComp = nodeComp,
                 level.t2m = level.t2m,
                 level.tp = level.tp,
                 nodesEvidence = nodesEvidence,
                 level.ev = level.ev,
                 perm = perm)
list2env(envlist, envir = parent.frame())
   ngrids<- length(nodes(baysnet))/2
  nodesEvents <- c(nodeComp,ngrids+nodeComp)
 # valueEvents <- c(quote(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m]),
                 #quote(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp]))
  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m])"),
                   paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])"))
  
  if(length(nodesEvidence)==1){
    #valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]))
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  } else if (length(nodesEvidence)>1){
   # valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]),
             #          quote(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]]))
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"),
                       paste0("(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]])"))
  }
  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  
  if (length(nodesEvidence) == 1) {
    # For LS:
    str2 <- paste0("(", names(baysnet)[nodesEvidence[1]],"==", valueEvidence[1], ")")
    #FOR LW:
    #str2 <- paste0("(", names(baysnet)[nodesEvidence[1]], "=", valueEvidence, ")")
    #str2 <- valueEvidence
    # str2 <- paste0("(", names(baysnet)[nodesEvidence[1]]," == ", valueEvidence[1], ")")
    probname <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp  ,"|",names(baysnet)[nodesEvidence[1]]," = ", level.ev[1], ")")
  } else if (length(nodesEvidence) > 1){
    proves <- c()
    j <- 1
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0(names(baysnet)[nodesEvidence[j]],"== ", valueEvidence[j])
    }
    
    text <- "(("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],") & (")
    }
    text <- paste0(text,proves[length(nodesEvidence)],"))")
    str2 <- text
    
    
    probname <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp  ,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,names(baysnet)[nodesEvidence[j]],"=",level.ev[j]," & ")
    }
    probname <- paste0(probname,names(baysnet)[nodesEvidence[length(nodesEvidence)]],"=",level.ev[length(nodesEvidence)],")")
    
  }
  
  
  str <- paste0("(", names(baysnet)[nodesEvents],"==", valueEvents, ")", collapse = " & ")
  #str
  
  #cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",",n = 1000000)")
  cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",")")
  
  cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'ls'",")")
  #cmd1
  #cmd3
  
  with <- eval(parse(text = cmd1))
  without <- eval(parse(text = cmd3))
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp ,")")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}
a <- PropagationCompSimple(baysnet = fitted, data.c = disc.df.t2m.tp, 
                      nodeComp = 87, 
                      level.t2m = 4,
                      level.tp = 1, nodesEvidence = c(81),
                      level.ev = c(4), perm = NULL)
a$with
environment(PropagationCompSimple)
PropagationExactGeneralPerm
