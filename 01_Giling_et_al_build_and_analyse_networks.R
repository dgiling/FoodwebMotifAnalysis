
#===========================================
# 01 Construct and analyse Jena foodwebs
#    Giling et al.
#===========================================

# This code constructs the matrices of the probability of interaction
# between species co-occurring on the experimental
# plots under a number of scenarios. It subsequently 
# calculates a selection of network properties for each network.

# The output is a table of mean values and z-scores
# across iterated webs for each food web


# 1. Required libraries and definitions
#---------------------------------------

## libraries
library(igraph)
library(RColorBrewer)
library(tidyverse)
library(plyr)
library(foreach)
library(doParallel)

## definitions
plant.nodes <- c(paste0("T000",1:9),paste0("T00",10:60))
static.nodes <- c("T0062", "T1204", "T1202", "T1203", "T1207", "T0061","T1208")
basal.nodes <- c(plant.nodes, static.nodes) 


# 2. Define functions
#------------------------

# Curveball algoritm (Strona et al. 2014 Nat. Comms.)
curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

# Motif counter (Borrelli 2015 Oikos)
motif_counter <- function(graph.lists) {
  require(igraph)
  if (!is.list(graph.lists)) {
    stop("The input should be a list of graph objects")
  }
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists),
                         ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1",
                              "d4", "d3", "s2", "s3", "d8", "d2", "d1", "d5", "d7", "d6")
  triad.df <- as.data.frame(triad.matrix)
  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3,
                                 s4 = triad.df$s4, s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2,
                                 d3 = triad.df$d3, d4 = triad.df$d4, d5 = triad.df$d5, d6 = triad.df$d6,
                                 d7 = triad.df$d7, d8 = triad.df$d8)
  rownames(motif.data.frame) <- names(graph.lists)
  return(motif.data.frame)
}

# 'Curving' - to apply curve_ball (modified from Borrelli 2015 Oikos)
curving <- function(web, n.rewires, tgs) {
  #mot <- motif_counter(list(graph.adjacency(adjmat)))
  mot <- web.properties(web, tgs) # modified here to calculate a range of properties (inc. motif_counter)
  newmat <- web

  for (i in 1:n.rewires) {
    newmat <-  curve_ball(newmat)
    row.names(newmat) <- row.names(web)
    colnames(newmat) <- colnames(web)
    #m <- motif_counter(list(graph.adjacency(newmat)))
    m <- web.properties(newmat, tgs) # modified here to calculate a range of properties (inc. motif_counter)
    mot <- rbind(mot, m)
  }
  return(mot[-1, ])
}


# web.properties - Function to calculate a set of network metrics
web.properties <- function(web, tgs) { # web is the adjacency matrix, tgs is a vector of trophic groups
 
  dis <- which(rowSums(web)==0 & colSums(web)==0)  # omit disconnected nodes
  if(length(dis)>0){
    web <- web[-dis, -dis]
    tgs <- tgs[-dis]
  }
  
  web.graph <- graph.adjacency(web)# create igraph object
  
  # basic properties - numbers of species, trophic groups etc.
  S = nrow(web) 
  L = sum(web)
  C = L / S^2
  n.plant = sum(tgs == "plant")
  n.static = sum(tgs == "resource")
  n.herb = sum(tgs == "herbivore")
  n.det = sum(tgs == "detritivore")
  n.omn = sum(tgs == "omnivore")
  n.pred = sum(tgs == "predator")
  p.plant = sum(tgs == "plant")/nrow(web) 
  p.static = sum(tgs == "resource")/nrow(web) 
  p.herb = sum(tgs == "herbivore")/nrow(web) 
  p.det = sum(tgs == "detritivore")/nrow(web) 
  p.omn = sum(tgs == "omnivore")/nrow(web) 
  p.pred = sum(tgs == "predator")/nrow(web) 
  
  # median degree in by group (generality)
  cons.gen = median(degree(web.graph, mode="in")[tgs %in% c("detritivore", "herbivore", "omnivore", "predator")])
  herb.gen = ifelse(n.herb>0, median(degree(web.graph, mode="in")[tgs == "herbivore"]), NA)
  det.gen  = ifelse(n.det>0, median(degree(web.graph, mode="in")[tgs == "detritivore"]), NA)
  omn.gen  = ifelse(n.omn>0, median(degree(web.graph, mode="in")[tgs == "omnivore"]), NA)
  pred.gen = ifelse(n.pred>0, median(degree(web.graph, mode="in")[tgs == "predator"]), NA)
  
  # motifs
  motifs <- motif_counter(list(web.graph))

  # return
  return(cbind(as.data.frame(S), as.data.frame(L), as.data.frame(C),
               as.data.frame(n.plant), as.data.frame(n.static), as.data.frame(n.herb), as.data.frame(n.det), as.data.frame(n.omn), as.data.frame(n.pred),
               as.data.frame(p.plant), as.data.frame(p.static), as.data.frame(p.herb), as.data.frame(p.det), as.data.frame(p.omn), as.data.frame(p.pred),
               as.data.frame(cons.gen), as.data.frame(herb.gen), as.data.frame(det.gen), as.data.frame(omn.gen), as.data.frame(pred.gen),
               motifs))
  
}

# What are the properties being calculated? Must match the return of the above function. Used for labelling results tables.
columns <- c("S", "L", "C",
             "n.plant", "n.static", "n.herb", "n.det", "n.omn", "n.pred",
             "p.plant", "p.static", "p.herb", "p.det", "p.omn", "p.pred",
             "cons.gen", "herb.gen", "det.gen", "omn.gen", "pred.gen",
             "s1","s2","s3","s4","s5","d1","d2","d3","d4","d5","d6","d7","d8")




# 3. Import data
#-----------------------

# Set working directory
dir <- "C:/Desktop/Analyses"

# Import food webs into a list
# These csv files are matrices populated with link type codes according to Hines et al. (accepted)
setwd(paste0(dir,"/prec_webs"))
files <- list.files()
files <- files[201:520] # subset to 2010 + 2012 webs

matrix.list <- list()
for (i in 1:length(files)) {
  m <- read.csv(files[i], row.names=1)
  matrix.list[[i]]<- m
}
names(matrix.list)<-gsub(".csv","",files)
n.webs <- length(matrix.list)

# Import node tables into list
# These csv files are tables of each node in the food web with identity and traits
setwd(paste0(dir,"/nodes"))
files <- list.files()
files <- files[201:520] # subset to 2010 + 2012 webs

node.list <- list()
for (i in 1:length(files)) {
  n <- read.csv(files[i], row.names=1)
  node.list[[i]]<- n
}
names(node.list)<-gsub(".csv","",files)

# Define identifiers
plot.list <- substr(names(matrix.list),13,17)
year.list <- substr(names(matrix.list),1,4)
time.list <- substr(names(matrix.list),6,11)

# Import properties of experimental plots
plot.prop <- read.csv(paste0(dir,"/plot_properties.csv"))
plot.prop$sowndiv.f <- as.factor(plot.prop$sowndiv)
  


# 4. Setup function with loops to implement different probability 
#    scenarios and calculate network properties
#------------------------------------------------------------------- 

# foreach web
calculator <- function(i) {

  results <- vector() # setup an empty vector for results
  
  # select network and node list
  prec.web <- matrix.list[[i]]
  nodes <- node.list[[i]]

  # set identifiers
  focal.plot=plot.list[i]
  focal.year=year.list[i]
  focal.time=time.list[i]
  focal.PSR.level=as.numeric(plot.prop$sowndiv.f[which(focal.plot==plot.prop$plot)])
  
  # calculate relative abundances/cover (within groups: pitfall, suction, plant cover)
  nodes.sum <- nodes %>% group_by(type, sample_method) %>%
    dplyr::summarise(count.sum = sum(count), cover.sum = sum(cover))
  nodes <- merge(nodes, nodes.sum, by=c("type", "sample_method"), all.x=T)
  nodes <- mutate(nodes, rel.abund = count / count.sum)
  nodes <- mutate(nodes, rel.cover = cover / cover.sum)
  nodes$relative <- NA
  nodes$relative[nodes$type=="cons"] <- nodes$rel.abund[nodes$type=="cons"]
  nodes$relative[nodes$type=="plant"] <- nodes$rel.cover[nodes$type=="plant"]
  nodes$relative[nodes$type=="static"] <- 1 # static resources always have relative abundance 1
  
  # make sure the network and node list are in the same order (important!)
  prec.web <- prec.web[order(row.names(prec.web)) , order(row.names(prec.web))]
  nodes <- nodes[order(nodes$taxa.id),]
  if (all(names(prec.web) != nodes$taxa.id)) stop("nodes do not match web") # check they are in the same order

  
  ## Apply probability scenarios
  
  # for behaviours of link types 2 & 3
  for (a in scenarios.a) {
    
    # for behaviours of link types 4 & 5
    for (b in scenarios.b) {
      
      # setup temporary results
      temp.results <- data.frame(matrix(ncol = length(columns)*2, nrow = n.iter))
      colnames(temp.results) <- c(columns, paste0(columns,".z"))
      
      # required indices
      plant.rows <- which(nodes$trophic.group == "plant")
      other.rows <- which(nodes$trophic.group != "plant")
      
      # encounter probability N
      #-------------------------
      N <- as.data.frame(nodes$relative %o% nodes$relative) # encounter prob based on relative abundance
      row.names(N) <- nodes$taxa.id
      colnames(N) <- nodes$taxa.id
     
      N[prec.web==0] <- 0
      N[prec.web==1] <- 1 # 1's are assumed to always encounter each other
      
      # encounter behaviour of link cat 2 & 3
      if (a == "fixed") {
        N[prec.web==2] <- 1
        N[prec.web==3] <- 1
      }
      if (a != "fixed") {
        aa <- as.numeric(a)
        N[plant.rows,][prec.web[plant.rows,]==2] <- N[plant.rows,][prec.web[plant.rows,]==2]*aa # plants as resource
        N[other.rows,][prec.web[other.rows,]==2] <- 1 # not plant as resource, always encounter
        N[plant.rows,][prec.web[plant.rows,]==3] <- N[plant.rows,][prec.web[plant.rows,]==3]*aa # plants as resource
        N[other.rows,][prec.web[other.rows,]==3] <- 1 # not plant as resource, always encounter
      }
      
      # encounter behaviour of link cat 4 & 5
      if (b == "fixed") {
        N[prec.web==4]   <- 1
        N[prec.web==5.1] <- 1
        N[prec.web==5.2] <- 1
        N[prec.web==5.3] <- 1
        N[prec.web==5.4] <- 1
      }
      if (b != "fixed") {
        bb <- as.numeric(b)
        
        ## Optional sensitivity analysis (Supplementary Note 2)
        ## Simple reduction in cons-cons encounter probability with increasing habitat complexity:
        #bb <- bb * c(1, 0.9, 0.8, 0.7, 0.6, 0.5)[focal.PSR.level] 
        
        N[prec.web==4]   <- N[prec.web==4]*bb
        N[prec.web==5.1] <- N[prec.web==5.1]*bb
        N[prec.web==5.2] <- N[prec.web==5.2]*bb
        N[prec.web==5.3] <- N[prec.web==5.3]*bb
        N[prec.web==5.4] <- N[prec.web==5.4]*bb
      }
      
      N[N>1] <- 1 # max encounter probability = 1
    
      # trait-based interaction probability J (termed 'T' in methods)
      #--------------------------------------------------------------
      J <- prec.web
      J[J>0] <- 1
      J[prec.web==1] <- 1
      J[prec.web==2] <- 1
      J[prec.web==3] <- 1
      J[prec.web==4] <- 1
      J[prec.web==5.1] <- 1
      J[prec.web==5.2] <- 1
      J[prec.web==5.3] <- 1
      J[prec.web==5.4] <- 1
      
      # interaction probability A
      #---------------------------
      A <- as.matrix(N * J)
       # hist((A[A!=0]))
      
      
      # number of iterations for each probability scenario
      for (j in 1:n.iter) { 
        
        # sample probabilities
        W <- matrix(rbinom(length(A),prob=A,size=1),nrow=nrow(A))
        row.names(W) <- nodes$taxa.id
        colnames(W) <- nodes$taxa.id
        
        ## SUBSET TO CONSUMER SUB-WEB
        ## Uncomment here to remove plant and/or basal nodes
        #rm <- which(row.names(W) %in% basal.nodes)
        #W <- W[-rm, -rm]
        #nodes.sub <- nodes[-rm,]
          ## or
        nodes.sub <- nodes
        
        ## calculate web properties
        prop <- web.properties(web=W, tgs=nodes.sub$trophic.group)
        
        ## null models for each iteration - based on curveball algoritm as implemented by Borrelli 2015
        rewire.df <- curving(web=W, n.rewires=n.rewires, tgs=nodes.sub$trophic.group)
        rewire.mean <- apply(rewire.df, 2, FUN=mean)
        rewire.sd <- apply(rewire.df, 2, FUN=sd)
        zscore <- (prop - rewire.mean) / rewire.sd
        
        # normalise z scores for motifs - modified code from Borrelli 2015
        motif.zs <- which(names(zscore) %in% c("s1","s2","s3","s4","s5","d1","d2","d3","d4","d5","d6","d7","d8"))
        zscore[,motif.zs] <- t(apply(zscore[,motif.zs], 1, function(x) {
           x/sqrt(sum(x[is.finite(as.numeric(x))]^2))
         }))

        # add to temp results
        temp.results[j,] <- c(prop, zscore)
        
      } # end j loop (number of probability iterations)
      
    # add results to table
    means <- apply(temp.results, 2, FUN=mean)
    sds <- apply(temp.results, 2, FUN=sd)
  
    results <- rbind(results, c(focal.plot, focal.year, focal.time, a, b, means, sds))
    
    } # end a loop
    
  } # end b loop

  data.frame(results) # return table
  
} # end foreach



 
# 5. Run analysis using foreach and dopar (for parallel computing)
#--------------------------------------------------------------------

# settings for analysis
n.iter <- 50
n.rewires <- 250
scenarios.a <- c("fixed", 100, 1000, 10000)
scenarios.b <- c("fixed", 100, 1000, 10000)

# run calculator using foreach for parallel
registerDoParallel()
Sys.time()
results <- foreach (i = 1:320, .combine = rbind, .packages=c('tidyverse', 'igraph','plyr')) %dopar% calculator(i)
Sys.time()
stopCluster()

colnames(results) <- c("plot", "year", "time.period", "alpha", "beta", paste0(columns,".mean"), paste0(columns,".z.mean"), paste0(columns,".sd"), paste0(columns,".z.sd"))
write.csv(results, paste0(dir,"/web_properties.csv"))




