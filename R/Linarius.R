#############################################################################################################################################################################
###                                                                            Linarius.R                                                                                ####
#############################################################################################################################################################################

#Alleles frequency & heterozygocy 
allele.count<- function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
for(i in plolev) cbind(alco,apply(as.matrix(xx[ploidy==i,])==0,2,sum),apply(as.matrix(xx[ploidy==i,])==1,2,sum))->alco
for(j in plolev) noms<-c(noms,paste(j,"x-absence", sep='',collapse=''),paste(j,"x-presence",sep='', collapse=''))
colnames(alco)<-noms
return(alco)
}

new.allele.count <- function(xx, ploidy = 2){
	plolev <- sort(unique(ploidy))
	print(plolev)
	alco <- vapply(plolev, count_alleles, matrix(0, nrow = ncol(xx), ncol = 2), xx, ploidy)
	noms <- as.vector(apply(x, 3, colnames))
	alco <- as.data.frame(alco)
	rownames(alco) <- colnames(xx)
	names(alco) <- noms
	return(alco)
}
#############################################################################################################################################################################
#sep ploidy al frec
allele.frec<-function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
NULL->glob
for(i in plolev) cbind(alco, colMeans(xx[ploidy == i, ]*100) -> alco
for(j in plolev) noms<-c(noms,paste(j,"x-frequence(%)", sep="",collapse=""))
colnames(alco)<-noms
for(k in plolev) glob<-c(glob,nrow(xx[ploidy==k,]))
cbind(alco,apply(alco*glob/sum(glob),1,sum))->alco
colnames(alco)<-c(noms,"Frequency")
return(alco)
}
#############################################################################################################################################################################
#sep ploidy al frec 
allele.hetero<-function(xx,ploidy=2) {
	plolev <- sort(unique(ploidy))
  alco <- NULL
  noms <- NULL
  for (i in plolev){ 
  	alco <- cbind(alco, 100 * (1 - ((1 - (apply(as.matrix(xx[ploidy == 
      i, ]) == 0, 2, sum)/(apply(as.matrix(xx[ploidy == i, 
      ]) == 0, 2, sum) + apply(as.matrix(xx[ploidy == i, ]) == 
      1, 2, sum)))^(1/i))^i + (apply(as.matrix(xx[ploidy == 
      i, ]) == 0, 2, sum)/(apply(as.matrix(xx[ploidy == i, 
      ]) == 0, 2, sum) + apply(as.matrix(xx[ploidy == i, ]) == 
      1, 2, sum))))))
 	}
  for (j in plolev) noms <- c(noms, paste(j, "x-heterozygocy(%)", 
      sep = "", collapse = ""))
  colnames(alco) <- noms
  return(alco)
}

#############################################################################################################################################################################
#generate Fake Data
datagen <- function(frec, ploidy){
	fakedata <- matrix(data = NA, nrow = length(ploidy), ncol = length(frec), 
      							 byrow = FALSE, dimnames = NULL)
  plolev <- unique(ploidy)

  lingen <- function(xx) {
  	allfrec <- matrix(data = NA, nrow = length(ploidy), ncol = 1, 
  					        byrow = FALSE, dimnames = NULL)
    for (i in plolev){ 
    	abs_freq <- (1 - xx)^i
    	probs <- c(abs_freq, 1 - abs_freq)
    	nsamp <- sum(ploidy == i)
    	allfrec[ploidy == i, ] <- sample(0:1, nsamp, replace = TRUE, prob = probs)
    }
    return(allfrec)
  }

  for (j in 1:length(frec)){
  	fakedata[, j] <- lingen(frec[j])
  }
  return(fakedata)
}


count_alleles <- function(ploidy, mat, ploidvec){
	mat <- mat[ploidvec == ploidy, ]
	presence <- colSums(mat)
	absence <- colSums(abs(mat - 1))
	datalist <- matrix(c(absence, presence), ncol = 2)
	noms <- paste0(ploidy, c("x-presence", "x-absence"))
	colnames(datalist) <- noms
	return(datalist)
}