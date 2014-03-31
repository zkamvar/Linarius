#############################################################################################################################################################################
###                                                                            Linarius.R                                                                                ####
#############################################################################################################################################################################

#Alleles frequency & heterozygocy 
allele.count<- function(xx,ploidy=2) {
plolev<-unique(ploidy)
NULL->alco
NULL->names
for(i in plolev) cbind(alco,apply(as.matrix(xx[ploidy==i,])==0,2,sum),apply(as.matrix(xx[ploidy==i,])==1,2,sum))->alco
for(j in plolev) names<-c(names,paste(j,"x-absence", sep='',collapse=''),paste(j,"x-presence",sep='', collapse=''))
colnames(alco)<-names
alco
}
#############################################################################################################################################################################
#sep ploidy al frec
allele.frec<-function(xx,ploidy=2) {
	plolev<-unique(ploidy)
NULL->alco
NULL->names
NULL->glob
for(i in plolev) cbind(alco,(1-(apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum)))^(1/i))*100)->alco
for(j in plolev) names<-c(names,paste(j,"x-frequence(%)", sep="",collapse=""))
colnames(alco)<-names
for(k in plolev) glob<-c(glob,nrow(xx[ploidy==k,]))
cbind(alco,apply(alco*glob/sum(glob),1,sum))->alco
colnames(alco)<-c(names,"Frequency")
alco
}
#############################################################################################################################################################################
#sep ploidy al frec 
allele.hetero<-function(xx,ploidy=2) {
	plolev<-unique(ploidy)
NULL->alco
NULL->names
for(i in plolev) cbind(alco,100*(1-((1-(apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum)))^(1/i))^i + (apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum))))))->alco
for(j in plolev) names<-c(names,paste(j,"x-heterozygocy(%)", sep="",collapse=""))
colnames(alco)<-names
alco
}
#############################################################################################################################################################################
#generate Fake Data
datagen<-function(frec,ploidy) {
fakedata<-matrix(data = NA, nrow = length(ploidy), ncol = length(frec), byrow = FALSE,  dimnames = NULL)
plolev<-unique(ploidy)
lingen<-function(xx){
	matrix(data = NA, nrow = length(ploidy), ncol = 1, byrow = FALSE,  dimnames = NULL)->allfrec
	 for(i in plolev)
	sample(c(0,1),length(ploidy[ploidy==i]),replace=T,prob=c((1-xx)^i,1-((1-xx)^i)))->allfrec[ploidy==i,]
	allfrec
	}
for(j in 1:length(frec)) lingen(frec[j])-> fakedata[,j]
	fakedata
}
