#Script to pull out 15kb upstream RNA sequences near UAAP operons
#This can then be run through hmmscan using HMMER to identify Guanidine-I riboswitches using HMM_guanidineI.hmm
#import programs
library(rentrez)
library(ape)
library(Biostrings)
library(stringr)
#set entrez key
set_entrez_key("entrez_key")
t=1
#load UAAP accession list
id <- read.csv("UAAP1_refseq_acc.txt",header=F)
acc_id=matrix(id[1,],ncol=1)
#retrieve genome accession number and positions
while(t<=nrow(id)){
GI <- entrez_link(dbfrom="protein", id = id[t,], db="nuccore",cmd="neighbor")
if(is.character(GI$links$protein_nuccore)==FALSE){
  t=t+1
}else{
acc <- entrez_fetch(id=GI$links$protein_nuccore, db='nuccore',rettype ='acc')
acc <- unlist(strsplit(acc,"\n"))[1]
fasta <- entrez_fetch(id=GI$links$protein_nuccore, db='nuccore',rettype ='fasta_cds_na')
tab <- write.table(fasta,file="fasta.fasta",sep='\n',quote=FALSE, col.names=FALSE,row.names = FALSE)
tab<- read.FASTA(file='fasta.fasta') 
cds = grep(paste("protein_id=",id[t,],sep=""),names(tab))
if(isEmpty(cds)==TRUE){
  t=t+1
}else{
locate_feature = grep("location=",unlist(strsplit(names(tab)[cds],"[]]")))
location<-unlist(strsplit(names(tab)[cds],"[]]"))[locate_feature]
genepos<-unlist(strsplit(location,"location="))[2]
if(grepl("join",genepos)==TRUE){
  t=t+1
}else{
if(grepl(">",genepos)|grepl("<",genepos)==TRUE){
  genepos<-str_replace(genepos, "([>])", "")
  genepos<-str_replace(genepos, "([<])", "")
}
if(grepl("complement",genepos)==TRUE){
  genestart<-as.numeric(unlist(strsplit(genepos,"[(.)]"))[4])
  geneend<-as.numeric(unlist(strsplit(genepos,"[(.)]"))[2])
  seq_end=toString(genestart+15000)
  seq_beg=toString(genestart)
} else{
genestart<-as.numeric(unlist(strsplit(genepos,"[..]"))[1])
geneend<-as.numeric(unlist(strsplit(genepos,"[..]"))[3])
seq_beg=toString(genestart-15000)
seq_end=toString(genestart)}

#extract fasta sequences 15kb in direction opposite to the direction UAAP is transcribed
upstream <- entrez_fetch(db="nuccore",id=acc,seq_start=seq_beg,seq_stop=seq_end,rettype="fasta")
upstream <- write.table(upstream,file="upstream.fasta",sep='\n',quote=FALSE, col.names=FALSE,row.names = FALSE)
upstream_dna<-readDNAStringSet('upstream.fasta')
#Transcribe into RNA and reverse complement
if(genestart > geneend){
  revcomp<-reverseComplement(upstream_dna)
  rna<-RNAStringSet(str_replace_all(toString(revcomp),"T","U"))
  names(rna)<-id[t,]
  writeXStringSet(rna,'upstream_DUF1989_rna.fasta',append=TRUE,format="fasta")
} else{
  rna<-RNAStringSet(str_replace_all(toString(upstream_dna),"T","U"))
  names(rna)<-id[t,]
  writeXStringSet(rna,'upstream_DUF1989_rna.fasta',append=TRUE,format="fasta")
}
#acc_id=rbind(acc_id,matrix(id[t,],ncol=1))
#write.table(acc_id,file='upstream_DUF1989_acc.fasta',sep="\t",col.names = F, row.names = F,quote = F,append=T)
t=t+1}}}
}

#Run through HMM Search
#Done using ubuntu based HMMer