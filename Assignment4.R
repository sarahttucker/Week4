#Sarah Tucker - Assignment 4
#R file to cycle through ab1 files in a Data folder, check their quality against BarcodePlatesStats.csv 
#If their quality is good ("Ok = TRUE" in this case), extract the primary sequence and write it to a FASTA file
#---

#To the person grading this, you likely have seqinr attached/loaded.. I  need you to NOT have that for the first
#chunk of code. I re-load in seqinr when I need to write to the fasta file later on.
detach("package:seqinr", unload = TRUE)

#load in libraries 
library(sangerseqR)
library(Biostrings)


#Assign the quality checking file to the DataQuality variable
DataQuality<-read.csv("./Data/BarcodePlateStats.csv")

#assign all ab1 data files within Data folder to DataFiles variable
DataFiles<-list.files("Data", "*.ab1")

#make empty vector to store names of the good quality abif files
givemedeath<-vector(mode = "character", length = sum(DataQuality$Ok == TRUE))

#IF the data quality is good (i.e. if the Ok column in the DataQuality file reads TRUE), 
#store the name of the file/dna sequence into the empty vector givemedeath
m = 1
for (i in 1:length(DataQuality$Ok)){
  if(DataQuality$Ok[i] == TRUE){
    givemedeath[m] <- as.character(DataQuality$Chromatogram[i])
    m = m + 1
  }
}


#Make two empty vectors to store the DNA sequences (FilteredFiles), and the DNA headers/file name (FileHeader)
FilteredFiles <- vector(mode = "character", length = length(givemedeath))
FileHeader <- vector(mode = "character", length = length(givemedeath))
o = 1


for (i in givemedeath){
    ITS<-read.abif(paste0("./Data/",i)) #read file from the data folder; all file names in givemedeath have
    #already been filtered for quality
    ITSseq <- sangerseq(ITS) # Extract the file
    SeqX<-makeBaseCalls(ITSseq) # Call the primary and secondary peaks
    
    PrimSeq<-primarySeq(SeqX, string = TRUE) #Pull out primary sequence only
    
    FileHeader[o]<-paste(i) #store the file name into file header vector
    FilteredFiles[o] <- paste(PrimSeq) #store the primary sequence into vector
    
    o = o + 1 #cycle through
  }

#load in seqinr library to allow more efficient output to a fasta file
library(seqinr)

#output the header and primary sequence into a fasta file
#after fasta file is made, add additional files to the end of the already made fasta file (GoodDnaTime)
for (e in 1:length(FileHeader)) {
  if (e==1){
    #for the first entry, make the fasta file and input the header and dna sequence
    write.fasta(FilteredFiles[e], file.out = "GoodDnaTime", as.string = TRUE, 
                names = paste(FileHeader[e], "; length=", nchar(FilteredFiles[e]),
                              "; type=DNA", sep = ""), open = "w")
  }
  else if (e > 1){
    #for all additional entries, add onto the existing fasta file
    write.fasta(FilteredFiles[e], file.out = "GoodDnaTime", as.string=TRUE,
                names = paste(FileHeader[e], "; length=", nchar(FilteredFiles[e]), 
                              "; type=DNA", sep = ""), open = "a")
  }
}


