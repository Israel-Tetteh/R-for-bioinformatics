############################### LESSON 1 ####################################

# Learning material👇
browseURL('https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html')

# We would make use of packages within the bioconductor open source software
browseURL('https://www.bioconductor.org/packages/release/BiocViews.html#___Software')

## Packages needed for this lesson.
 #- Biological Sequences Retrieval and Analysis
if(!require('seqinr',quietly = TRUE)) install.packages('seqinr')
 #- 	Access the Bioconductor Project Package Repository
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.21")
#! To install packages from Bioconductor.
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

## FASTA format : Format for storing biological data(protein & DNA) sequences
#! Example:
      #> A06852 183 residues
      # MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEE
      # PQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSN
      # LNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIAR
      # LC

## Retrieving genome sequence data via NCBI.
browseURL('https://www.ncbi.nlm.nih.gov/') 
 # get the e Dengue DEN-1 virus genome with accession number: NC_001477
 # download as FASTA file using the send to drop down and open using your text editor.

## Read the downloaded file in to R
list.files()
dengueseq <- read.fasta(file = "Dengue_sequence (1).fasta" )[[1]]
dengueseq[1:50] # get first 50 nucleotides.

## Simple statistics
 #- Get total base pairs count.
length(dengueseq)

 #- Base composition in sequence.
table(dengueseq)

 #- GC content iof DNA
GC(seq = dengueseq) * 100 # how would you have computed this manually? in R
 
 #- DNA words
 count(seq = dengueseq ,wordsize = 1) # single nucleotides
 count(seq = dengueseq ,wordsize = 2) # two nucleotides
 # Can you find which two nucleotide has the highest frequency?🙃
 # Write your sequence into a fasta output.
 
 #- Sliding window analysis of GC content.
 setseq <- seq(from = 1 , to = length(dengueseq), by =  2000) # interested sequence length
 
 gc_countstore <- numeric(length = length(setseq))

 for (i in seq_along(setseq)) {
   
   gc_countstore[i] <- GC(seq = dengueseq[setseq[i] : (setseq[i] + 1999)]) 
   
 } 
print(gc_countstore) 
gc_countstore * 100 # express as percentage

# Plot GC content.
plot(x = setseq,y = gc_countstore,type = 'b',xlab = 'Nuceotide Start position',ylab = 'GC content')

# Can you do this in ggplot2?🙃


## Retrieving Genome sequence from R using the seqinr function.
seqinr::choosebank("refseqViruses")
a <- seqinr::query(listname = "Dengue1" ,query = 'AC=NC_001477')
getSequence(a) # Get sequence
getAnnot(a) # Get annotation
getName(a) # Get accession name
seqinr::closebank() # close database


## Pairwise Sequence Alignment.
browseURL('https://www.uniprot.org/') #  curated database which focuses on protein sequences
#! look out for  Q9CD83 & A0PQ23

## Retrieving Uniprot protein sequence using seqinr
seqinr::choosebank()
seqinr::choosebank("swissprot") # select database to query from
 # Extract Q9CD83 = leprosy  
leprae <- query("leprae", "AC=Q9CD83")
 # Get sequence
lepraeseq <- getSequence(leprae)[[1]]
print(lepraeseq)
 # Extract A0PQ23 = Ulcer
ulcerans <- query("ulcerans", "AC=A0PQ23")
 # Get sequence
ulceransseq <- getSequence(ulcerans)[[1]]
print(ulceransseq)
closebank() # close database


## Comparing two sequences using dot plot.
seqinr::dotPlot(seq1 = lepraeseq ,seq2 = ulceransseq)


## Pairwise global alignment of DNA sequences using the Needleman-Wunsch algorithm
  #  A global alignment is an alignment of the full length of two sequences
  #   A local alignment is an alignment of part of one sequence to part of another sequence.

#! Example; Using +2 for match, -2 for gap and -1 for mismatch
#! G A A T T C
#! G A T T - A

# Let's generate a similar matrix using biostrings.
if(!require('Biostrings')) install.packages('Biostrings')
BiocManager::install("pwalign") # install
# Generate substituition matrix
subs_matrix <- pwalign::nucleotideSubstitutionMatrix(match = +2 ,mismatch = -1 ,baseOnly = T)
print(subs_matrix)


