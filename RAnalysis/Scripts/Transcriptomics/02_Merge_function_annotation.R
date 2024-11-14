# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# LOAD PACKAGES :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(dplyr)
library(kableExtra)
library(pander)
library(data.table)
library(stringr)
library(devtools)





# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
setwd("C:/Users/samjg/Documents/Github_repositories/Airradians_CellularMolecular_OA/RAnalysis") # personal computer
setwd("C:/Users/samuel.gurr/Documents/Github_repositories/Airradians_CellularMolecular_OA/RAnalysis") # personal computer







# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# LOAD DATA (with edits) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# count matrix from prepDE.py script 
# NOTE: aligned to the Airradians draft and unannotated genome!
raw.countmatrix  <- read.csv(file="Output/Transcriptomics/raw_count_matrix_all.csv", header=T) %>% 
                             dplyr::rename(transcript_id = X)
raw.countmatrix[is.na(raw.countmatrix)] <- 0 # replace all occurances of NA with 0 in the cell NOT THE WHOLE ROW!
nrow(raw.countmatrix) # 49900 total transcripts


filtered.countmatrix  <- read.csv(file="Output/Transcriptomics/Filtered_count_matrix/filtered_count_matrix_all.csv", header=T) %>% 
                             dplyr::rename(transcript_id = X)
filtered.countmatrix[is.na(filtered.countmatrix)] <- 0 # replace all occurances of NA with 0 in the cell NOT THE WHOLE ROW!
nrow(filtered.countmatrix) # 10001 total transcripts


# Airradians reference file 
seq_ref_Airradians  <- read.csv(file="Data/Transcriptomics/hpc_out/diamond/GCF_041381155.1_Ai_NY_genomic_ID_REFERENCE.csv", header=F) %>% 
                            tidyr::separate(V1, into = c("NCBI_refseq", "transcript_id", "gene_id"), sep = " ") %>% 
                            dplyr::mutate(gene_id = gsub(".*gene-", "", gene_id), # convert gene-LOC138315196 to just LOC138315196
                                          sseqid  = paste0(NCBI_refseq,"_cds",gene_id)) # the diamond sseqid name to get cvirg KEGG and cgig IDs



# Merge the count matricx with the sequence IDs
# NOTE: this is necessary for the assessement below where sseqid is neeed to merge with the diamon blastx results

raw.countmatrix.REFS <- merge(raw.countmatrix, seq_ref_Airradians, by = 'transcript_id') # 47704 - there are some omitted accession IDs



# Get the GO terms fromthe Cvirgnica genome - the bst annoation of a Atlantic mollusc to date
# call the Cvirginica database of protein names and transcript ID calls
Cvirg_seqID      <-  as.data.table(read.delim2(file = "Data/Transcriptomics/metadata/seq_id.txt", header =F)) %>% 
                              `colnames<-`("fullID")
nrow(Cvirg_seqID) # 66625
Cvirg_GOterms    <-  read.csv(file = "Data/Transcriptomics/metadata/Cviginiva_GOterms.csv", header =T) %>% 
                              dplyr::select(c('GeneID','Annotation_GO_ID', 'Length')) %>% 
                              dplyr::group_by(GeneID) %>% # tif you add GO column here we get duplicates, some of the same gene ID calls (of diff length) have GO term or do not, weird!
                              dplyr::summarise(
                                meanLength = mean(Length)) %>% 
                              unique() # there are many redundant rows here
subset(Cvirg_GOterms,duplicated(GeneID)) # no duplicates, BUT need to filter in the GO terms here 
Cvirg_GOterms2 <- merge(Cvirg_GOterms,
                       ( unique(read.csv(file = "Data/Transcriptomics/metadata/Cviginiva_GOterms.csv", header =T) %>% 
                          dplyr::select(c('GeneID','Annotation_GO_ID')) %>% 
                          dplyr::filter(!Annotation_GO_ID == "")) ), 
                       by = 'GeneID')
nrow(Cvirg_GOterms2) #19667


# Now with the AirradiansIDs and the Cvirgnica GO term refernce - load in the blastx we execuate on SEDNA 
# of the Airradians protein database to the Cvirginina and the Cgigas queries to merge their annotation  (GO and KEGG)
# diamond result to obtain accession IDs of annotated genes Cvirg and Cgigas for gene ID, GO, and KEGG ID information 
#(1) Airradians protein database (...pep.fna file) with Cvirginica nucleotide query
blastx_Airr_Cvirg <- as.data.table(read.delim2(file="Data/Transcriptomics/hpc_out/diamond/AirrProDB_CvirgNQuery/blastx_AirrProDB_CvirgNQuery_out", header=F)) %>% 
                              `colnames<-`(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

#(2) Airradians protein database to the Cgigas KEGG query
blastx_Airr_Cgig  <- as.data.table(read.delim2(file="Data/Transcriptomics/hpc_out/diamond/AirrProDB_CgigNKEGGQuery/blastx_AirrProDB_CgigNKEGGNQuery_out", header=F)) %>% 
                              `colnames<-`(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# where do we go from here? We have to estimate which diamond (blast) was best for obtaining 
# gene annotation for the Airradians transcripts - diagnostics below to reveal which obtained most gene relatedness

  


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  WHICH BLASTX IS BETTER SUITED? :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

nrow(raw.countmatrix) # 49900 total unique transcrips calls in A irradians count matrix

# how many unique trnascript IDs of Airradians were covered by oyster blastx(s)?
# Cvirginica
length(unique(blastx_Airr_Cvirg$sseqid)) # 19598 - Airradians transcripts - in blast x Airradiads Prot database  to Cvriginica nucleotide query
(length(unique(blastx_Airr_Cvirg$sseqid)) / nrow(raw.countmatrix.REFS))* 100 # 41.08% of genes!
# C gigas
length(unique(blastx_Airr_Cgig$sseqid)) # 19667 - Airradians transcripts - in Cgigas protein database to Airradians nucleotide query
(length(unique(blastx_Airr_Cgig$sseqid)) / nrow(raw.countmatrix.REFS))* 100 #  41.22% of genes!

# both are very similar, go with the Atlantic sepcies C virginica 

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# PROCEED WTH C VIRGINICA ANNOTATION      :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# call the best hits with highest bitscore/evalue 

# by bitscore (highest is the best hit) use 'which.max'
bybitscore  <- blastx_Airr_Cvirg[,.SD[which.max(bitscore)],by=sseqid] # max bitscore
length(unique(bybitscore$sseqid)) # 19598
length(unique(bybitscore$sseqid))  == length(unique(blastx_Airr_Cvirg$sseqid))# TRUE
nrow(bybitscore %>% dplyr::filter(sseqid %in% raw.countmatrix.REFS$sseqid)) # 19351
mean(as.numeric(bybitscore$bitscore)) # 454.1144
sd(as.numeric(bybitscore$bitscore))/(sqrt(length(bybitscore$bitscore))) # 3.829952
mean(as.numeric(bybitscore$pident)) # 51.13289
sd(as.numeric(bybitscore$pident))/(sqrt(length(bybitscore$pident))) # 0.1140987
# by evalue (lowest is the best hit) - use 'which.min'
byevalue    <- blastx_Airr_Cvirg[,.SD[which.min(evalue)],by=sseqid] # min evalue
length(unique(byevalue$sseqid)) # 19598
length(unique(byevalue$sseqid))  == length(unique(blastx_Airr_Cvirg$sseqid))# TRUE
nrow(byevalue %>% dplyr::filter(sseqid %in% raw.countmatrix.REFS$sseqid)) # 19351

# calla dataframe for the two sseqids of blatx dataframes by e value and bitscore
# what does this do? if only one column output than the two are the exact same,
#  if two than bitscore (highest) and evalue (lowest) call different transcript IDs
head(as.data.table(c(byevalue$sseqid, bybitscore$sseqid)), header=F) # one column  meaning they are the exact same!

# lets go with evalue as the 'gold stnadard' 
# 'byevalue' gives us the Airradians trnascript ID (i.e. evm.model.Contig....' alonside 
# for each of the corresponding C virginica IDs (i.e. XM_....') to obtaitn KEGG and GO annotation based
# on sequence relatedness
head(byevalue)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MASTER SEQ ID FOR CVRIGNICA             :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



# Now lets call the C virginica transcriptome and edit to fit our needs 
# seq ID reference fr C virginica data 
Cvirg_seqID_editted <- as.data.frame(Cvirg_seqID[Cvirg_seqID$fullID %like% "XM_", ]  %>% # call all mRNA samples - accession always starts with XM
  dplyr::mutate(TranscriptID = (str_match(fullID, ">\\s*(.*?)\\s* PREDICTED:")[,2])) %>% # remove excess ID information
  dplyr::mutate(ProteinID = sub('.*Crassostrea virginica ', '',(gsub("\\s\\(LOC.*|\\sLOC111.*", "", perl=TRUE, fullID))) ) %>% # parse out the protein ID
  dplyr::mutate(GeneID = paste('L', (gsub('),.*', '',(gsub(".*\\s\\(L", "", fullID)))), sep = '')) %>%  # parse out the gene ID
  dplyr::select(-fullID)) # remove the full ID

nrow(Cvirg_seqID_editted) # 60201
nrow(Cvirg_GOterms2) # 19667 - only rows with a GO term present, meanLength of all unique gene IDs
Cvirg_seqIDMASTER <- unique(merge(Cvirg_seqID_editted,Cvirg_GOterms2, by = 'GeneID'))
nrow(Cvirg_seqIDMASTER) # 36573

# write csv
write.csv(Cvirg_seqIDMASTER, file = "Data/Transcriptomics/metadata/seq_id_Cvirginica_master.csv", row.names = FALSE)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MERGE THE BEST BLAST HITS WITH ANNOTATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# objective - merge the data for a master metadata spreadsheet omitted of all annotation (for Cvrginica) that is orrelevant 
# to trnascripts that did not hit to the Airradians transcripts - this will be used for gene function analysis of DEGs!

#  read 'Cvirg_seqIDMASTER' output above 
# file contains the Cvirgnica transcript ID, protein ID, gene ID and GO term annotation
Cvirg_seqID      <-  read.csv(file = "Data/Transcriptomics/metadata/seq_id_Cvirginica_master.csv", header =T) %>% 
                              dplyr::rename(Cvirginica_TranscriptID = TranscriptID)
# # lern how many unique A irradians transcript IDs we have in the raw count matrix 
Airr.ID         <- as.data.frame(raw.countmatrix.REFS$transcript_id) %>% 
                         `colnames<-`("sseqid")
nrow(unique(Airr.ID)) == nrow(Airr.ID) # TRUE 
nrow(Airr.ID) # 47704 - the number of transcripts TOTAL in the raw count matrix1
# merge the Cvirginica seIDs (all cvirginica IDs) with the blastx table we made contianing Airradians hits!
Cvirg_ID.bitscore <- merge(Cvirg_seqID, 
                         #(byevalue   %>% 
                          ( bybitscore %>% 
                            dplyr::select(sseqid, qseqid) %>% 
                            `colnames<-`(c("Airradians_TranscriptID", "Cvirginica_TranscriptID"))), by="Cvirginica_TranscriptID",  all=T) %>% 
                  `colnames<-`(c("blastxEval_CvirgTranscriptID", 
                                 "blastxEval_CvirgProteinID",
                                 "blastxEval_CvirgGeneID", 
                                 "meanLength",
                                 "blastxEval_CvirgGOterms",
                                 "sseqid"))


# we can now do a final merge
# here was have all Airradians Transcript IDs that had the highest 
# evalue hit to the Cvirginica protein database
# merged are the protein names, geneID, GOterms from the Cvirginica database
# to facilitate functional analsiss of DEGs in the Airradians data 
Airr_Cvirg_master_seq_ID  <- merge(seq_ref_Airradians,Cvirg_ID.bitscore,by="sseqid") 
# merge2  <- merge(merge1, Cvirg_ID.bitsc,by="Airradians_TranscriptID", all=T)
nrow(Airr_Cvirg_master_seq_ID) # 34988
(nrow(Airr_Cvirg_master_seq_ID) / nrow(raw.countmatrix))*100 # 70.11 % of genes in our count matrix are represented

# write csv
write.csv(Airr_Cvirg_master_seq_ID, file = "Data/Transcriptomics/metadata/seq_id_AirrCvirg_MERGED_master.csv", row.names = FALSE)
# 'Airr_Cvirg_master_seq_ID' IS OUT MAIN TAKEAWAY FROM THIS CHUNK

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MERGE ANNOTATION WITH THE COUNT MATRIX   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# note: in the following script we will upload the 3 cpm read matrix 
# and merge bytranscript ID to the master data frame with the gene annotation 
# however!!!! the 'byevalue' table of Cvirginica database to Airradians blast query results DOES NOT have 
# all the the Airradians transcripts needed for 

# load the filtered count matrix (decided upon 3 CPM in 50% of samples) review the Counts_Filtered.Rmd script 
read_matrix_raw  <- read.csv(file="Output/Transcriptomics/raw_count_matrix_editted.csv", header=T) %>% 
                             dplyr::rename(Airradians_TranscriptID = 'X')
# load the filtered count matrix (decided upon 3 CPM in 50% of samples) review the Counts_Filtered.Rmd script 
read_matrix_3CPM <- read.csv(file="Output/Transcriptomics/Filtered_count_matrix/filtered_count_matrix_3cpm33perc.csv", header=T) %>% 
                             dplyr::rename(Airradians_TranscriptID = 'X') # make sure you name the transcript ID column the SAME NAME as the master 'Airr_Cvirg_master_seq_ID'



# FILTERED COUNT MATRIX - read_matrix_3CPM 
# merge the functional annotation (from above) to the filtered matrix
MERGE_filtered_reads  <- full_join(read_matrix_3CPM, Airr_Cvirg_master_seq_ID, by = "Airradians_TranscriptID")
MASTER_filtered_reads <- (subset(MERGE_filtered_reads, !is.na(S32)))
nrow(MERGE_filtered_reads) # 20601 - there are some duplicates - based on the amount oF GO terms present in the merge
# call the duplicates and get a single to add back into the matrix 
n_occur <- data.frame(table(MASTER_filtered_reads$Airradians_TranscriptID))
n_occur[n_occur$Freq > 1,]
n_occur_dups <- as.data.frame(n_occur %>% dplyr::filter(Freq > 1))
select_dups_to_add <- MASTER_filtered_reads[MASTER_filtered_reads$Airradians_TranscriptID %in% n_occur$Var1[n_occur$Freq > 1],] %>% 
                        group_by(Airradians_TranscriptID) %>%
                        slice(which.max(is.na(blastxEval_CvirgGOterms)))
# call and ommit hte duplicates
MASTER_filtered_reads_omdups <- MASTER_filtered_reads %>% dplyr::filter(!Airradians_TranscriptID %in% n_occur_dups$Var1)
nrow(MASTER_filtered_reads_omdups) # 9429 with duplicates omitted
# add the, now single form duplicates, back into the matrix
MASTER_filtered_reads_2.0 <- rbind(MASTER_filtered_reads_omdups, select_dups_to_add)
nrow(MASTER_filtered_reads_2.0) # 9557
# diagnostics - must be TRUE
nrow(MASTER_filtered_reads_2.0) == nrow(read_matrix_3CPM)
unique(sort(MASTER_filtered_reads_2.0$Airradians_TranscriptID) == sort(read_matrix_3CPM$Airradians_TranscriptID)) # nust be TRUE!!
# write csv
MASTER_filtered_reads <- write.csv(MASTER_filtered_reads_2.0, file = "RAnalysis/Output/Transcriptomics/Filtered_counts_matrix/filter_3cmpm50perc_WITH_ANNOTATION.csv", row.names = FALSE)



# RAW COUNT MATRIX - read_matrix_raw 
# merge the functional annotation (from above) to the filtered matrix
MERGE_raw_reads  <- full_join(read_matrix_raw, Airr_Cvirg_master_seq_ID, by = "Airradians_TranscriptID")
MASTER_raw_reads <- (subset(MERGE_filtered_reads, !is.na(Ai13)))
nrow(MASTER_filtered_reads) # 2553 - there are some duplicates - based on the amount oF GO terms present in the merge
# call the duplicates and get a single to add back into the matrix 
n_occur <- data.frame(table(MASTER_raw_reads$Airradians_TranscriptID))
n_occur[n_occur$Freq > 1,]
n_occur_dups <- as.data.frame(n_occur %>% dplyr::filter(Freq > 1))
select_dups_to_add <- MASTER_raw_reads[MASTER_raw_reads$Airradians_TranscriptID %in% n_occur$Var1[n_occur$Freq > 1],] %>% 
  group_by(Airradians_TranscriptID) %>%
  slice(which.max(is.na(blastxEval_CvirgGOterms)))
# call and ommit hte duplicates
MASTER_raw_reads_omdups <- MASTER_raw_reads %>% dplyr::filter(!Airradians_TranscriptID %in% n_occur_dups$Var1)
nrow(MASTER_raw_reads_omdups) # 26306 with duplicates omitted
# add the, now single form duplicates, back into the matrix
MASTER_raw_reads_2.0 <- rbind(MASTER_raw_reads_omdups, select_dups_to_add)
nrow(MASTER_raw_reads_2.0) # 26595
# diagnostics - must be TRUE
nrow(MASTER_raw_reads_2.0) == nrow(read_matrix_raw)
unique(sort(MASTER_raw_reads_2.0$Airradians_TranscriptID) == sort(read_matrix_raw$Airradians_TranscriptID)) # nust be TRUE!!
# write csv
MASTER_raw_reads <- write.csv(MASTER_raw_reads_2.0, file = "RAnalysis/Output/Transcriptomics/Raw_counts_matrix/raw_count_matrix_WITH_ANNOTATION.csv", row.names = FALSE)
