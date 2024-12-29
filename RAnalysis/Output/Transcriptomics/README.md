## About Transcriptomics output 

### WGCNA


#### subfolders (WGCNA)


	- challenge
	
	*About*: Our limited sequencing budget was alloted to the full reexposure of the low and high pCO2 histories. 
	Contains 30 indivisuals (columns) as 15 per hisory and 5 per pCO2 exposure (to low, moderate, and high pCO2)
	

		- contents (non subfolders of challenge) 

		contains outputs from the WGCNA workflow including pre/post cut clusterTrees (outlier detection), 
		dendrograms, scale topology (thredholding for eigenge), cluster eigengenes (grouping modules), 
		**and lastly!...** the outputs as heatmaps, ModuleMembership and SignificantModules_annotated csv files. Ive output a .R stepBystep 
		to load these files
	
		* Frontoaded_Activated: 
		
		A few significnat modules following out criteria of frontloading and activated were further investigate for 
		GO (gProfiler) and KEGG enrichment (clusterProfiler). Conetnts of this folder include the subsets of the modulmebership file 
		fr the target modules of interest, further trncated for the target genes deemed frontloaded and activated. Pdf files display 
		the target genes of these modules. 
		
		gProfiler - geneID (LOC...) for frontloaded and activated gene sets copied and pasted into the online GO enrichment interface 
		downloaded the cumulative csv files for (1) enriched terms and (2) gene identities within enriched terms, named for each term dubbed TRUE in the #1
		
		clusterProfiler - KEGG pathway analysis in R for the frontloaded and activated gene sets
		
		* KEGG
		
		clusterProfiler run for significnat modules for the challenge experiment, note that this 
		is not as meaningful or a prior as the frontloaded and activated gene sets - focus on those for the manuscipt described above 
		
		* ModuleExpression_Treatment
	
	- cohort
	
	*About*: Each of the three pcO2 histories exposures to *matched* condition, meaning low x low, moderate x moderate, and high x high. 
	Contains 15 individuals (columns) as 5 per history 

	
		- contents (non subfolders of cohort) 

		contains outputs from the WGCNA workflow including pre/post cut clusterTrees (outlier detection), 
		dendrograms, scale topology (thredholding for eigenge), cluster eigengenes (grouping modules), 
		**and lastly!...** the outputs as heatmaps, ModuleMembership and SignificantModules_annotated csv files. Ive output a .R stepBystep 
		to load these files
	
		* gProfiler

		geneID (LOC...) for frontloaded and activated gene sets copied and pasted into the online GO enrichment interface 
		downloaded the cumulative csv files for (1) enriched terms and (2) gene identities within enriched terms, named for each term dubbed TRUE in the #1
		
		* clusterProfiler 
		
		KEGG pathway analysis in R for the frontloaded and activated gene sets
			


### DESeq2

Alternative method for investigating the cohort experiment. Within this folder are all output files 
listing DEGs for three binary models low v mod, low v high, and mod v high models. I pursued overlap 
between these models to identify genes driving response, using venn R packages to visualize 


### Filtered_count_matrices & Raw_count_matrices

	- all: 
	
	each sample, output from the HPC
	
	- challenge
	
	Our limited sequencing budget was alloted to the full reexposure of the low and high pCO2 histories. 
	Contains 30 indivisuals (columns) as 15 per hisory and 5 per pCO2 exposure (to low, moderate, and high pCO2)
	
	
	- cohort
	
	Each of the three pcO2 histories exposures to *matched* condition, meaning low x low, moderate x moderate, and high x high. 
	Contains 15 individuals (columns) as 5 per history 

**About filtered and transformed count matrices**

	- used *'cpm'* in ```edgeR``` to filter the challenge and cohort matrices of low expressed genes in a thredhl of individuals. 
	Our criteria was >1 x 10 ^3 mean reads per gene in 1/2 of the challenge matrix and 1/3 of the cohort matrix **so that in a theorical gene** 
	**not expressed by a pCO2 hisory (N = 3 2 in challenge and 3 in cohort) remains for analysis (e.g. low x low has 0 reads, moderate x moderate and high x high have > 0, gene remains)**
	
	- matrices were variance stabalized transfomed 'vst'
	
	