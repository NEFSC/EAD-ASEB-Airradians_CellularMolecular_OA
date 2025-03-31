readme.md

# Archive

## What is in this folder? & How is the data organized?

**About**

Raw data for this project must be archived on NCEI. Below are each of the datasets in the this folder and details on their contents. 


### 'Seawater_chemistry_master.csv'

* contains only the seawater chemistry data measured on days 0, 7, and 14 of the experiment period. **These are the data reported as a tbale in the manuscript**

* *Note:* additional data are reported for the summary means of the F1 and F2 lifetime exposures prior to this experiment 
*however* these data are best suited for the manuscript(s) on these data specifically, here we only archive the exposure challenge 
when organismal data (hemolypmh and expression profiling) was recorded/reported

* columns reflect the NCEI archives of the bay scallop project thus far, including experiment metadat (tank ID, treatment ID, date) 
dipping probe measurements in the tank (temperature, salinity), and instrumental chemistry (temp_spec, pH, TA, TCO2). Again, these same columns 
are those that are archived for manuscripts published on this grant thus far (as of early 2025)


### Experiment_metadata_master.csv 

* contains the experiment metadata for destructive sampling of scallops at day 0 and 14. Important details here are the scallop ID (numeric 1-90) and the tank replicate 
in which each scallop was taken. 

* Column **Note** importantly indicates whether rows pertain to the scallop ID and corresponding replicate_tank prior to the experiment exposure *or* when destructively sampled. 
One can see infer the experiment design from this metadata. Scallops were introdoced to these tanks (inthe same OA system and conditions!, no change in condition!)and labelled on April 24th 2023
and were switched around to start the experiment. 


### Hemolymph_data_master.csv

* contains the raw flow cytometry data for (1) SYBRgreen + PI count live and dead, (2) mitosox greent mean FL1 and (3) the monmoer and Jaggregate for JC-10

* includes the same date and scallop ID information as other archive files to relate these data to seawater cehmistry, etc. 