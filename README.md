This set of scripts are to analyse the expression on hexaploid wheat (but may be useful for other organisms)



Tables used for `plot_genes_report.R`


This tables are loaded with the function `loadGeneInformation(dir=../TablesForExploration)`

## $canonicalTranscripts

File : **CanonicalTranscript.rds**

* **Gene**
* **transcript**
* **Chr**
* **Start**
* **End**
* **mrna_length**
* **geneconf** If the gene is High Confidence or Low confidence (HC/LC)
* **size_cds**
* **exon_no**
* **exon_length**
* **scaled_1per_position**
* **scaled_5per_position**
* **X3UTR_length**
* **X5UTR_length**
* **intron_length** Calculated from `mrna_length-exon_length`
* **chr_group** Calculated from 4th character from Chr
* **genome** Calculated from 5th character from Chr

Only the genes that are present in the ``$meanTpms`` table are retained for the reports. 
The expectation of genes per position is based ONLY on the genes that have ``geneconf=='HC'`` 

## $meanTpms

File: ``MeanTpms.rds``

* **value** mean TPM
* **factor** Tissue/stress/age.
* **gene** 
* **samples** The number of samples used for the average
* **subset** A name for the analysis where this values where obtained. 
* **min_mean_tpm** The minimum TPM value for at least 1% of the sample used to produce the file. 

The file contains the mean TPM per each factor (obtained by replicates) for each gene. 
When the factor is ``all_means_filtered`` the number of samples represent the number of tissues/factors that whent into that category. So when ``samples == 1``  the gene is uniquely expressed. 
## $triads

File: ``Triads.rds``
* **group_id**
* **clust**
* **description**
* **general_description**
* **Central**
* **A.dominant**
* **B.dominant**
* **D.dominant**
* **A.suppressed**
* **B.suppressed**
* **D.suppressed**
* **value**
* **factor**
* **gene**
* **samples**
* **chr_group**
* **triad_sum**
* **normalised_triad**
* **Distance**
* **P.rank**
* **min_triad_sum**s
* **dataset**

## $triadMovement

File: ``TriadMovement.rds``

* **group_id**
* **factor_count**
* **central_total_distance**
* **central_mean_distance**
* **central_max_distance**
* **central_sd_distance**
* **central_max_over_mean**
* **factor_total_distance**
* **factor_mean_distance**
* **factor_max_distance**
* **factor_sd_distance**
* **factor_max_over_mean**
* **category**
* **sum_mean_tpm**
* **total_categories**
* **categories**
* **Central**
* **A.dominant**
* **A.suppressed**
* **B.dominant**
* **B.suppressed**
* **D.dominant**
* **D.suppressed**
* **dataset**

## $gene_universe

File: ``universe_table.csv``

* **gene**
* **dataset**

## $ontologies

File: ``OntologiesForGenes.rds``

* **Gene**
* **ID**
* **ontology**

## $id_names

File: ``id_names_merged.txt``



* **V1**
* **V2**

This table gives a human readable version of the ontologies. 

Beware, this table doesn't have headers. So when R reads it it automatically assign this names. x§

## $WGCNA

File: ``WGCNA_table.csv``
* **Gene**
* **ModuleLabels**
* **ModuleColors**
* **Number_genes_in_module**
* **set**