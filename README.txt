This readme contains descriptions for all scripts and notebooks contained in this repository.
Note that some paths in certain scripts are hard coded, so if you want to run these scripts, you will have to modify
the respective paths. 


1. Repository structure

/origin_rfc/origin_identification contains files associated with identifying the origins described in the article, 
i.e processing of the assembly data, classifier building and further analyses of the respective species/genomes.

/origin_rfc/metagenome analyses contains files associated with building the kraken2 database and using it to search 
the respective metagenomes for the origin species.


2. File descriptions

2.1 /origin_rfc/origin_identification/

find_origins.py: Takes a directory containing context files created by 'genview-visualize' and creates alignment blocks for all sequences containing the respective gene. GEnView is available at https://github.com/EbmeyerSt/GEnView.

summary_to_tabular.py: Takes summary files created by 'find_origins.py' and transforms them to tabular format.

summary_data_annotated_feb23.xlsx: Contains alignment blocks for all searched genes in all assemblies in tabular format. Labels for 'origin' have been added manually for classifier training.

origin_classifier.ipynb: Jupyter notebook containig data processing, feature engineering/selection and classifier training as well as classifier predictions based on 'summary_data_annotated_feb23.xlsx'. 

gani_classification.py: Uses ANIcalculator_v1 to calculate gANI between selected genomes.

visualize_true_origins.py: Takes GEnView output files and file with GEnView accessions for a restricted number of selected genomes and creates comparative genomics visualization (i.e main text figure 2)

2.2 /origin_rfc/metagenome_analysis

2.2.1 Genome selection for kraken2 database construction and testing

extract_new_origins.py: Use specified genview database to extract all origin species containing the respective gene with >=90% sequence similarity.

extract_origin_genus_references.py: Parse NCBI assembly_summary.txt file to extract accession numbers for origin genus assemblies.

filter_library.py/filter_origin_asms.py: Filter kraken2 bacterial standard library and origin species assemblies for contigs annotated as plasmids or shorter than 5000bp.

convert_headers_to_kraken2.py: Rewrite headers from NCBI assembly format to kraken2-compatible format. Subselect status of genomes to include (i.e 'complete', 'reference', 'all').

generate_mock_communities.py: Use custom assembly summary files to create X mock communities from the genomes specified in the file using art_illumina.

run_kraken2_genus_metrics.py: Run specified mock communities against specified kraken2 database and extract classification metrics for new origin species

kraken2_classification_analysis.ipynb: Evaluate kraken2 classification of new origins in mock communities consisting of all species from the same genus, but not the origin. Contains first analysis, analysis of misclassification and addition/removal of genomes followed by second round of analysis.

2.2.2 Metagenome analysis using kraken2 database.

run_kraken2_metagenomes.py: Run kraken2 on selected metagenome files. 

merge_results.py: Collect results from all runs against kraken2 metagenomes and combine them in a single dataframe.

metagenome_analysis.ipynb: Processing, analysis and visualization of origin abundances in metagenome samples 

abundance_all_metagenomes.csv: Contains kraken2-based abundances of new origin species in all analysed metagenome samples.
