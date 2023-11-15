import sys
import os
import logging
import numpy as np
import pandas as pd

"""Collect information from kraken2 metagenome runs and metadata and merge results into a single dataframe"""

def merge():

    """Parse all metadata and kraken2 result files to obtain a dataframe containing the counts for each investigated origin species"""

    #Paths to metadata/kraken2 output files
    metadata={'HMP':'/storage/stefan/new_origins_august23/metagenomes/metagenome_results/HMP_all/metadata/hmp_HHS_summarized_meta.tsv', \
              'GS':'/storage/shared/metagenomes/global_sewage/metadata/global_sewage_metadata.csv',\
              'Pig_and_poultry':'/storage/stefan/new_origins_august23/metagenomes/pig_and_poultry/metadata/pig_and_poultry_metadata.txt', \
                'bovine_faecal': '/storage/stefan/new_origins_august23/metagenomes/metagenome_results/bovine_fecal/metadata/bovine_faeces_metadata.txt', \
                'daisy_lake':'/storage/stefan/new_origins_august23/metagenomes/daisy_lake/metadata/daisy_lake_metadata.txt', \
                'sudbury_lake':'/storage/stefan/new_origins_august23/metagenomes/sudbury_lake/metadata/sudbury_lake_metadata.txt', \
                'forest_soil':'/storage/stefan/new_origins_august23/metagenomes/forest_soil/metadata/forest_soil_metadata.txt', \
                'aquaculture':'/storage/stefan/new_origins_august23/metagenomes/aquaculture/metadata/aquaculture_metadata.txt', \
                'indian_lake':'/storage/stefan/new_origins_august23/metagenomes/indian_lake/metadata/Indian_lake.txt', \
                'salt_marsh':'/storage/stefan/new_origins_august23/metagenomes/metagenome_results/salt_marsch/metadata/salt_marsh_metadata.txt', \
                'agriculture_soil':'/storage/stefan/new_origins_august23/metagenomes/metagenome_results/finish_agriculture_soil/metadata/agriculture_soil_metadata.txt', \
                'tara_oceans':'/storage/stefan/new_origins_august23/metagenomes/metagenome_results/tara_oceans/metadata/tara_ocean_metadata.txt', \
                'pune':'/storage/stefan/new_origins_august23/metagenomes/pune/metadata/pune_metadata.txt', \
                'desert_soil':'/storage/stefan/new_origins_august23/metagenomes/desert_soil/metadata/desert_soil_metadata.txt'} #For pune river, the metadata are not accesible via MG-RAST. For now, assume everything to be river water, ask Fanny for metadata

    data_dirs={'HMP':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/HMP_all'], \
               'GS':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/global_sewage_PRJEB40798', \
                    '/storage/stefan/new_origins_august23/metagenomes/metagenome_results/global_sewage_PRJEB40816'], \
                'Pig_and_poultry':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/pig_and_poultry'], \
                'bovine_faecal':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/bovine_fecal/'], \
                'daisy_lake':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/daisy_lake'], \
                'sudbury_lake':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/sudbury_lake'], \
                'forest_soil':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/forest_soil'], \
                'aquaculture':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/aquaculture'], \
                'indian_lake':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/indian_lake'], \
                'salt_marsh':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/salt_marsch'], \
                'agriculture_soil':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/finish_agriculture_soil'], \
                'tara_oceans':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/tara_oceans'], \
                'pune':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/pune'], \
                'desert_soil':['/storage/stefan/new_origins_august23/metagenomes/metagenome_results/desert_soil']}
    
    origins=['Morganella morganii', 'Pseudochrobactrum asaccharolyticum', 'Proteus terrae', 'Providencia rettgeri', 'Providencia stuartii', \
             'Kluyvera cryocrescens', 'Atlantibacter hermannii', 'Citrobacter amalonaticus', 'Aeromonas media']
    
    summary_dict={'sample_acc':[], 'n_reads':[], 'n_classified_reads':[], 'source':[], 'source_class':[], 'n_bacterial_reads':[], 'species':[], 'n_origin_reads':[], 
                  'sample_loc':[], 'project':[]}
    
    for k, v in metadata.items():
        logging.info(f'Parsing results for {k}')

        #Get metadata for respective project 
        for path in data_dirs[k]:
            logging.info(f'Parsing from {path}...')

            meta=[l for l in open(metadata[k], 'r')]

            #Parse kraken2 report files
            for f in [f'{path.rstrip("/")}/{f}' for f in os.listdir(path) if f.endswith('_report.txt')]:
                logging.info(f'File: {f}')
                line_dict={l.split('\t')[5].strip():int(l.split('\t')[1].strip()) for l in open(f, 'r') if any(o==l.split('\t')[5].strip() for o in origins) \
                       or l.endswith('unclassified\n') or l.endswith('root\n') or l.endswith('Bacteria\n')}

                for o in origins:
                    summary_dict['sample_acc'].append(os.path.basename(f).split('_')[0])
                    summary_dict['n_reads'].append(line_dict['unclassified']+line_dict['root'])
                    summary_dict['n_classified_reads'].append(line_dict['root'])
                    summary_dict['n_bacterial_reads'].append(line_dict['Bacteria'])
                    summary_dict['species'].append(o)
                    summary_dict['n_origin_reads'].append(line_dict[o])

                    if k=='HMP':
                        summary_dict['project'].append('Human Microbiome Project')
                        try:
                            if 'duplicates_marked' in f:
                                match=[l for l in meta if os.path.basename(f).split('_')[0].split('.')[0] in l.split('\t')[1]][0]
                                summary_dict['sample_acc'][-1]=os.path.basename(f).split('_')[0].split('.')[0]
                            elif 'duplicates_removed' in f:
                                match=[l for l in meta if os.path.basename(f).split('_')[0] in l.split('\t')[1]][0]

                        except Exception as e:
                            logging.exception(f'{e} - FILE: {f}')

                        summary_dict['source'].append(match.split('\t')[4])
                        summary_dict['source_class'].append('Human feces')
                        summary_dict['sample_loc'].append('Unknown')

                    #The negative samples for GS will have to be removed!
                    elif k=='GS':
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(';')[-2]][0]
                        summary_dict['project'].append('Global sewage')
                        summary_dict['source'].append(match.split(';')[14])
                        summary_dict['source_class'].append('Wastewater')
                        summary_dict['sample_loc'].append(match.split(';')[4])
                    
                    elif k=='Pig_and_poultry':
                        summary_dict['project'].append(f'Pig and poultry feces')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append(f'{match.split(",")[24]} {match.split(",")[26]}')
                        summary_dict['source_class'].append('Animal')
                        summary_dict['sample_loc'].append(match.split(',')[31])
                    
                    elif k=='bovine_faecal':
                        summary_dict['project'].append(f'{k}')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Cow feces')
                        summary_dict['source_class'].append('Animal')
                        summary_dict['sample_loc'].append('North America')

                    elif k=='sudbury_lake':
                        summary_dict['project'].append(f'Sudbury lake')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Freshwater/sediment')
                        summary_dict['source_class'].append('Freshwater')
                        summary_dict['sample_loc'].append('North America/Canada')

                    elif k=='daisy_lake':
                        summary_dict['project'].append(f'Daisy Lake')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Freshwater/sediment')
                        summary_dict['source_class'].append('Freshwater')
                        summary_dict['sample_loc'].append('North America/Canada')
                    
                    elif k=='indian_lake':
                        summary_dict['project'].append(f'Indian Lake')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Freshwater/sediment')
                        summary_dict['source_class'].append('Freshwater')
                        summary_dict['sample_loc'].append(match.split(',')[26])
                        
                    elif k=='aquaculture':
                        summary_dict['project'].append(f'Aquaculture')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Aquaculture')
                        summary_dict['source_class'].append('Saltwater/Animal')
                        summary_dict['sample_loc'].append('China')

                    elif k=='forest_soil':
                        summary_dict['project'].append(f'Forest Soil')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Forest soil')
                        summary_dict['source_class'].append('Soil')
                        summary_dict['sample_loc'].append(match.split(',')[27])

                    elif k=='salt_marsh':
                        summary_dict['project'].append(f'Salt marsh')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Salt water/Sediment')
                        summary_dict['source_class'].append('Saltwater')
                        summary_dict['sample_loc'].append('United Kingdom')

                    elif k=='agriculture_soil':
                        summary_dict['project'].append(f'Finnish agriculture soil')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Agricultural soil')
                        summary_dict['source_class'].append('Soil')
                        summary_dict['sample_loc'].append(match.split(',')[30])

                    elif k=='pune':
                        summary_dict['project'].append(f'Pune river')
                        summary_dict['source'].append('Polluted river water')
                        summary_dict['source_class'].append('Freshwater')
                        summary_dict['sample_loc'].append('India')

                    elif k=='desert_soil':
                        summary_dict['project'].append(f'Desert soil')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Desert rock')
                        summary_dict['source_class'].append('Soil')
                        summary_dict['sample_loc'].append('North America')

                    elif k=='tara_oceans':
                        summary_dict['project'].append(f'Tara oceans')
                        match=[l for l in meta if os.path.basename(f).split('_')[0]==l.split(',')[0]][0]
                        summary_dict['source'].append('Ocean water')
                        summary_dict['source_class'].append('Saltwater')
                        summary_dict['sample_loc'].append('Ocean')

                    """
                    else:
                        print('Add random stuff')
                        summary_dict['project'].append('x')
                        summary_dict['source'].append('x')
                        summary_dict['source_class'].append('x')
                        summary_dict['sample_loc'].append('x')
                    """

                    

        #Check that all lists in summary dict are equally long, if not log warning
        len_check=[False if len(summary_dict['sample_acc'])!=len(v) else True for k, v in summary_dict.items()]
        if False in len_check:
            logging.warning('Arrays are of unequal length, something is wrong!')

    for k, v in summary_dict.items():
        print(k, len(v))

    merged=pd.DataFrame.from_dict(summary_dict)

    return merged

def main():

    outdir=sys.argv[1] #target output directory

    #Set logging configurations
    logging.basicConfig(filename=f'{outdir.rstrip("/")}/run_kraken2.log', filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',\
                         level=logging.INFO)

    df=merge()
    df.to_csv(f'{outdir.rstrip()}/abundance_all_metagenomes.csv', sep='\t')
    logging.info('All files parsed!')
    print('Done!')

if __name__=='__main__':
    main()





