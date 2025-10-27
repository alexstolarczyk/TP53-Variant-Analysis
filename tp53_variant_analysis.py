#!/usr/bin/env python3
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Change this to change the analysis of the script
gene = "TP53"

def get_gene_info(gene):
    """
    Returns the Ensembl ID, chromosome, start coordinates, and end coordinates of a gene.

    Input:
    gene (str): Gene ID to search

    Returns: 
    gene_info: dictionary containing 'Ensembl ID', 'Chromosome', 'start', and 'end' coordinates
    for the inputted gene. Received from ensembl Rest API
    """
    url = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens/' + gene

    headers = {'Content-Type': 'application/json'}

    response = requests.get(url, headers = headers)

    if response.status_code == 200:
        data = response.json()
        gene_info = {'Ensembl_ID': data['id'],
                     'Chromosome': data['seq_region_name'],
                     'start': data['start'],
                     'end': data['end']}
        print("get_gene_info completed successfully")
    else:
        print("Request failed with error code", response.status_code)

    return gene_info
    

def get_variant_info(gene_id, gene):
    """
    Searches the Ensembl Rest API for information about variants of the inputted gene and creates a .csv file
    containing infomation about the variants.

    Input:
    gene_id (str): Ensembl ID of a gene to search
    gene (str): The gene symbol of the gene being searched (for example: 'TP53')

    Output files:
    .csv file containing the variant IDs, Consequence Type, and Start Coordinates for each variant identified in the 
    Ensembl API search. File will be named {gene}_Variant_Data.csv

    Returns:
    df: Pandas data frame containing the variant IDs, Consequence Type, and Start Coordinates for each variant identified
    in the Ensembl API search.
    """
    # Searching the Ensembl API for a list of variants that occur in a given gene
    url = 'https://rest.ensembl.org/overlap/id/' + gene_id + '?feature=variation'

    headers = {'Content-Type': 'application/json'}

    response = requests.get(url, headers = headers)

    variant_ids = []
    consequences = []
    variant_coords = []

    # Saving the variant ids, consequence types, and start coordinates to separate list variables
    if response.status_code == 200:
        data = response.json()
        for entry in data:
            variant_ids.append(entry['id'])
            consequences.append(entry['consequence_type'])
            variant_coords.append(entry['start'])
    
    # Converting the list variables with the data into a data frame and csv file
        variant_data = {"Variant_IDs": variant_ids, "Consequence_Type": consequences, "Start_Coordinates": variant_coords}
        df = pd.DataFrame(variant_data)
        df.to_csv(gene + '_Variant_Data.csv', index = False)
        print("get_variant_info completed successfully")
    else:
        print('Request failed with error code', response.status_code)  

    return df      

def make_variant_histogram(df, gene, chromosome, start, end):
    """
    Makes a histogram showing the amount of known variants at different genomic coordinate regions of the gene of interest.

    Inputs:
    df: Pandas Data frame that contains known variants of the gene being studied and the start coordinate location of those variants
    gene (str): The gene being studied
    chromosome (int): the chromosome in which the gene being studied is located
    start (int): The start coordinate of the gene being studied
    end (int): The end coordinate of the gene being studied

    Output files: 
    Variant CSV File: .csv file showing the raw data behind the analysis. File shows 50 bins of start coordinates and the count of variants within each region.
    Variant Histogram: Histogram showing the count of variants across the gene being studied. Files will be named
    {gene}_Variant_Counts.csv and {gene}_Variant_Histogram.png

    """
    # Creating the .csv file
    counts_data = pd.value_counts(df['Start_Coordinates'], bins = 50)
    counts_dataframe = counts_data.to_frame()
    counts_df_final = counts_dataframe.reset_index()
    counts_df_final.columns = ['Start_Coordinates', 'Count']

    counts_df_final.to_csv(gene + "_Variant_Counts.csv", index=False)
 
    # Creating the histogram

    sns.histplot(data=df, x="Start_Coordinates", bins = 50)
    plt.xlabel(f"Genomic Coordinates (from {start} - {end} bp)")
    plt.ylabel("Variant Count")
    plt.title(f"Variant Frequency across {gene} Located on Chromosome {chromosome}")
    plt.savefig(f"{gene}_Variant_Histogram.png", dpi=300, bbox_inches='tight')
 
    print("make_variant_histogram completed successfully")

def make_consequence_csv(df, gene):
    """
    Makes a csv file containing the number of times each variant consequence type (effect the variant will have on function)
    shows up in the data

    Inputs:
    df: pandas dataframe that contains consequence data about variants of a certain gene
    gene (str): The gene ID being studied

    Output File:
    .csv file with the list of variant consequences and counts. File will be named {gene}_Consequence_Counts.csv
    """
    consequences_data = pd.value_counts(df['Consequence_Type'])
    consequences_dataframe = consequences_data.to_frame()
    consequences_dataframe_final = consequences_dataframe.reset_index()
    consequences_dataframe_final.columns = ['Consequence_Type', 'Count']
    consequences_dataframe_final.to_csv(gene + '_Consequence_Counts.csv', index=False)

    print('make_consequence_csv completed successfully')

def sort_consequences(df, gene):
    """
    Sort the variant consequence types into 4 categories for simplicity and better understanding of the data. The four categories
    are Noncoding, UTR (untranslated region), Coding, and Splice, based on the effects and location of the variants. The
    consequence types were taken from data about TP53 variants. Any variant that does not fit a written consequence category
    will be sorted into "Other."

    Inputs:
    df: pandas dataframe containing a list of variants with variant consequence data
    gene: The gene ID being studied

    Output file: 
    .csv file which will be the same data as the inputted data frame except with one additional column containing the
    sorted consequence category.

    returns:
    df: original pandas dataframe except with one additional column containing the sorted consequence category.
    """
    Sorted_Consequences = []

    for consequence in df['Consequence_Type']:
        if consequence in ['intron_variant', 'intergenic_variant',
                           'regulatory_region_variant']:
            Sorted_Consequences.append("Noncoding")
        elif consequence in ['3_prime_UTR_variant', '5_prime_UTR_variant']:
            Sorted_Consequences.append("UTR")
        elif consequence in ['missense_variant', 'frameshift_variant',
                            'coding_sequence_variant', 'synonymous_variant',
                            'stop_gained', 'inframe_deletion',
                            'inframe_insertion', 'stop_lost',
                            'protein_altering_variant', 'start_lost']:
            Sorted_Consequences.append('Coding')
        elif consequence in ['splice_region_variant', 'splice_polypyrimidine_tract_variant',
                             'splice_acceptor_variant', 'splice_donor_variant',
                             'splice_donor_region_variant', 'splice_donor_5th_base_variant']:
           Sorted_Consequences.append('Splice')
    
        else: 
            Sorted_Consequences.append('Other')

    # Add a new column to the data frame with the Sorted Consequence data
    df = df.assign(Sorted_Consequences = Sorted_Consequences)
    
    df.to_csv(gene + '_Variant_Data.csv', index=False)
    print("sort_consequences completed successfully")

    return df

def make_categorized_histogram(df, gene, chromosome, start, end):
    """
    Makes a histogram showing the amount of known variants at different genomic coordinates across a gene while highlighting
    the variant consequence types of each variant (Noncoding, UTR, Coding, or Splice)
    
    inputs:
    df: pandas dataframe that contains variants and their variant consequence types and start coordinates
    gene: the gene ID being studied
    chromosome: the chromosome in which the gene being studied is located
    start: the start coordinate of the gene being studied
    end: the end coordinate of the gene being studied
    
    output file:
    .png file containing a histogram with variant counts across a gene. Variant consequence type is highlighted by color.
    file will be named {gene}_Variant_Histogram_Categorized.png
    """
    sns.histplot(data=df, x="Start_Coordinates", bins=50,
                 hue="Sorted_Consequences", multiple="stack").legend_.set_title('Consequence Type')
    plt.xlabel(f"Genomic Coordinates (from {start} - {end} bp)")
    plt.ylabel("Variant Count")
    plt.title(f"Variant Frequency across {gene} Located on Chromosome {chromosome}")
    plt.savefig(f'{gene}_Variant_Histogram_Categorized.png', dpi=300, bbox_inches='tight')

    print('make_categorized_histogram completed successfully')


# Pipeline

# Get the information for the inputted gene
gene_info = get_gene_info(gene)

# Return a data frame with data about all variants for the inputted gene
df = get_variant_info(gene_info['Ensembl_ID'], gene)

# Create a histogram showing the frequency of variants at different bp positions across the gene
make_variant_histogram(df, gene, gene_info["Chromosome"],
                    gene_info['start'], gene_info['end'])  

# Create a .csv file showing the number of variants in each consequence category
make_consequence_csv(df, gene)

# Create an additional column in the dataframe categorizing the variant consequences into 4 categories
# Noncoding, UTR, Coding, and Splice
df = sort_consequences(df, gene)

# Create histogram showing the frequency of variants with the consequence type highlighted in different colors
make_categorized_histogram(df, gene, gene_info['Chromosome'],
                           gene_info['start'], gene_info['end'])
