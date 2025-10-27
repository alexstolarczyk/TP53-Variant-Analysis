# TP53-Variant-Analysis
This is a python script that collects information about TP53 variants and creates graphs to help identify mutational hotspots within the gene.

tp53_variant_analysis.py obtains information from the Ensembl Rest API about every known variant in the TP53 gene. Variant ID's, Consequence Types, and Start Coordinates for TP53 variants are saved into a data frame. The script runs multiple analyses on this data and outputs multiple .csv and .png files.


# Output files: 
TP53_Variant_Data.csv: Variant ID's, Consequence Types, Start Coordinates, and Sorted Consequence Types for each variant in the TP53 gene.

TP53_Variant_Counts.csv: The TP53 gene was split into 50 different regions(bins) and the number of variants that occur in each of those regions was counted. This csv file shows the genomic regions and the number of variants within each region.

TP53_Consequence_Counts.csv: Contains a list of every Variant Consequence category and the count of occurrences in the TP53 Variant dataset.

TP53_Variant_Histogram.png: Histogram that plots the count of variants at 50 different regions of the TP53 gene.

TP53_Variant_Histogram_Categorized.png: Histogram that plots the count of variants at 50 different regions of the TP53 gene. Bars are stacked and color coded based on the consequence type of the variants (Coding, UTR/Untranslated Region, Noncoding, and Splice)

# Running the script:
python3 final_analysis.py

** The script can be used to create data for genes other than TP53 by changing the value of the 'gene' variable. However, the sort_consequences() function was created using the TP53 consequence types. If a different gene is used, some variants may be sorted into the 'other' category. This can be avoided by adding any of these new consequence types to the respective category within the sort_consequences() function.
