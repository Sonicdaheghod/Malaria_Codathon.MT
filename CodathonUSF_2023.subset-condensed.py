#!/usr/bin/env python
# coding: utf-8

# In[21]:


#Credit for code:
#http://alimanfoo.github.io/2018/04/09/selecting-variants.html
#https://scikit-allel.readthedocs.io/en/stable/io.html
#https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
#https://www.geeksforgeeks.org/how-to-convert-numpy-array-to-dictionary-in-python/


# # The goal of this project with the team I am in  (Justin Gibbons, Indu Parameswaran, Shenica Jerome) is to determine a specific variant of a gene that caused patients to contract Malaria from a dataset consisting of patients with Malaria.
# 
# ![image.png](attachment:image.png)
# 

# # 1) Import and read through gzip file

# In[22]:


pip install scikit-allel


# In[23]:


pip install "setuptools<58"


# In[24]:


pip install pyvcf


# # 2) Organize Data 
# 
# ## The subset dataself is a chuck of the whole dataset, makes the code easier to run

# In[25]:


import vcf

vcf_file = vcf.Reader(filename=r"random_sample_n100.vcf.gz")

#we have to iterate through the file to get the SNP from each chromosome point in the dataset

for record in vcf_file:
    # Check if the variant is a SNP (single nucleotide variant) via iteration through dataset
    if len(record.REF) == 1 and len(record.ALT) == 1 and len(record.ALT[0]) == 1:
        print(f"Chromosome: {record.CHROM}")
        print(f"Position: {record.POS}")
        print(f"Reference allele: {record.REF}")
        print(f"Alternate allele: {record.ALT[0]}")
        print()
        


#heterozygous
#convert to binary


# # 3) Mark Data 0 or 1
# 
# ## Mark each data point 1 or 0 where 1 = reference allele ( Plasmodium falciparum) 0 = alternative allele

# In[26]:


#O(n^x): Polynomial complexity notation

for record in vcf_file:
    for sample in record.samples:
        print(sample["GT"])
        
print(record.samples)


# ## 3a) Get different titles from dataset

# In[27]:


import allel
callset = allel.read_vcf(r"random_sample_n100.vcf.gz")
sorted(callset.keys())


# ## 3c) convert 0/1 marking to martix for only GT label

# ##  - filter subtelomericHypervariable
# ## -1 is reference, 1 is alternative
# 
# ## 4) Get homozygous data [0,0], [1,1]

# In[28]:


## here, needed to attach np to array part of code

gt_zarr = callset['calldata/GT']
        
import numpy as np

arrays = gt_zarr

#get homozygous data

for array in arrays:
    if np.array_equal(array, np.array([1, 1])) or np.array_equal(array, np.array([0, 0])):
        print(array)
        
arrays


# # 5) Graph Data

# ## 5a) Create Dictionary

# In[29]:


#array into dictionary
#key:value
#key is the sample name, value is the genotype

#I define sample and GT as the callset function of respective titles

keys = callset["samples"]
vals = callset['calldata/GT']

Malaria_dict = dict(zip(keys, zip(*vals)))
list(Malaria_dict.items())


# ## 5b) Encode Genotype to -1,1,or N/A

# ### (negative)-1 for Reference homozygous genotype (0,0), 1 for Alternative homozygous genotype (1,1), N/A for Heterozygous

# In[30]:


import allel
import numpy as np

# Load the VCF file
callset = allel.read_vcf("random_sample_n100.vcf.gz")

# Access the genotype data
genotypes = callset['calldata/GT']

# Convert the genotypes array into a set
genotypes_set = set(tuple(map(tuple, gt)) for gt in genotypes)

# Print the first 10 values in the set
print(list(genotypes_set))


# In[31]:


# Iterate over the genotypes
def iteration():
    for record in genotypes:
        for sample in record:
            print(sample)

iteration()


# In[32]:


# Define the genotype-value mappings
genotype_value_mapping = {
    (0, 0): -1,
    (1, 1): 1,
    (0, 1): "N/A",
    (1, 0): "N/A",
    (-1, -1): "N/A"  
}
# Define GT variable
GT = callset['calldata/GT']
GT
# Create an empty dictionary
genotype_dict = {}

# Iterate over the genotypes and update the dictionary
for i, sample_name in enumerate(keys):
    sample_genotypes = []
    for genotype in GT[i]:
        value = genotype_value_mapping.get(tuple(genotype), "NA")
        sample_genotypes.append(value)
    genotype_dict[sample_name] = sample_genotypes

# Print the genotype dictionary

print("Sample and Genotype as Dictionary:")

for sample, genotypes in genotype_dict.items():
    print(f"{sample}, {genotypes}")


# ## 6) Organize Dictionary into Pandas Dataframe

# In[33]:


#top label = sample name
#vertical axis = labelling genotypes

import pandas as pd

df = pd.DataFrame(genotype_dict)
df


# In[34]:


#add column for chromosome positions 

position = callset["variants/POS"]
position

df_position = pd.DataFrame(position, columns=["Chromosome Position"])
df_position



# In[35]:


#Chromosome and Position 
chromosome = callset["variants/CHROM"]
chromosome
df_chromosome = pd.DataFrame(chromosome,columns=["Chromosome name"])

df_chromosome


# In[36]:


#concatenate chromosome, position, samples and the genotype value marking into single pd dataframe

frames = [df_chromosome,df_position,df ]
result = pd.concat(frames, axis=1)

result


# ## 7) File to csv

# In[38]:


#reference code: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html

compression_opts = dict(method='zip',
                        archive_name='subset2.csv')  

result.to_csv('subset2.zip', index=False,
          compression=compression_opts)  

