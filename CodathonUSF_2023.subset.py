#!/usr/bin/env python
# coding: utf-8

# In[1]:


#make a table of data
#http://alimanfoo.github.io/2018/04/09/selecting-variants.html


# # The goal of this project with the team I am in  (Justin Gibbons, Indu Parameswaran, Shenica Jerome) is to determine a specific variant of a gene that caused patients to contract Malaria from a dataset consisting of patients with Malaria.
# 
# ![image.png](attachment:image.png)
# 

# # 1) Import and read through gzip file

# In[2]:


pip install scikit-allel


# In[3]:


pip install "setuptools<58"


# In[4]:


pip install pyvcf


# # 2) Organize Data 
# 
# ## The subset dataself is a chuck of the whole dataset, makes the code easier to run

# In[5]:


import vcf

vcf_file = vcf.Reader(filename=r"C:\Users\Megan Tran\Desktop\Megan's USB\College\Code\Codathon\subsampled_vcf.gz")

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

# In[6]:


#O(n^x): Polynomial complexity notation

for record in vcf_file:
    for sample in record.samples:
        print(sample["GT"])
        
print(record.samples)


# ## 3a) Get the GT info seperately

# In[7]:


# import scikit-allel
import allel
# check which version is installed
print(allel.__version__)


# ## 3b) Get different titles from dataset

# In[8]:


callset = allel.read_vcf(r"C:\Users\Megan Tran\Desktop\Megan's USB\College\Code\Codathon\subsampled_vcf.gz")
sorted(callset.keys())


# ## 3c) convert 0/1 marking to martix for only GT label

# In[9]:


GT = callset['calldata/GT']
GT


# ##  - filter subtelomericHypervariable
# ## -1 is reference, 1 is alternative
# 
# ## 4) Get homozygous data [0,0], [1,1]

# In[10]:


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

# In[24]:


#isolate sample name
sample = callset["samples"]

sample


# In[25]:


#isolate position chromosomes
variant = callset["variants/POS"]
variant


# In[26]:


#array into dictionary
#key:value
#key is the sample name, value is the genotype

#I define sample and GT as the callset function of respective titles

keys = sample
vals = GT

Malaria_dict = dict(zip(keys, zip(*vals)))
list(Malaria_dict.items())


# ## 5b) Encode Genotype to -1,1,or N/A

# ### (negative)-1 for Reference homozygous genotype (0,0), 1 for Alternative homozygous genotype (1,1), N/A for Heterozygous

# In[14]:


import allel
import numpy as np

# Load the VCF file
callset = allel.read_vcf("subsampled_vcf.gz")

# Access the genotype data
genotypes = callset['calldata/GT']

# Convert the genotypes array into a set
genotypes_set = set(tuple(map(tuple, gt)) for gt in genotypes)

# Print the first 10 values in the set
#should be able to print all
print(list(genotypes_set)[:10])


# In[15]:


# Iterate over the genotypes
def iteration():
    for record in genotypes:
        for sample in record:
            print(sample)

iteration()


# In[47]:


# Define the genotype-value mappings
genotype_value_mapping = {
    (0, 0): -1,
    (1, 1): 1,
    (0, 1): "N/A",
    (1, 0): "N/A",
    (-1, -1): "N/A"  
}

# Create an empty dictionary
genotype_dict = {}

# Iterate over the genotypes and update the dictionary
for i, sample_name in enumerate(keys):
    sample_genotypes = []
    for genotype in GT[i]:
        value = genotype_value_mapping.get(tuple(genotype), "N/A")
        sample_genotypes.append(value)
    genotype_dict[sample_name] = sample_genotypes

# Print the genotype dictionary

print("Sample and Genotype as Dictionary:")

for sample, genotypes in genotype_dict.items():
    print(f"{sample}, {genotypes}")


# ## 6) Organize Dictionary into Pandas Dataframe

# In[48]:


#top label = sample name
#vertical axis = labelling genotypes

import pandas as pd

df = pd.DataFrame(genotype_dict)
df


# In[76]:


#add column for chromosome positions 

position = callset["variants/POS"]
position

df_position = pd.DataFrame(position, columns=["Chromosome Position"])
df_position



# In[78]:


#Chromosome and Position 
chromosome = callset["variants/CHROM"]
chromosome
df_chromosome = pd.DataFrame(chromosome,columns=["Chromosome name"])

df_chromosome


# In[87]:


#concatenate chromosome, position, samples and the genotype value marking into single pd dataframe

frames = [df_chromosome,df_position,df ]
result = pd.concat(frames, axis=1)

result


# In[86]:


#My attempt:
# shift column 'Name' to first position
# first_column = df.pop("Chromosome name")
  
# # insert column using insert(position,column_name,
# # first_column) function
# df.insert(0, "Chromosome name", first_column)
  
# print()
# print("Moving last column to first:")
# display(result)

#Correct code:

print(df_chromosome.columns)

