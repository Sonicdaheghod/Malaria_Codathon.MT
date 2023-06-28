#!/usr/bin/env python
# coding: utf-8

# In[64]:


#make a table of data
#http://alimanfoo.github.io/2018/04/09/selecting-variants.html


# # The goal of this project with the team I am in  (Justin Gibbons, Indu Parameswaran, Shenica Jerome) is to determine a specific variant of a gene that caused patients to contract Malaria from a dataset consisting of patients with Malaria.
# 
# ![image.png](attachment:image.png)
# 

# # 1) Import and read through gzip file

# In[20]:


pip install scikit-allel


# In[1]:


pip install "setuptools<58"


# In[33]:


pip install pyvcf


# # 2) Organize Data 
# 
# ## The subset dataself is a chuck of the whole dataset, makes the code easier to run

# In[16]:


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

# In[15]:


#O(n^x): Polynomial complexity notation

for record in vcf_file:
    for sample in record.samples:
        print(sample["GT"])
        
print(record.samples)


# ## 3a) Get the GT info seperately

# In[21]:


# import scikit-allel
import allel
# check which version is installed
print(allel.__version__)


# ## 3b) Get different titles from dataset

# In[25]:


callset = allel.read_vcf(r"C:\Users\Megan Tran\Desktop\Megan's USB\College\Code\Codathon\subsampled_vcf.gz")
sorted(callset.keys())


# ## 3c) convert 0/1 marking to martix for only GT label

# In[26]:


callset['calldata/GT']


# ##  - filter subtelomericHypervariable
# ## 4) Get homozygous data [0,0], [1,1]

# In[82]:


## here, needed to attach np to array part of code

gt_zarr = callset['calldata/GT']
        
import numpy as np

arrays = gt_zarr

for array in arrays:
    if np.array_equal(array, np.array([1, 1])) or np.array_equal(array, np.array([0, 0])):
        print(array)
        
arrays


# # 5) Graph Data
