#!/usr/bin/env python
# coding: utf-8

# # More Molecules
# 
# Below are some **interaction matrices** that have been set up for the student. **Examine** how they are constructed, especially for the **cyclic molecules**. Cut and paste these into the HückelSolver Jupyter notebook and **explore** the electronic properties of these &pi; systems.

# 
# ```{figure} images/molecules.png
# ---
# width: 550px
# name: molecules
# ---
# *Various $\pi$ systems*
# ```
# 

# In[1]:


import numpy


# In[2]:


## Benzene

Name = "Benzene"
Connections = numpy.array([[0,1,0,0,0,1],
                           [1,0,1,0,0,0],
                           [0,1,0,1,0,0],
                           [0,0,1,0,1,0],
                           [0,0,0,1,0,1],
                           [1,0,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 6


# In[3]:


## Furan

Name = "Furan"
Connections = numpy.array([[2,   0.8, 0,   0,   0.8],
                           [0.8, 0.2, 1,   0,   0  ],
                           [0,   1,   0,   1,   0  ],
                           [0,   0,   1,   0,   1  ],
                           [0.8, 0,   0,   1,   0.2]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 6


# In[4]:


## Cyclopentadiene

Name = "Cyclopentadiene"
Connections = numpy.array([[0,1,0,0,1],
                           [1,0,1,0,0],
                           [0,1,0,1,0],
                           [0,0,1,0,1],
                           [1,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 6


# In[5]:


## Cumyl Cation

Name = "Cumyl cation"
Connections = numpy.array([[0,1,0,0,0,0,0],
                           [1,0,1,0,0,0,1],
                           [0,1,0,1,0,0,0],
                           [0,0,1,0,1,0,0],
                           [0,0,0,1,0,1,0],
                           [0,0,0,0,1,0,1],
                           [0,1,0,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 6


# In[6]:


## Naphthalein

Name = "Naphthalein"
Connections = numpy.array([[0,1,0,0,0,0,0,0,1,0],
                           [1,0,1,0,0,0,0,0,0,0],
                           [0,1,0,1,0,0,0,0,0,0],
                           [0,0,1,0,1,0,0,0,0,1],
                           [0,0,0,0,0,1,0,0,0,1],
                           [0,0,0,0,1,0,1,0,0,0],
                           [0,0,0,0,0,1,0,1,0,0],
                           [0,0,0,0,0,0,1,0,1,0],
                           [1,0,0,0,0,0,0,1,0,1],
                           [0,0,0,1,1,0,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 10


# In[7]:


## Azulene

Name = "Azulene"
Connections = numpy.array([[0,1,0,0,0,0,0,0,1,0],
                           [1,0,1,0,0,0,0,0,0,0],
                           [0,1,0,1,0,0,0,0,0,0],
                           [0,0,1,0,1,0,0,0,0,0],
                           [0,0,0,1,0,0,0,0,0,1],
                           [0,0,0,0,0,0,1,0,0,1],
                           [0,0,0,0,0,1,0,1,0,0],
                           [0,0,0,0,0,0,1,0,1,0],
                           [1,0,0,0,0,0,0,1,0,1],
                           [0,0,0,0,1,1,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 10

