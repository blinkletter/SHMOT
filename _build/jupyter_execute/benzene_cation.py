#!/usr/bin/env python
# coding: utf-8

# # Example
# 
# How is SHMOT **useful**? It is mostly a teaching tool. Calculating the coefficients and electronic properties will build a better understanding of **conjugation**, **aromaticity** and **reactivity**. Let us use the example of the **benzyl cation**.
# 
# ## The Benzyl Cation
# 
# Perhaps nothing is as **electron withdrawing** as a carbocation. The benzyl cation would be an extreme example of an electron-withdrawing group.
# 

# ### Step 1: Number the Atoms
# 
# We will **draw** the benzyl cation and **number** the atoms.
# 
# ```{figure} images/benzylcation.png
# ---
# width: 100px
# name: benzylcation
# ---
# *The benzyl cation*
# ```
# 
# Is the **positive charge** really on the atom as shown? HMOT will reveal all.
# 

# ### Step 2: Create The Interaction Matrix
# 
# I will **copy and paste** the *Python* code from the previous chapter for the benzyl cation. I could easily type it in myself but since **it already exists** I will use it.
# 

# In[1]:


import numpy

## Benzyl Cation

Name = "benzyl cation"
Connections = numpy.array([[0,1,0,0,0,0,0],
                           [1,0,1,0,0,0,1],
                           [0,1,0,1,0,0,0],
                           [0,0,1,0,1,0,0],
                           [0,0,0,1,0,1,0],
                           [0,0,0,0,1,0,1],
                           [0,1,0,0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 6


# ### Step 3: Create the HuckelSolver Object
# 
# We will use the HuckelSolver tool as before. Boom! Done!

# In[2]:


from shmo import HuckelSolver

Result = HuckelSolver(z,n)


# ### Step 4: Display  the Result
# 
# We will use my `print_shmo.py` library and the `Print_It_All()` function to access the HMOT results.

# In[3]:


from print_shmo import Print_It_All

Print_It_All(Result, Name)


# ## MO Diagram
# 
# First, lets **construct** a MO diagram using the coefficients and energies. **Observe** that the LUMO has the largest coefficients at the benzyl carbon and the *ortho* and *para* positions. This where **nucleophiles** might add. Nucleophiles can add to the benzyl cation where the largest part of the empty LUMO is located.
# 
# But we are interested in **how** an EWG affects **electrophilic aromatic substitution**. The HOMO shows the coefficients are greatest in the *meta* and *ortho* positions. We know that EWG **deactivate** the *ortho* and *para* sites and that would **leave** the *meta* carbon as the location for electrophilic substitution.
# 
# ```{figure} images/benzylcationMO.png
# ---
# width: 500px
# name: benzylcationMO
# ---
# *The benzyl cation Hückel molecular orbital diagram*
# ```
# 
# ## Bond Orders
# 
# In **benzene**, all the bonds would have a &pi; order of 0.5 and a total bond order of 1.5. I read the **data** for the benzyl cation off the tables above and made this image. You can see that there is no longer perfect symmetry but the bonds are still all **partially** double bonds. The bond to the benzyl cation is drawn as a single bond but clearly has significant double bond character due to **conjugation**.
# 
# ```{figure} images/benzylcationbondorder.png
# ---
# width: 200px
# name: benzylcationbondorder
# ---
# *The bond orders in the benzyl cation HMOT model*
# ```
# 
# ## Electron Density
# 
# I read the reported values for electron density and created this **diagram**. You can see the **greatest** electron density where an electrophilic aromatic substitution could take place is the *meta* position.
# 
# ```{figure} images/benzylcationelectrondensity.png
# ---
# width: 200px
# name: benzylcationelectrondensity
# ---
# *The electron density in the benzyl cation HMOT model*
# ```
# 
# ## Partial Charges
# 
# I read the charges from the result data and created the following **diagram**. We can see that the **largest** charge is indeed resident on the benzyl carbon. With the positive charge also shared at the ortho and para positions. An electron-deficient carbon is **unlikely** to attach an electrophile. The *meta* carbon is all that is left.
# 
# ```{figure} images/benzylcationcharges.png
# ---
# width: 200px
# name: benzylcationcharges
# ---
# *The electron density in the benzyl cation HMOT model*
# ```

# ## A Challenge
# 
# **Repeat** the calculation above but **add** two electrons to model the idea of a benzyl **anion**. That would be a very effective **electron-donating** group, would it not? What is the HOMO now? Would that HOMO **enable** electrophilic addition at the *ortho* and *para* positions as we would expect? How are **charges** and electron density affected? What are the **bond orders**? Draw new **diagrams** to collate this information and make conclusions.
# 
# Perhaps a **better model** for an **electron-donating** group would be a methyl group (toluene) or oxide (phenol anion). Perhaps a better model for an **electron-withdrawing** group would be a carbonyl group (benzaldehyde). Model those and examine the *ortho*,*meta*,*para* changes in the benzene ring.
# 
# We are **limited** in the choices of structure as we want to keep our Hückel systems **simple**. We could model nitrobenzene if we wanted to, but the polar parameter adjustments would be more complicated. The references in [chapter 3](Polar_Atoms_HMOT) will give you **details** about that if you want to go there.
# 
# ## Be Prepared
# 
# We will be discussing the electronic properties of **naphthaline** and **azulene** in class. The interaction matrices are all set up for you in [chapter 4](SHMO_Github_more). Use them and **calculate** all the properties demonstrated in this worksheet. Make some diagrams and predict the regiochemistry in both cases. Are these molecules **non-polar** – or are the surprisingly **polar**?

# ## Documentation
# 
# The following *Python* commands were demonstrated in this workbook. They arte gathered here with links to the *Python* documentation.
# 
# The documentation for the **NumPy** library can be found [here](https://numpy.org/doc/stable/).
# 
# - [numpy.array()](https://numpy.org/doc/stable/reference/generated/numpy.array.html)
# - [shmo.HuckelSolver()](https://github.com/randlet/SHMO)
# 
