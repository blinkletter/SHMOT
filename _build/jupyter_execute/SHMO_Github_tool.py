#!/usr/bin/env python
# coding: utf-8

# # HückelSolver
# 
# Never reinvent the wheel. I am willing to bet that every *Python* program that you plan to make has **already** been written. Randal Taylor has **written** a *Python* class library that provides tools for **calculating** parameters of SHMO systems. The code can be obtained [here](https://github.com/randlet/SHMO)
# 
# A *class* is a **data object** that includes **functions** within it that can operate on the data that it contains. Imagine a **container** that can hold a number but also includes functions for incrementing and decrementing the number by a set value. This could all be wrapped into a **class** rather than a variable and two separate functions.
# 
# When you call the class function, you **create** an object with that class (it has all the properties of the class). This is the heart of "**object oriented programming**". We don't care about that, we're just going to use it.
# 

# ## Using a HuckelSolver Object
# 
# To start, ensure that the `shmo.py` file is located in the **same directory** as this notebook.
# 
# We **create** a HuckelSolver **object** by calling the `HuckelSolver()` class function. We provide the Hückel interaction **matrix** and the **number** of $\pi$ electrons. The command `z = HuckelSolver(Array, Number)` will create the object `z` that contains all the Huckel **parameters** for the system defined by the interaction matrix `Array` with `Number` $\pi$ electrons. `Array` may be any "array-like" object. It can be a list of lists, a *NumPy* array or a *NumPy* matrix. 
# 
# Once the HuckelSolver object is created from the data, the following **properties** may be queried.
# 
# - `data`
#     - This is the **connection matrix** that was input into the HuckelSolver class function that created the object being queried. 
# - `.num_electrons`
#     - An integer for the **number** of electrons that were used in the creation of the HuckelSolver object.
# - `.energies`
#     - a *NumPy* array of the **eigenvalues**. These are the relative orbital energies.
# - `.eigen_vectors`
#     - a list of *NumPy* arrays. Each item in the list is a *NumPy* array with orbital **coefficients** for the corresonding energy in `.energies`. For example, `.eigen_vectors[2]` is the array of coefficients for the orbital with the energy in `.energies[2]`. Each **array** is ordered according to atoms. So `.eigen_vectors[2][0]` is the coefficient on the first (0<sup>th</sup>) atom in the array for the orbital in position 2 of the list (which is the 3<sup>rd</sup> orbital, remember all *Python* lists and arrays begin at position "0").
# - `.charge_densities`
#     - A *NumPy* array of the $\pi$ electron density on each atom.
# - `.net_charges`
#     - A *NumPy* array of the net charge (due to $\pi$ electrons) on each atom.
# - `.bond_orders`
#     - a *NumPy* matrix with bond orders. Look to the cells that **correspond** to bonds to get the order of the $\pi$ bond between the two referenced atoms (e.g., cell 2,3 will contain the bond order between atoms 2 and atom 3). Cells that do not reference bonds contain values but are not important.
# - `.energy_eigens`
#     - This is a **combination** of `.energies` and `.energies`. It is a list of lists for each orbital. Each sublist contains the orbital energy and the corresponding array of eigenvectors.
# - `.populated_levels`
#     - This is a subclass that contains **three types of data**. It is a list of data for each orbital. Each orbital is represented by an object that contains the energy (of the given molecular orbital), eigenvector (orbital coefficients at each atom) and the number of electrons in the orbital.
# 
# 

# ### Step 1: Import Your Tools
# There is a *Python** package for just about every endeavour. We will use the *NumPy* package to provide **array objects** and we will use the *shmo* package that provides the **HuckelSolver class function**. *NumPy* is included in most *Python* distributions. `shmo.py` can be downloaded from Randal Taylor's github site at <https://github.com/randlet/SHMO>

# In[1]:


from shmo import HuckelSolver
import numpy


# ### Step 2: Create The Interaction Matrix
# Using the **rules** of Simplified Hückel Molecular Orbital Theory (SHMO) we can create a 2-dimensional *NumPy* array that represents the **bond connectivity** for our $\pi$ system. We can also provide **adjustments** for heteroatoms as described in the [previous chapter](Polar_Atoms_HMOT).
# 
# Below are **two options**. Execute a single choice of the blocks below and then proceed to step 3. Note that in these examples, I am using position zero as the first atom. This **matches** the way *Python* addresses arrays.
# 
# ```{figure} images/acrolein_butadiene.png
# ---
# width: 300px
# name: images/acrolein_butadiene
# ---
# *Butadiene and Acrolein*
# ```
# 

# In[2]:


## Acrolein

Name = "Acrolein"
Connections = numpy.array([[1,  1,   0,  0],
                           [1,  0.1, 1,  0],
                           [0,  1,   0,  1],
                           [0,  0,   1,  0]])



z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 4


# In[3]:


## Butadiene

Name = "Butadiene"
Connections = numpy.array([[0,1,0,0],
                           [1,0,1,0],
                           [0,1,0,1],
                           [0,0,1,0]])

z = -Connections    # The HuckelSolver uses negative numbers for some reason.
n = 4


# ### Step 3: Create the HuckelSolver Object
# By inputing the **interaction matrix** and the **number of electrons**, the HuckelSolver() function will create an object with all the Hückel calculation results contained within. In the example below, we create a HuckelSolver **object** called `Result` using the interaction matrix, `z`, and the number of $\pi$ electrons, `n`. We now have all the **results**. We are done. One line of code performed all the hard work for us. And the tool was written by **someone else** - sweet!

# In[4]:


Result = HuckelSolver(z,n)


# ### Step 4: Display and Manipulate the Result
# We now have all the information we need to **interpret** a Hückel system. Let us use some code to **print out** the results. I created a **library** that contains a single function called `Print_It_All()`. If we **pass** the HuckelSolver object to this function, it will **print** all the results that we are interested in. The file `print_shmo.py` contains the code needed. Make sure it is in the **same directory** as this notebook. 

# In[ ]:





# In[5]:


from print_shmo import Print_It_All

Print_It_All(Result, Name)


# ## Summary
# We essentially wrote a **two-line program**. The **first** line used the interaction matrix and number of electrons to create the HuckelSolver object with all the solutions contained inside it. Then the **second** line printed out that information. Below is an example of the code that was executed above using the butadiene example.
# 
# ```python
# import shmo, numpy, print_shmo
# 
# Name = "Butadiene"
# Connections = numpy.array([[0,1,0,0], [1,0,1,0], [0,1,0,1], [0,0,1,0]])
# z = -Connections
# n = 4
# 
# Result = shmo.HuckelSolver(z,n)          # solve the matrix
# print_shmo.Print_It_All(Result, Name)    # print out the result 
# ```
# We set up the data and then **called** the two **tools** that were written for us in advance. *Python* is very simple when you let other people do all the work. **Take advantage** of tools that are available and give credit where it is due.

# ## A Challenge
# 
# **Download** the Jupyter notebook for this chapter and **edit** it to remove all the stuff you dont want (my boring words, for example). Create a **lean and mean** SHMOT machine. You now will have a tool for all your future SHMOT calculations. Just **change** the data and **execute** the notebook.

# ## Documentation
# 
# The following *Python* commands were demonstrated in this workbook. They arte gathered here with links to the *Python* documentation.
# 
# The documentation for the **NumPy** library can be found [here](https://numpy.org/doc/stable/).
# 
# - [numpy.array()](https://numpy.org/doc/stable/reference/generated/numpy.array.html)
# - [shmo.HuckelSolver()](https://github.com/randlet/SHMO)
