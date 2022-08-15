#!/usr/bin/env python
# coding: utf-8

# # HMOT with Polar Atoms
# 
# We have set up the **interaction matrix** for butadiene as follows.
# 
# $$
# \begin{bmatrix}
#     x        & 1      & 0      &  0 \\
#     1        & x      & 1      &  0 \\
#     0        & 1      & x      &  1 \\
#     0        & 0      & 1      &  x \\
# \end{bmatrix}
# $$
# 
# If we included an **oxygen atom** in the structure, would we have the same energies and geometries in our orbitals? No, we would not. We can **change the matrix** to give us more accurate results.
# 
# ## Atoms Other than Carbon
# 
# The atomic *p*-orbitals were all considered to have an energy value of $\alpha$. **Electronegative atoms** will have a different energy for *p*-atomic orbitals compared to carbon. They will also have a **different overlap** with adjacent orbitals due to the different energy level and different size. The coulomb and overlap integrals will be **different**.
# 
# If we declare the original $\alpha$ of a **carbon** atomic *p*-orbital to be $\alpha_C$, the new **coulomb integral** will become $\alpha_X = \alpha_C + h_X\beta$. The **energy** will be changed by a factor of $\beta$. The **overlap integral** will be changed by a different factor, $k_X$. It becomes $\beta_X=k_X \beta$. The values of $h_X$ and $k_X$ for a given atom, X, have been determined and are available in **tables** in textbooks and websites. [^ref1][^ref2][^ref3]
# 
# [^ref1]: "A brief review and table of semiempirical parameters used in the Hueckel molecular orbital method", William P. Purcell and Judith A. Singer, *Journal of Chemical & Engineering Data*, **1967**, *12*, 235-246. DOI:10.1021/je60033a020 \[[Link](https://pubs.acs.org/doi/abs/10.1021/je60033a020)\]
# 
# [^ref2]: "A Pariser-Parr-Pople-based set of Hueckel molecular orbital parameters", F. A. Van-Catledge, *The Journal of Organic Chemistry*, **1980**, *45*, 4801-4802. DOI:10.1021/jo01311a060 \[[Link](https://pubs-acs-org.proxy.library.upei.ca/doi/abs/10.1021/jo01311a060)\]
# 
# [^ref3]: "Heteroatom parameters in Hückel Theory", Timothy Hughbanks
# obtained on June 15, 2022 at \[[Link](https://www.chem.tamu.edu/rgroup/hughbanks/courses/673/handouts/Huckel_Heteroatom_parameters.pdf)\]
# 
# [^ref4]: Ian Flemming, *Molecular Orbitals and Organic Chemical Reactions*, **2010**, Wiley. ISBN 9780470746585.
# 
# Below is a short table of **selected** $h_X$ and $k_X$ values (from Flemming, pg. 59[^ref4]). Note that in the "number of electrons" column, one electron implies the atom is participating in a **covalent** &pi; bond and two electrons implies that the atom is donating a **lone pair** into the &pi; system. Also note that the first entry is for *sp<sup>2</sup>* carbon and the values of $h_X$ and $k_X$ will not alter the matrix.
# 
#  | Atom               | # of e<sup>–</sup> | $h_X$ | $k_X$ |
#  | :---               |    :---:   | :---  | :---  |
#  |   =C               | 1          | 0     | 1.0   |
#  |   =N               | 1          | 0.5   | 1.0   |
#  |   –N:              | 2          | 1.4   | 0.9   |
#  |   =O               | 1          | 1.0   | 1.0   |
#  |   –O:              | 2          | 2.0   | 0.8   |
#  |   –CH<sub>3</sub>  | 2          | 2.0   | 0.7   |
# 
# ## Acrolein 
# 
# We will explore this idea by changing one of the carbon atoms in butadiene to an **oxygen**. We now have **acrolein**. We will now **repeat** all the steps for butadiene, but we need to change the interaction matrix to account for the **heteroatom**.
# 
# 
# ### Step 1
# 
# Step one is to **number** the atoms in your &pi; system. We have numbered the carbons of acrolein like so.
# 
# ```{figure} images/Acrolein.png
# ---
# width: 150px
# name: fig1_1A
# ---
# *Butadiene with assigned numbering*
# ```
# 

# ### Step 2
# 
# We now must create the interaction matrix and then translate it into a **data struture** that we can manipulate in *Python*. We will again be using a *NumPy* **array**. 
# 
# Let us **start** with a simple matrix that defines the connections. We have a matrix of zeros with *x* along the diagonal and values of one where there are bonds. The 4<sup>th</sup> atom is the **oxygen**.
# 
# $$
# \begin{bmatrix}
#    x        & 1      & 0      &  0 \\
#    1        & x      & 1      &  0 \\
#    0        & 1      & x      &  1 \\
#    0        & 0      & 1      &  x \\
# \end{bmatrix}
# $$
# 
# The polar atom will **affect** the position in the matrix that it **occupies** (the *x* at position 4,4) and the **adjacent** bonds. Recall that we defined $x = \frac{\alpha – E}{\beta}$. And now we must change the value of $\alpha$ by $h_X \beta$. We will have&hellip;
# 
# $$
# \begin{align*}
# x & = \frac{\alpha + h_X \beta – E}{\beta} \\
# \therefore x & = \frac{\alpha – E}{\beta} + h_X 
# \end{align*}
# $$
# 
# You can see that the **result** of the heteroatom is to **change** the value of *x*. The new value is $x + h_X$.
# 
# For the **adjacent** atoms we had a value of $\frac{\beta}{\beta}$ or one. Now the atoms adjacent to position 4,4 in the matrix (the oxygen) are **set** at $\frac{k_X\beta}{\beta}$ or $k_X$.
# 
# The **interaction matrix** will now look like this&hellip;
# 
# $$
# \begin{bmatrix}
#    x        & 1      & 0        &  0        \\
#    1        & x      & 1        &  0        \\
#    0        & 1      & x        &  k_X      \\
#    0        & 0      & k_X      &  x + h_X  \\
# \end{bmatrix}
# $$
# 
# Results are **more accurate** when we apply a 10% change to the coulomb integral adjacent to the heteroatom. So your **matrix** could look like the one below. We will use this **modification** in our calculations.
# 
# $$
# \begin{bmatrix}
#    x        & 1      & 0      &  0 \\
#    1        & x      & 1      &  0 \\
#    0        & 1      & x + \frac{h_X}{10}      &  k_X \\
#    0        & 0      & k_X      &  x + h_X  \\
# \end{bmatrix}
# $$
# 
# So lets do this.

# In[1]:


import numpy

hX = 1.0   # coulomb integral adjustment for oxygen (carbonyl)
kX = 1.0   # overlap integral adjustment for oxygen (carbonyl)

x=0      # set diagonal values to zero

m = numpy.array([[x, 1, 0,       0   ],
                 [1, x, 1,       0   ],
                 [0, 1, x+hX/10, 1*kX],
                 [0, 0, 1*kX,    x+hX]])   
print(m)


# ### Step 3
# 
# The interaction matrix for acrolein is now **complete**. Observe that the diagonal is **no longer** a slash of zeros. The adjacent locations might also have changed, but the factor they were **multiplied** with was 1.0, so they **appear** unchanged in this case.
# 
# We will now **calculate** the eigenvalues and eigenvectors. **Boom, we're done!** The rest of this workbook is just about **interpreting** the results.

# In[2]:


l,v = numpy.linalg.eig(m)
l=-l
print(l)
print()
print(v)


# The variable `l` is the **list** of eigenvalues and `v` is the **array** of eigenvectors. The sign was **flipped** on the eigenvalues so that the lowest energy orbitals will have negative values.

# ## Processing the Results
# 
# We will **sort** the energies (eigenvalues) and columns of coefficients (eigenvectors) as we did **previously**.

# In[3]:


order = numpy.argsort(l) # create an array that represents the order of the values
l = l[order]             # use the order to sort the eigenvalues
v = v[:, order]          # # use the order to sort the eigenvectors

# set print precision to 3 decimal places (will only appl to NumPy arrays)
numpy.set_printoptions(precision=3)   
# suppress scientific notation so small values read as 0 rather than 2E-29
numpy.set_printoptions(suppress=True) 

print("orbital energy values")
print(l)
print("orbital coefficents")
print(v)


# ## MO Diagram
# 
# We can use the eigenvalue results to construct an **energy diagram** for the &pi; MO system of acrolein. We can also use the coefficients for each atomic *p*-orbital to make graphical representations of the **orbitals** themselves. For example, the energy of the HOMO (the second lowest MO in this case) is at position 1 (the second position) in the two arrays. We can read the values from the output above or we can **display** it directly as shown below.

# In[4]:


MO = 2                    # The MO you want to grab

position = MO-1           # location in array (0 is the first position)

energy = l[position]
print(f"The energy for orbital #{position+1} is alpha{energy:+.3f}*beta \n")
print(f"The coefficents for orbital #{position+1} are: {v[:, position]}")


# With this information, we can **construct** the molecular orbitals as was demonstrated in the previous butadiene example. **Compare** the SHMOT results with the surfaces calculated for the molecular orbitals using much more expensive math. I think that the simple diagram derived from the eigenvectors is much more **informative**. What's your opinion?

# ```{figure} images/Acrolein_MOs.png
# ---
# width: 600px
# name: Acrolein_MOs
# ---
# *A molecular orbital diagram of acrolein*

# ## Calculated Results
# 
# Let us calculate the **electron density** and partial **charges** at each atom. 

# ### Electron Density
# 
# As you may expect, the oxygen will have a **slight negative** charge and the adjacent carbon will be **slightly positive**. Also observe how resonance is revealed in the math and the terminal carbon is **sharing the charge** with the carbonyl carbon.

# In[5]:


Filled = [2,2,0,0]                                               # locations of electrons

Coefficients_Squared = v**2                                      # Square all the orbital coefficients to get matrix of probabilities
Orbital_Contributions = Coefficients_Squared * Filled            # Create matrix of probabilites * electron occupancy
Electron_Densities = numpy.sum(Orbital_Contributions, axis = 1)  # Sum the values for each atom to total up electron density to get array of e-density at each atom
Charges = 1 - Electron_Densities                                 # Subtract the e-densities from 1 to get an array of charges
Total_Electrons = numpy.sum(Electron_Densities)                  # Sum the array of densities. The value should total to # or electrons
Total_Charge = numpy.sum(Charges)                                # Sum the array of charges to get value of the molecular charge

axis0,axis1 = numpy.shape(v)                                     # axis of matrix is # of atoms
for carbon in range(axis0):                                      # Print the results
    print(f"Atom {carbon+1}: electron density = {Electron_Densities[carbon]:.2f}, charge = {Charges[carbon]:+.2f}")
print("----------")
print(f"Total electrons = {Total_Electrons:.1f}")
print(f"Total charge = {Total_Charge:.1f}")


# 

# ## Summary
# 
# We have explored constructing an interaction matrix for a polar molecule and calculating the **molecular orbitals** using SHMOT. In the example of **acrolein** we observe that the molecule has an uneven distribution of charge, as expected.

# ## Documentation
# 
# The following *Python* commands were demonstrated in this workbook. They arte gathered here with links to the *Python* documentation.
# 
# The documentation for the **NumPy** library can be found [here](https://numpy.org/doc/stable/).
# 
# - [numpy.array()](https://numpy.org/doc/stable/reference/generated/numpy.array.html)
# - [numpy.argsort()](https://numpy.org/doc/stable/reference/generated/numpy.argsort.html)
# - [numpy.linalg.eig()](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html)
# - [numpy.shape()](https://numpy.org/doc/stable/reference/generated/numpy.shape.html)
# - [numpy.sum()](https://numpy.org/doc/stable/reference/generated/numpy.sum.html)
# - [numpy.set_printoptions()](https://numpy.org/doc/stable/reference/generated/numpy.set_printoptions.html)
