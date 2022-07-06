#!/usr/bin/env python
# coding: utf-8

# # HMOT with *Python*
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
# Below are instructions for obtaining the eigenvalues and eigenvectors using the `numpy.linalg.eig()` function.
# 
# ## Applying the Method
# 
# All of the *Python* commands in this document are already entered and **ready to go** in the Jupyter notebook that accompanies this tutorial. You can **download** the notebook for this page using the link at the top right corner of the page.
# 
# ### Step 1
# 
# Step one is to **number** the atoms in your &pi; system. We have numbered the carbons of butadiene **like so**&hellip;
# 
# ```{figure} images/1-butadiene.png
# ---
# width: 150px
# name: fig1_1
# ---
# *Butadiene with assigned numbering*
# ```
# 

# ### Step 2
# 
# We now must translate the interaction matrix into a **data struture** that we can manipulate in *Python*. We will be using a *NumPy* **array**. 
# 
# Recall that we has set $x = \frac{\alpha – E}{\beta}$ and that will still be true here. However, the `eig()` function that we will be using can only work with numbers. The letter "*x*" will cause errors. We can avoid trouble by setting a **numeric** value for x. In practice, this will be to set the value of *x* to zero. The result of this statement is that we are effectively setting the start value of &alpha; to **zero**. We are saying that each atomic *p*-orbital starts at the same **relative energy** and that will be true if they are all carbon atoms.
# 
# For the individual atoms we state that&hellip;
# 
# $$
# \begin{align*}
# 0 & = \frac{\alpha – E}{\beta} \\
# \therefore E & = \alpha 
# \end{align*}
# $$
# 
# As a **quick rule**, we could say that we construct an interaction matrix by **starting** with a matrix of all zeros and inserting the value of one for all **connections** between carbon atoms. It's that easy. We now have the following **diagonal matrix**&hellip;
# 
# $$
# \begin{bmatrix}
#    0        & 1      & 0      &  0 \\
#    1        & 0      & 1      &  0 \\
#    0        & 1      & 0      &  1 \\
#    0        & 0      & 1      &  0 \\
# \end{bmatrix}
# $$
# 
#  The eigenvalues of this matrix will correspond to the **relative energy** values above and below &alpha;. That energy will be the eigenvalue multiplied by &beta;. 
# 
# The code below will **buid** the matrix as a *NumPy* array.

# In[1]:


import numpy

x=0      # set diagonal values to zero

m = numpy.array([[x,1,0,0],
                 [1,x,1,0],
                 [0,1,x,1],
                 [0,0,1,x]])   
print(m)


# ### Step 3
# 
# We will now calculate the eigenvalues and eigenvectors. Boom, **we're done**! The rest of this workbook is just about **interpreting** the results.
# 
# 
# 
# 

# In[2]:


l,v = numpy.linalg.eig(m)
l=-l
print(l)
print()
print(v)


# ```{note}
# In the code above you will see that I **flipped** the signs of the eigenvalues. This is because the energy was **defined** earlier as $E = \alpha -x \cdot \beta$. Take note of the negative sign for the relative energy term, $-x \cdot \beta$. A positive value is actually the **lowest energy**. I felt that it is easier to **interpret** if I apply that negative sign to the eigenvalues. This is why using a notebook like this to **document your procedure** is important. Now you know where the negative sign came from.
# ```
# 
# ## Processing the Results
# 
# The `eig()` function **found** the eigenvalues and their corresponding eigenvectors. It did not care about the **order** though. The code below will **sort** the eigenvalues in order from **lowest** energy to **highest** and then **sort** the eigenvectors **to match**.

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


# ### Commentary
# 
# You will see that each eigenvalue **corresponds** to the column below it, which contains the eigenvector. The **first** eigenvalue (`l[0]`) corresponds to the **first** column of the eigenvector matrix (`v[:,0]`). Observe that the lowest energy molecular orbital has **coefficients** that are all the same sign. There is **no node**, as expected.
# 
# ```{note}
# *Python* **starts** all arrays at position "0". So atoms #1 in butadiene is described by **row \#0** in the matrix. Keep this in mind. Should we have numbered the diagram starting at zero? You can, just **docummnet your choices** so the reader knows what you are doing.
# ```

# ## MO Diagram
# 
# We can use the eigenvalue results to **construct** an energy diagram for the &pi; MO system of butadiene. We can also use the coefficients for each atomic *p*-orbital to make **graphical representations** of the orbitals themselves. For example, the energy of the HOMO (the second lowest MO in this case) is at position 1 (the second position) in the two result arrays. We can **read** the values from the output above or we can display it directly as shown below.

# In[4]:


MO = 2                    # The MO you want to grab

position = MO-1           # location in array (0 is the first position)

energy = l[position]
print(f"The energy for orbital #{position+1} is alpha{energy:+.3f}*beta \n")
print(f"The coefficents for orbital #{position+1} are: {v[:, position]}")


# With this information, we can construct the following **molecular orbital**. Observe that I scaled the *p*-orbitals to show a size that is **proportional** to their contribution. The actual molecular orbital will have a similar shape when calculated and visualized in much higher quality software. **Compare** it to the estimate given by SHMOT.

# 
# ```{figure} images/Butadiene_MO_2.png
# ---
# width: 350px
# name: Butadiene_MO_2
# ---
# *A molecular orbital of butadiene*
# ```
# 

# ### Rinse & Repeat
# 
# Now that you have **made one** molecular orbital, you can **make them all**. Just change the MO choice (MO = 1, MO = 2, etc&hellip;) in the code block above and use the energy and the coefficients to build your **diagram**. Here is my own version. You see that we **start** with four atomic *p*-orbitals at an energy value of &alpha;. When they are **combined**, they form four molecular &pi; orbitals. The **relative energies** are factors of &beta; and the geometries of the orbitals are given by the **coefficients** for each atomic orbital within a &pi; orbital.
# 
# ```{figure} images/Butadiene_MO_all.png
# ---
# width: 450px
# name: Butadiene_MO_all
# ---
# *The complete &pi; molecular orbital diagram of butadiene*
# ```

# ## Calculated Results
# 
# We can make a visual diagram as above by **inspecting** the energies and the coefficients (the eigenvalues and eigenvectors). We can also get quantitative **information** from the eigenvector matrix about electron density, bond orders, free valence values and more. Its just **math**.
# 
# The first important information we need to state is orbital **occupancy**. Butadiene has four electrons in the &pi; system. The **first two** orbitals are occupied and the final two are not, as seen in the diagram above. We can **express** this as an array of values.

# In[5]:


Filled = [2,2,0,0]


# We now know where the **electrons** are and we have calculated the geometry of each **molecular orbital**. This will allow use to calculate **properties** that arise from electrons and molecular orbitals together.

# ### Electron Density
# For each atom, the electron density is the **square** of the atomic *p*-orbital **coefficient** multiplied by **number** of electrons in each orbital. 
# 
# Thus, the *&pi;*-electron **charge density** (a misnomer), *q<sub>r</sub>* at atom *r* is given by
# 
# $$
# q_r = \sum_i n_i c_{i,r}^2
# $$
# 
# where *n<sub>i</sub>* is the **number** of electrons in the *i<sup>th</sup>* molecular orbital and *c<sub>i,r</sub>* is the **coefficient** of atom *r* in the *i<sup>th</sup>* molecular orbital.
# 
# 
# For butadiene, the electron density will be **even across** all atoms (in the Hückel approiximation) and so all atoms should have the **same charge** of zero. 

# In[6]:


Coefficients_Squared = v**2                                      # Square all the orbital coefficients to get matrix of probabilities
Orbital_Contributions = Coefficients_Squared * Filled            # Create matrix of probabilites * electron occupancy
Electron_Densities = numpy.sum(Orbital_Contributions, axis = 1)  # Sum the values for each atom to total up electron density to get array of e-density at each atom
Charges = 1 - Electron_Densities                                 # Subtract the e-densities from 1 to get an array of charges
Total_Electrons = numpy.sum(Electron_Densities)                  # Sum the array of densities. The value should total to # or electrons
Total_Charge = numpy.sum(Charges)                                # Sum the array of charges to get value of the molecular charge

axis0,axis1 = numpy.shape(v)                                     # axis of matrix is # of atoms
for carbon in range(axis0):                                      # Print the results
    print("Atom {}: electron density = {:.2f}, charge = {:+.2f}"
          .format(carbon, Electron_Densities[carbon], Charges[carbon]))
print("----------")
print("Total electrons = {:.1f}".format(Total_Electrons))
print("Total charge = {:.1f}".format(Total_Charge))


# ### More properties
# 
# There are more **electronic** properties that can be obtained from the orbital coefficients. 
# 
# #### Bond Orders
# 
# Bond order is derived from the **sum** of the **products** of adjacent **coeficients** and number of **electrons**. The bond order between atoms *r* and *s* is *p<sub>r,s</sub>*. It can be expressed as&hellip;
# 
# $$
# p_{r,s} = \sum_i n_i c_{i,r} c_{i,s}
# $$
# 
# where *n<sub>i</sub>* is the number of electrons in the *i<sup>th</sup>* molecular orbital and *c<sub>i,r</sub>* & *c<sub>i,r</sub>* are the coefficients of adjecent atoms in the bond in that *i<sup>th</sup>* molecular orbital. **Sum it up** and you get the &pi; bond order. If *p<sub>r,s</sub>* = 1 then we have a perfect double bond.  
# 
# Applying a square function to the whole matrix was easy in the case of electron density, but now we would need to **write** a small amount of code to **iterate** through the coefficients and perform this calculation. I'm not going to do that because we will soon be using a package of code that does all this for us. It was written by **someone else** and I will not be reinventing that wheel.
# 
# #### Delocalization Energy
# 
# Delocalization energy, *E<sub>&pi;</sub>*, is the difference between the total electronic energy in our system and the total energy if we **separated** groups of atoms. For example, the electronic energy of butadiene is the sum of the products of the energy level of the *i<sup>th</sup>* molecular orbital, *&epsilon;<sub>i</sub>*, and the number of electrons in that orbital, *n<sub>i</sub>*.
# 
# $$
# E_{\pi} = \sum_i n_i \epsilon_{i}
# $$
# 
# In butadiene, there are two electrons in the first MO (*&pi;<sub>1</sub>*) and two in the second (*&pi;<sub>2</sub>*). The **electronic energy** (in the &pi; system) is therefore&hellip;
# 
# $$
# \begin{align*}
# E_{\pi} & = n_1 \epsilon_{1}+n_2 \epsilon_{2} \\
# E_{\pi} & = 2(\alpha - 1.618\beta) +  2(\alpha - 0.618\beta) \\
# E_{\pi} & = 4\alpha - 4.472\beta) \\
# \end{align*}
# $$
# 
# Now we need to **compare it** to a system where there has been a change in the &pi; system. One idea is to imagine that **one atom** is converted to *sp<sup>3</sup>*. We would then need to do a **new** Hückel analysis of the new system and **compare** the new delocalization to the former.
# 
# 
# #### Free Valence
# The maximum valence according to HMOT for a carbon atom is 4.732. **Subtract** bond orders for each carbon atom.
# Note $F_r = 3 + \sqrt{3}$ for tertiary carbon (three connections), $F_r = 3 + \sqrt{2}$ for secondary and $F_r = 3 + \sqrt{1}$ (or $4.0$) for primary.  In all our SHMO systems we will see each carbon connected to three other atoms so we will use $F_r = 4.732$. Heteratoms (and the corresponding changes in polarity) will change these numbers but we will ignore that. 
# 
# So if we **add up** all the &pi; bond orders associated with a given atom and then add the value of 1.0 for each &sigma; bond we will get the **valence** at an atom. Subtract that value from 4.732 to obtain the **free valence** value for that atom. Larger free valence values indicate a more reactive site in the molecule (due to electronic structure. Physical structure and sterics will also play an important role.)
# 
# #### Freeloading
# 
# Free valence, delocalization energy and bond orders are simple to calculate but will require a small amount of **coding**. We will need to iterate through the eigenvectors and eigenvalues. I would do this except that there are many tools **already made** that we can use. Later we will switch to using the HückelSolver tool.
# 
# ## Summary
# 
# We have explored constructing an **interaction matrix** as a *NumPy* array and using the eig() unction to obtain the eigenvalues and eigenvectors. We have **demonstrated** this with butadiene and used the results to calculate the **electron density** and **partial charge** at each atom.
# 

# 
