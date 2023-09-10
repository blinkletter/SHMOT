# HMO Theory

Hückel molecular orbital theory (HMOT) is a **simple method** to approximate the energy and shape of **molecular orbitals**. It can be simplified further by using it only for conjugated &pi; systems.

Let us **begin** with a series of atomic *p*-orbitals that are **conjugated** together. Our first example will be **butadiene**.

## The Interaction Matrix

All atomic orbitals will **interact** in a molecule. However some interactions are much, **much more important** than others. Orbitals on **adjacent** atoms will interact far more than orbitals on **distant** atoms. Orbitals that are **aligned** along an axis of symmetry will interact more strongly that orbitals that are **not aligned**. Orbitals of **similar** size and energy will interact better compared to when energies and sizes are **mismatched**.

Let us consider only the atomic *p*-orbitals in the **conjugated &pi; system** of butadiene. We will imagine the atoms of the molecule are in the *xy* plane with the &pi;-system orbitals aligned with the *z*-axis. Only the atomic *p{sub}`z`* orbitals will be **contributing** to the &pi; system that we are modeling. 

First, let us **assign** numbers to the atoms of butadiene.

```{figure} images/1-butadiene.png
---
width: 150px
name: fig1_butadiene
---
*Butadiene with assigned numbering*
```

Now we can construct the **interaction matrix**. We will have an "*A*" value for each interacting **pair** of atoms. *A{sub}`1,2`* describes the interaction between the *p{sub}`z`* atomic orbitals of atoms 1 and 2 (likely **strong**) while *A{sub}`1,4`* describes the iteraction between the *p{sub}`z`* orbitals of atoms 1 and 4 (likely very **weak**).

Below is the interaction matrix for butadiene.

$$
\begin{bmatrix}
    A_{11}       & A_{12} & A_{13} &  A_{14} \\
    A_{21}       & A_{22} & A_{23} &  A_{24} \\
    A_{31}       & A_{32} & A_{33} &  A_{34} \\
    A_{41}       & A_{42} & A_{43} &  A_{44} \\
\end{bmatrix}
$$

### The Definitions

We will now set some definitions to **simplify** the matrix. The first **statement** is that atoms that are **not adjacent do not interact**.  This is a severe assumption, but it makes the math much easier. It's mostly true (almost). This assumption, more than any others, is why we deem HMOT an **approximation** rather than an accurate calculation. 

So let us set *A{sub}`i,j`* to zero for all pairs that are not connected. For example, *A{sub}`1,3` = 0*. 

$$
\begin{bmatrix}
    A_{11}       & A_{12} & 0      &  0 \\
    A_{21}       & A_{22} & A_{23} &  0 \\
    0            & A_{32} & A_{33} &  A_{34} \\
    0            & 0      & A_{43} &  A_{44} \\
\end{bmatrix}
$$

All the remaining orbitals are interacting. There will be iteractions between **adjacent** orbitals and interactions of each orbital with **itself**. Both of these interactions are defined by the following equation.

$$
A_{i,j} = H_{i,j} – E\cdot S_{i,j}
$$

*H* is the "coulomb integral" and it is given a **value** of &alpha; when *i = j* and a value of &beta; when *i &ne; j*.  The value of &alpha; can be considered to be the energy of the isolate atomic *p*-orbital. The **value** of &beta; is the change in the energy due to overlap between adjacent orbitals. We will soon see that the calculated energy of the MOs will be above and below the **value** of &alpha; by some multiple of &beta;.

*S* is the "overlap integral" and we define it as having a **value of one** when *i = j* and **zero** when *i &ne; J*. This is another extreme **simplification** that keeps HMOT in the approximation category.

*E* is the **energy** of the system and this is what we are **solving** for. There will be four values of *E*, one for each molecular orbital.

So at each **individual** atom (*i = j*), we will have $A_{1,1} = \alpha – E$. When pairs of atoms are **adjacent**, we will have $A_{1,2} = \beta$. Now let us enter those values into the interaction matrix.

$$
\begin{bmatrix}
    \alpha – E       & \beta      & 0          &  0          \\
    \beta            & \alpha – E & \beta      &  0          \\
    0                & \beta      & \alpha – E &  \beta      \\
    0                & 0          & \beta      &  \alpha – E \\
\end{bmatrix}
$$

Now divide everything by &beta; and we get the following result.

$$
\begin{bmatrix}
   \frac{\alpha – E}{\beta} & 1          & 0          &  0 \\
    1                       & \frac{\alpha – E}{\beta}  & 1      &  0 \\
    0                       & 1          & \frac{\alpha – E}{\beta}  &  1 \\
    0                       & 0          & 1      &  \frac{\alpha – E}{\beta}  \\
\end{bmatrix}
$$

We can now **solve** for the values of $\frac{\alpha – E}{\beta}$. We can set $x = \frac{\alpha – E}{\beta}$ in the matrix. Each value of *x* will give us a **value** for *E* in terms of &alpha; and &beta;

Now the matrix is as follows.

$$
\begin{bmatrix}
    x        & 1      & 0      &  0 \\
    1        & x      & 1      &  0 \\
    0        & 1      & x      &  1 \\
    0        & 0      & 1      &  x \\
\end{bmatrix}
$$

### Solving the Matrix

The **eigenvalues** of the matrix will give us the energies of each molecular orbital. The eigenvalues are the **roots** of the **determinant** of the matrix. In a 4&times;4 matrix, like the one above, there will be four possible solutions for *x*. Once we have a value for *x*, we can get the **energy** as follows.

$$
\begin{align*}
x & = \frac{\alpha – E}{\beta} \\
\therefore E & = \alpha -x \cdot \beta
\end{align*}
$$

The **energies** will be &alpha;, the original starting point, raised or lowered by a **multiple** of &beta;.

This is not a linear algebra course but, if you recall your **linear algebra**, you will know that a diagonal matrix represents a set of secular equations and the eigenvalues are the solutions to each equation. Once we have the **eigenvalues**, we can solve the equations and determine the values of the coefficients at each carbon atom in each orbital. 

We can skip all that with the knowledge that the **eigenvectors** of the matrix are these sets of **coefficients**. There will be an **eigenvalue** for each eigenvector. There will be an **energy** and a set of **coefficients** for each molecular orbital. Determining the eigenvectors and eigenvalues of a 4&times;4 matrix is possible using established algorithms. Larger matrices require numerical methods and, fortunately, **there is an app for that**.

We will use the **tools** in the *NumPy* package to obtain the eigenvalues and their corresponding eigenvectors.

## Applying *Python*

Within the *NumPy* package is a **function** that can calculate the eigenvalues and eigenvectors. The tool is called `eig()` and it is within the *linalg* sublibrary of *NumPy*. You can **learn** more about it [here](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html). It is a function that will return a list of **eigenvalues** and a list of **eigenvectors**. We can call the function as follows.

```
eigenvalues, eigenvectors = numpy.linalg.eig(matrix)
```

In the **next chapter** I will present a Jupyter notebook in which we use this tool to calculate (estimate) the **molecular orbitals** in a molecule.

## Summary

We have briefly explored the background for **Hückel molecular orbital theory** and constructed an **interaction matrix** for butadiene. In the end, we learned this secret: just make a square **matrix** with x along the diagonal, the value of one where there is a bond **connection** and zero everywhere else. Then **solve** for x. 