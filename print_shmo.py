import numpy

def Print_It_All(Result, name = "No Name"):

    numpy.set_printoptions(precision=3)   # set print precision to 3 decimal places (only applies to numpy arrays and matices)
    numpy.set_printoptions(suppress=True) # suppress scientific notation

    print("Results for", name, "\n")
    
    print("Energies and Eigenvalues")
    print("from 'Result.energies' and 'numpy.array(Result.eigen_vectors).T)'")
    print("-"*40)
    print(Result.energies)
    print("-"*40)
    print(numpy.array(Result.eigen_vectors).T) # Transforms the list of arrays into a 2-dimensional array with columns matching each orbital energy.
    print()

    print("Energies and Probabilities (Magnitudes)")
    print("from 'Result.energies' and 'numpy.array(Result.eigen_vectors).T)**2'")
    print("-"*40)
    print(Result.energies)
    print("-"*40)
    print((numpy.array(Result.eigen_vectors).T)**2) 
    print()

    print("Electron Densities")
    print("from 'Result.charge_densities'")
    print("-"*40)
    print(Result.charge_densities)
    print()
    print("Net Atomic Charges")
    print("from 'Result.net_charges'")
    print("-"*40)
    print(Result.net_charges)
    print()
    print("Net Molecular Charge: {:.1f}".format(numpy.sum(Result.net_charges)))
    print()

    print("Bond Orders")
    print("from 'Result.bond_orders'")
    print("-"*40)
    print(Result.bond_orders)
