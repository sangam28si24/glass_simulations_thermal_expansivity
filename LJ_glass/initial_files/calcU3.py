import numpy as np

def main():
    cache = np.load("alpha_cache/alpha_cache_g0.npz")
    evecs = cache["valid_evecs"]
    inv_evals = 1.0 / cache["valid_evals"]
    
    N = evecs.shape[0] // 3
    bonds = np.loadtxt("bonds_g0.txt")
    U3_Hinv = np.zeros((N, 3))
    
    for row in bonds:
        i, j = int(row[0]), int(row[1])
        rij = row[2:5]
        D2, D3 = row[6], row[7]
        psi_i = evecs[i*3:i*3+3, :]  # Shape: (3, 11997)
        psi_j = evecs[j*3:j*3+3, :]  # Shape: (3, 11997)
        d_psi = psi_j - psi_i
        C = (d_psi * inv_evals) @ d_psi.T
        f_ij = -(D3 * np.dot(rij, C @ rij) * rij + 
                 D2 * np.trace(C) * rij + 
                 2 * D2 * (C @ rij))
        
        U3_Hinv[i] += f_ij
        U3_Hinv[j] -= f_ij

    output_file = "U3_Hinv_g0.txt"
    np.savetxt(output_file, U3_Hinv, fmt="%.12e", header="U''' : H^-1 vector (Nx3)")  

if __name__ == "__main__":
    main()