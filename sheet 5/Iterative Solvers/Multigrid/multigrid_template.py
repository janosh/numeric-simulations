import numpy as np
import matplotlib.pyplot as plt

# This function improves the solution Ax=b with one Gauss-Seifel step, with the
# result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
def gaussseidel_step(x, b, N):
   #TODO: fill in your code
   #x[i,j] = ...
   
# This function calculates the resdiuum vector res = b - Ax, for input vectors
# of length N. The output is stored in res.
def calc_residuum(x, b, N, res):
    #TODO: fill in your code
    #res[i,j] = ...

# This function calculates the norm of the vector of length N,
# defined as the usual quadratic vector norm.
def norm_of_residual(res, N):
    sum = np.sum(res**2)
    return np.sqrt(sum)


#This function restricts the NxN mesh stored in 'fine[ ]' to the NNxNN mesh stored in 'coarse[ ]'
def do_restrict(N, fine, NN, coarse):
    #TODO: fill in your code
    #coarse[i,j] = ... 

# This function interpolates the the NNxNN mesh stored in 'coarse[ ]' to NxN mesh stored in 'fine[ ]'
def do_prolong(NN, coarse, N, fine):
    #TODO: fill in your code
    #fine[i,j] = ...
    

# This function carries out a V-cycle using the multigrid approach, down to N=4.
def do_v_cycle(x, b, N):
    gaussseidel_step(x, b, N)
        
    if N > 4:
        NN = N / 2;
            
        #allocate some storage for the residual and error on the current and a coarsened mesh */
        res = np.zeros((N,N))
        err = np.zeros((N,N))
        res_coarse = np.zeros((NN,NN))
        err_coarse = np.zeros((NN, NN))
        
        #calculate the residuum
        calc_residuum(x, b, N, res)
        
        #restrict the residuum
        do_restrict(N, res, NN, res_coarse)
            
        #now multiply the residuum on the coarser mesh (our b)
        #with a factor of 4 to take care of the (h) -> (2h) change, and the fact that we do not change A
        res_coarse[:,:] *= 4
            
        #call the V-cycle again, recursively
        do_v_cycle(err_coarse, res_coarse, NN)
            
        #now interpolate the error to the finer mesh
        do_prolong(NN, err_coarse, N, err)
            
        #finally add the error to the current solution
        x[:,:] += err
    
    gaussseidel_step(x, b, N)

N = 256
steps = 2000
L = 1.0
h = L / N
eta = 0.1 * L
rho0 = 10.0
  
res = np.zeros((N,N))

#now set-up the density field
x = (np.arange(-N/2,N/2)+0.5)*h
mesh = np.meshgrid(x,x)
r2 = mesh[0] * mesh[0] + mesh[1] * mesh[1]

rho = rho0* np.exp(-r2/(2*eta**2))
sum = np.sum(rho)

#initialize the starting values for phi[] and b[]
rho -= sum/(N**2)
b = 4 * np.pi*h**2*rho
phi = np.zeros((N,N))

#open a file for outputting the residuum values, and then do 2000 Jacobi steps
residuals = np.zeros(steps)
f = file("res_multigrid.txt", "w")
f.write("%d\n"%steps)

for i in np.arange(steps):
    do_v_cycle(phi, b, N)
    
    calc_residuum(phi, b, N, res)
    
    r = norm_of_residual(res, N)
    
    print("iter=%d:  residual=%g\n"%(i, r))
    f.write("%g\n"%r)
    residuals[i] = r

f.close()

plt.semilogy(residuals, label="multigrid")
plt.show()
