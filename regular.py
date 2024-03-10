import numpy as np
import matplotlib.pyplot as plt

epsilon = 0.00000001
left = 1.01
right = 4

def compute_gradient (r, n):
    angle = 2*np.pi/n 
    gradx , grady = 0, 0
    i, flip = 1, -1
    while 1-r**2*np.sin(angle/2*i)**2>=0 and i<= (n-1)//2:
        gradx += flip * np.sin(angle/2*i)*np.sqrt(1-r**2*np.sin(angle/2*i)**2)  
        grady += flip * np.cos(angle/2*i)*np.sqrt(1-r**2*np.sin(angle/2*i)**2)  
        flip *= -1   
        i += 1
    i, flip = -1, 1
    while 1-r**2*np.sin(angle/2*i)**2>=0 and i >= -(n-1)//2:
        gradx += flip * np.sin(angle/2*i)*np.sqrt(1-r**2*np.sin(angle/2*i)**2)
        grady += flip * np.cos(angle/2*i)*np.sqrt(1-r**2*np.sin(angle/2*i)**2)  
        flip *= -1   
        i -= 1
    return  -4*gradx

def compute_second_deriv (r,n):
    angle = 2*np.pi/n 
    der = 0
    i, flip = 1, -1
    while 1-r**2*np.sin(angle/2*i)**2>=0 and i<= (n-1)//2:
        der += flip * np.sin(angle/2*i)*np.cos(angle/2*i)**2/np.sqrt(1-r**2*np.sin(angle/2*i)**2)  
        flip *= -1   
        i += 1
    i, flip = -1, 1
    while 1-r**2*np.sin(angle/2*i)**2>=0 and i >= -(n-1)//2:
        der += flip * np.sin(angle/2*i)*np.cos(angle/2*i)**2/np.sqrt(1-r**2*np.sin(angle/2*i)**2) 
        flip *= -1   
        i -= 1
    return  r**3*der

def compute_asymetric_partial_deriv (r, n, k):
    #we are computing the partial deriv w.r.t. theta_i and theta_{i+k}
    angle = np.pi / n * (-k)        #{theta_i - theta_k}/2
    if 1-r**2*np.sin(angle)**2 < 0:
        return False
    ans = r * np.sqrt(1-r**2*np.sin(angle)**2)*np.sin (angle) + r**3 * np.cos(angle)**2 * np.sin(angle)/np.sqrt(1-r**2*np.sin(angle)**2)
    return ans if k%2 ==0 else -ans

def compute_hessian (r, n, mode = "normal"):
    n = int(n)
    hessian = np.zeros((n,n))
    val = compute_second_deriv(r, n)
    for i in range(n):
        hessian[i,i] = val
    k = 1
    while True:
        val = compute_asymetric_partial_deriv(r, n, k)
        if val == False:
            break
        for i in range(n):
            hessian[i, (i+k)%n] = val
        #print(hessian)
        k += 1
    k = -1
    while True:
        val = compute_asymetric_partial_deriv(r, n, k)
        if val == False:
            break
        for i in range(n):
            hessian[i, (i+k)%n] = - val
        k -= 1
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    #print("The eigenvalues of the Hessian matrix are: ")
    #print(eigenvalues)
    # print("The eigenvectors of the Hessian matrix are: ")
    # print(eigenvectors[0])
    if mode == "eigenvalues":
        return eigenvalues
    return hessian, eigenvalues

def eigenvalues_to_file(filename):

    with open('localminmax.txt', 'r') as file:
        # Open a file for writing
        with open(filename, 'w') as output_file:
            
            for line in file:
                
                n, r = map(float, line.split())
                result = compute_hessian(r, n, mode="eigenvalues")
                #print(f"Result for pair ({r}, {n}): {result}")
                output_file.write(f"Result for pair ({r}, {n}): {result}\n")
                if result[0] > 0:
                    print("STOP, we have found an example with all positive eigenvalues!!!!!") 
                    return


def second_derivs_to_file (filename):
    pairs =[]
    with open('localminmax.txt', 'r') as file:
        for line in file:
            n, r = map(float, line.split())
            pairs.append((r, n))

    # Open a file for writing
    with open(filename, 'w') as output_file:
      # Apply a function to each pair and write the result to the file
        for pair in pairs:
            result = compute_second_deriv(pair[0], pair[1])
            output_file.write(f"Result for pair ({pair[0]}, {pair[1]}): {result}\n")

def compute_second_derivs_R (r,n):
    ans = 0
    for i in range(n):
        for j in range(i+1,n):
            angle = np.pi / n * (j-i)
            if 1-r**2*np.sin(angle)**2 < 0:
                continue
            flip = 1 if (i-j)%2 == 0 else -1
            ans += flip * np.sin(angle)/np.sqrt(1-r**2*np.sin(angle)**2)
    return ans

def plot_function(n):
    # Generate x values
    r_values = np.linspace(left, right, 1000)

    # Compute y values using the provided function
    y_values = [compute_gradient(r, n) for r in r_values]

    # Plot the function
    plt.plot(r_values, y_values, label='Function')

    plt.xlabel('R values')
    plt.ylabel('Value of dg/dx')
    plt.title(f'Plot for n={n}')

    plt.grid(True)
    plt.legend()
    plt.show()

def find_all_zeros(n):
    zeros = []
    l = left
    l_grad = compute_gradient(l, n)
    for r in np.linspace(left, right, 1000):
        r_grad = compute_gradient(r, n)
        if r_grad*l_grad <0:
            zeros.append(find_zero(l, r, l_grad, r_grad, n))
        l = r
        l_grad = r_grad
    return zeros

def find_zero(l, r, l_grad, r_grad, n):
    while r_grad > epsilon or r_grad<-epsilon:
        m = (l+r)/2
        m_grad = compute_gradient(m, n)
        if m_grad*l_grad <=0:
            r = m
            r_grad = m_grad
        elif m_grad*r_grad <=0:
            l = m
            l_grad = m_grad
    return r

def find_local_minmax():
    #return a list of pairs (n,r) such that the regular polygon with n vertices  
    #on a circle of radius r has the gradient equal to zero in each vertex
    pairs = []
    for n in range(9, 21, 2):
        r_list = find_all_zeros(n)
        for r in r_list:
            pairs.append((n, r))
    with open('localminmax.txt', 'w') as file:
        for item in pairs:
          file.write(f"{item[0]} {item[1]}\n")
    return pairs

#plot_function(7)
#print(find_all_zeros(11))
#find_local_minmax()
#print(compute_second_deriv(1.0333264755140672, 39))
#second_derivs_to_file("second_derivs.txt")
#print(compute_gradient(1.3,9))
#print(compute_asymetric_partial_deriv(1.03333, 11, -1))
#print('\n'.join([' '.join([f'{element:.2f}' for element in row]) for row in compute_hessian(1.0632816909765221,13)[0]]))
#eigenvalues_to_file("eigenvalues.txt")
#print(compute_second_derivs_R(1.0525576609534186, 21))