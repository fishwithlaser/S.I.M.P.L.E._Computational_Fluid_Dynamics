# Thomas Kosciuch
# thomas.kosciuch@mail.utoronto.ca
# Nov 28th 2018
# updated Nov. 29th 2018

from tqdm import tqdm




def boundsA(x_nodes,y_nodes,debug='n'):
    """Package builds the a-matrix for finite-difference discretization"""
    """return  w_bound, e_bound, n_bound, s_bound """
    dim = x_nodes * y_nodes         # specified to reduce calculations
    w_bound = [0] * y_nodes         # creating 0-matricies to find
    e_bound = [0] * y_nodes         # all boundary nodes
    s_bound = [0] * x_nodes
    n_bound = [0] * x_nodes

    if debug == 'y':  # Debugger notifies progress            
        print("finding boundaries")
    for i in range (x_nodes):            # boundary nodes are populated
        s_bound[i] = i                   # used to populate matricies 
        n_bound[i] = (dim - x_nodes + i)
    for i in range (y_nodes):            
        w_bound[i] = (i * x_nodes)
        e_bound[i] = (i+1) * x_nodes -1

    if debug == 'y':   # Debugger prints boundary nodes
        xB = "n:" + str(n_bound) + "\n"+ " s:" + str(s_bound) 
        yB = "e:" + str(e_bound) + "\n"+ " w:" + str(w_bound) 
        print("boundaries: \n",  xB, "\n",  yB)
    return w_bound, e_bound, s_bound, n_bound


def CDiff(x_nodes,y_nodes,D_x, D_y, F_x, F_y, Boundary, debug="n"):
    """
    Requires:
    x_nodes = number of nodes along x 
    y_nodes = number of nodes along y
    D_x = Diffusivity along x, type = float
    D_y = Diffusivity along y, type = float
    F_x = Force along X, type = matrix (x by y)
    F_y = Force along Y, type = matrix (x by y)
    Optional:
    debug, "y" if you want debug, off by default
    """
    ## setting up stuff (simple math)
    dim = x_nodes * y_nodes         # specified to reduce calculations
    w_bound = [0] * y_nodes         # creating 0-matricies to find
    e_bound = [0] * y_nodes         # all boundary nodes
    s_bound = [0] * x_nodes
    n_bound = [0] * x_nodes

    if debug == 'y':  # Debugger notifies progress            
        print("finding boundaries")
    for i in range (x_nodes):            # boundary nodes are populated
        s_bound[i] = i                   # used to populate matricies 
        n_bound[i] = (dim - x_nodes + i)
    for i in range (y_nodes):            
        w_bound[i] = (i * x_nodes)
        e_bound[i] = (i+1) * x_nodes -1

    if debug == 'y':   # Debugger prints boundary nodes
        xB = "n:" + str(n_bound) + "\n"+ " s:" + str(s_bound) 
        yB = "e:" + str(e_bound) + "\n"+ " w:" + str(w_bound) 
        print("boundaries: \n",  xB, "\n",  yB)
 
    dim = x_nodes* y_nodes
    
    a_w = a_e = a_s = a_n      = [0] * dim             
    sp_p = su_p                = 0
    sp_W = su_W = sp_E = su_E  = [0] * dim
    sp_S = su_S = sp_N = su_N  = [0] * dim

    #Boundaries
    phi_W = Boundary[0]
    phi_E = Boundary[1]
    phi_S = Boundary[2]
    phi_N = Boundary[3]

    ## answer matricies ##
    a = [0] * dim                      # make a 
    for i in range(dim):               # 0-array
        a[i] = [0] * dim
    B                          = [0] * dim
    # POPULATING USING CENTRAL DIFFERENCE
    print("populating central difference")
    pbar = tqdm(range(dim))              # load progress bar
    for i in range(dim):
        pbar.update(1)                   # progress indicator
        if i not in  w_bound:       # populating non-boundary conditions
            a[i][i - 1]        = a[i][i - 1] - (D_x + (F_x[i - 1] / 2))
            a[i][i]            = a[i][i] + (D_x + (F_x[i - 1] / 2))
        if i not in  e_bound:
            a[i][i + 1]        = a[i][i + 1] - (D_x - (F_x[i + 1] / 2))
            a[i][i]            = a[i][i] + (D_x - (F_x[i + 1] / 2))
        if i not in  s_bound:
            a[i][i - x_nodes]  = a[i][i - x_nodes] - (D_y + (F_y[i - x_nodes] / 2))
            a[i][i]            = a[i][i] + (D_y + (F_y[i - x_nodes] / 2))
        if i not in  n_bound:
            a[i][i + x_nodes]  = a[i][i + x_nodes] - (D_y - (F_y[i + x_nodes] / 2))
            a[i][i]            = a[i][i] + (D_y - (F_y[i + x_nodes] / 2))
        if i in      w_bound:      # populating boundary conditions
            sp_W[i]            = -(2*D_x + F_x[i])
            su_W[i]            = (2*D_x + F_x[i])* phi_W
            a[i][i]            = a[i][i] - sp_W[i]
            B[i]               = B[i] + su_W[i]
        if i in      e_bound:
            sp_E[i]            = -(2*D_x - F_x[i])
            su_E[i]            = (2*D_x - F_x[i])* phi_E
            a[i][i]            = a[i][i] - sp_E[i]
            B[i]               = B[i] + su_E[i]
        if i in      s_bound:
            sp_S[i]            = -(2*D_y + F_y[i])
            su_S[i]            = (2*D_y + F_y[i])* phi_S
            a[i][i]            = a[i][i] - sp_S[i]
            B[i]               = B[i] + su_S[i]
        if i in      n_bound:
            sp_N[i]            = -(2*D_y - F_y[i])
            su_N[i]            = (2*D_y - F_y[i]) * phi_N
            a[i][i]            = a[i][i] - sp_N[i]
            B[i]               = B[i] + su_N[i]
    return a, B



def UDiff(x_nodes,y_nodes,D_x, D_y, F_x, F_y, Boundary, debug="n"):
    """
    Requires:
    x_nodes = number of nodes along x 
    y_nodes = number of nodes along y
    D_x = Diffusivity along x, type = float
    D_y = Diffusivity along y, type = float
    F_x = Force along X, type = matrix (x by y)
    F_y = Force along Y, type = matrix (x by y)
    Optional:
    debug, "y" if you want debug, off by default
    """
    ## setting up stuff (simple math)
    dim = x_nodes * y_nodes         # specified to reduce calculations
    w_bound = [0] * y_nodes         # creating 0-matricies to find
    e_bound = [0] * y_nodes         # all boundary nodes
    s_bound = [0] * x_nodes
    n_bound = [0] * x_nodes

    if debug == 'y':  # Debugger notifies progress            
        print("finding boundaries")
    for i in range (x_nodes):            # boundary nodes are populated
        s_bound[i] = i                   # used to populate matricies 
        n_bound[i] = (dim - x_nodes + i)
    for i in range (y_nodes):            
        w_bound[i] = (i * x_nodes)
        e_bound[i] = (i+1) * x_nodes -1

    if debug == 'y':   # Debugger prints boundary nodes
        xB = "n:" + str(n_bound) + "\n"+ " s:" + str(s_bound) 
        yB = "e:" + str(e_bound) + "\n"+ " w:" + str(w_bound) 
        print("boundaries: \n",  xB, "\n",  yB)
 
    dim = x_nodes* y_nodes
    
    a_w = a_e = a_s = a_n      = [0] * dim             
    sp_p = su_p                = 0
    sp_W = su_W = sp_E = su_E  = [0] * dim
    sp_S = su_S = sp_N = su_N  = [0] * dim

    #Boundaries
    phi_W = Boundary[0]
    phi_E = Boundary[1]
    phi_S = Boundary[2]
    phi_N = Boundary[3]

    ## answer matricies ##
    a = [0] * dim                      # make a 
    for i in range(dim):               # 0-array
        a[i] = [0] * dim
    B                          = [0] * dim
    # POPULATING USING CENTRAL DIFFERENCE
    print("populating upwind difference")
    pbar = tqdm(range(dim))
    for i in range(dim):
        pbar.update(1)
        if i not in w_bound:       # populating non-boundary nodes
            a[i][i - 1]        = a[i][i - 1] - (D_x + max((F_x[i - 1]),0))
            a[i][i]            = a[i][i] + (D_x + max((F_x[i - 1]),0))
        if i not in e_bound:
            a[i][i + 1]        = a[i][i + 1] - (D_x + max(-(F_x[i + 1]),0))
            a[i][i]            = a[i][i] + (D_x + max(-(F_x[i + 1]),0))
        if i not in s_bound:
            a[i][i - x_nodes]  = a[i][i - x_nodes] - (D_y + max((F_y[i - x_nodes]),0))
            a[i][i]            = a[i][i] + (D_y + max((F_y[i - x_nodes]),0))
        if i not in n_bound: 
            a[i][i + x_nodes]  = a[i][i + x_nodes] - (D_y + max(-(F_y[i + x_nodes]),0))
            a[i][i]            = a[i][i] + (D_y + max(-(F_y[i + x_nodes]),0)) 
        if i in w_bound:           # populating boundary nodes
            sp_W[i]            = -(2*D_x + max((F_x[i]),0))
            su_W[i]            = (2*D_x +  max((F_x[i]),0))* phi_W
            a[i][i]            = a[i][i] - sp_W[i] 
            B[i]               = B[i] + su_W[i]
        if i in e_bound:     
            sp_E[i]            = -(2*D_x + max(-(F_x[i]),0))
            su_E[i]            = (2*D_x + max(-(F_x[i]),0))* phi_E
            a[i][i]            = a[i][i] - sp_E[i]
            B[i]               = B[i] + su_E[i]
        if i in s_bound:     
            sp_S[i]            = -(2*D_y + max((F_y[i]),0))
            su_S[i]            = (2*D_y +  max((F_y[i]),0))* phi_S
            a[i][i]            = a[i][i] - sp_S[i]
            B[i]               = B[i] + su_S[i]
        if i in n_bound:     
            sp_N[i]            = -(2*D_y + max(-(F_y[i]),0))
            su_N[i]            = (2*D_x + max(-(F_y[i]),0))* phi_N
            a[i][i]            = a[i][i] - sp_N[i]
            B[i]  
    return a, B


