

def matrixA(x_nodes,y_nodes,debug='n'):
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
        del xB, yB
    del dim

x_nodes = 5
y_nodes = 5

matrixA(x_nodes,y_nodes, debug='y')
