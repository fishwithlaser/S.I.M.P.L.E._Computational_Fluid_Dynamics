import matplotlib                  # enables graphing
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import buildMatrix as matrix
import numpy as np

def Contour(Ans, B , x_length=1, y_length=1,name="figure", titleTxt="",debug='n'):
    x_nodes = len(Ans[1][:])
    y_nodes = len(Ans[:][1])
    dx = float(x_length) / x_nodes
    dy = float(y_length) / y_nodes
    xpos = [0] * x_nodes                      # calculating position along x
    ypos = [0] * y_nodes                      # calculating position along y
    xpos[0] = dx / 2                          # populating the first x-node
    ypos[0] = dy / 2                          # populating the first y-node
    for i in range(x_nodes-1):
        xpos[i+1] = xpos[i] + dx
    for i in range(y_nodes-1):
        ypos[i+1] = ypos[i] + dy
    # plotting line
    x_line = [0] * x_nodes
    y_line = [0] * x_nodes
    for i in range(x_nodes):
        x_line[i] = xpos[i]
        y_line[i] = y_length-ypos[i]

    if debug == "y":
        print("building mesh")
    X, Y = np.meshgrid(xpos, ypos)
    # PLOTTING CONTOUR PLOTS
    # note: these will save in a directory called "plot"
    plt.figure()
    # Graphing gradient function with a bilinear interpolation
    im = plt.imshow(Ans, interpolation='bilinear', origin='lower',
                     extent=(0, 1, 0, 1))   # plottin interpolation
    CBI = plt.colorbar(im, shrink=0.9)      # creating colour-bar
    plt.axis([0, x_length, 0, y_length])    # plotting boundaries
    CS = plt.contour(X, Y, Ans, colors='k') # plotting contours
    plt.clabel(CS, fontsize=9, colors='k')  # populating lable
    plt.plot(x_line,y_line)                 # plotting the transcects
    plt.title(titleTxt)
    plt.savefig("plots/"+name+"Ans.png", bbox_inches='tight')
    plt.gcf().clear()
    print("figure saved")

def Error(Ans, x_length = 1, name="figure"):
    x_nodes = len(Ans[1][:])                  
    dx = x_length / x_nodes
    xpos = [0] * x_nodes                      # calculating position along x
    y_nodes = len(Ans[:][1]) 
    a_line = [0] * x_nodes
    x_line = [0] * x_nodes
    for i in range(x_nodes-1):
        xpos[i+1] = xpos[i] + dx
    for i in range(x_nodes):
        a_line[i] = Ans[y_nodes-i-1][i]
        x_line[i] = xpos[i]
    plt.figure()
    yaxisMax = (max(max(Ans)) + max(max(Ans))*0.05)
    plt.axis([0, x_length, 0, yaxisMax])
    titleTxt = ""
    ylabText = "Temperature (" + '\u00b0 C' + ")"
    plt.ylabel = ylabText
    plt.plot(x_line,a_line)
    plt.title(titleTxt)
    plt.savefig("plots/"+name+"Err.png" , bbox_inches='tight') 


