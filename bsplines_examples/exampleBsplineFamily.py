
__author__  = "Margaux BRULIARD"
__date__ ="19.06.2018"
__purpose__ = "Example of a B-Splines family"



###################### MODULES #################################
from bsplines import Bspline
import numpy as np
from matplotlib import pyplot as plt



##################### user functions ###########################

def splinesPlot(T, p, xmin, xmax, nx=100, show_fig=2):
    """
    		Plots all BSplines of degree p with knots vector T
    """
    

    x = np.linspace(xmin,xmax,nx)

    # number of the BSplines in the Schoenberg space
    N = len(T) - p - 1

    # create BSplines family with the knots sequence T 
    # of degree p
    bsp = Bspline(T,p)

    y = np.zeros((N,nx), dtype=np.double)
    for i in range(0,N):
        # evaluation of the i^th B-spline over x
        y[i]=bsp(x, i=i)
        plt.plot(x,y[i], label='$N_{}$'.format(i+1))
    plt.legend(loc=9, ncol=4)
    
    if show_fig==1:
        plt.savefig('rapport/img/bsplineFamily_degree2.png', transparent = True)
        print("The figure has been saved in 'rapport/img/bsplineFamily_degree2.png'")
    else:
        plt.show()    

    
    
#==============================================================================
# MAIN FUNCTION
#==============================================================================
if __name__ == "__main__":
    x_min = -2.
    x_max = 2.

    T = [0,0,0, 0.25, 0.5, 0.75, 1, 1, 1]
    # spline degree [here quadratic]
    p = 2
    print ("\nProgram Execution .... 'exampleBsplineFamily'")
    choice_fig = input ("Do you want to save the graph (1) or do you prefer to plot it (2)  ?")
    choice_fig = int(choice_fig)

    splinesPlot(T, p, xmin=0., xmax=3., nx=100, show_fig=choice_fig)

    print("\n")





