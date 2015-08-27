#####################
#  utility script to extract the number of degrees of freedom from a TF log file assuming SNES_VIEW has been called
#####################
import numpy as np
def getNdofs():
    """
    extract all matrix sizes from log file assuming SNES_VIEW has been called
    """
    N=[]
    filename = "terraferma.log-0"
    logfile = open(filename,"r")
    for line in logfile:
        if "cols=" in line:
            cols = line.split()[1].replace('cols=','')
            N.append(int(cols))
            
    logfile.close()
    Ndofs = np.unique(np.array(N))  
    return Ndofs

