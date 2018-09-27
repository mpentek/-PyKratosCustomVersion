from __future__ import print_function, absolute_import, division
import scipy
import sys
sys.path.append("..")
#print(sys.path)
import matplotlib.pyplot as plt
from numpy import *
from pyKratos import *

def Steepest_descent( model_part, strategy, objFunc, modNodes, objNodes, objVar, maxIter, epsilon, delta ):

    # Initialize nodes to be modified
    numModNodes = len(modNodes) 
    numDofs = numModNodes * 2
    x0 = []    
    for node in modNodes:
        x0.extend( node.coordinates )

    # Displacement before modification
    for i in range(len(objNodes)):
        print( 'w1_initial', objNodes[i].variables[0][objVar[i]] )
    #w1_start = objNodes.variables[0][objVar]
    #print( 'w1_initial', w1_start )

    # objective function
    def objective_function( x0 ):
        
        for u in range(len(model_part.Nodes)): # setzte alle Verschiebungen auf 0 zurück
            model_part.Nodes[u].variables[0]["DISPLACEMENT_Y"] = 0
            model_part.Nodes[u].variables[0]["DISPLACEMENT_X"] = 0

        k = 0
        for node in modNodes:
            node.coordinates = [ x0[k], x0[k+1] ] # setze Koordinaten der versch Punkte neu
            k += 2

        strategy.Solve()    # löse erneut

        ## objFunc
        var = []
        for i in range(len(objNodes)):
            var.append( objNodes[i].variables[0][objVar[i]] )
    
        F = objFunc( var ) # Wert der Lösungsfunktion

        ## penaltyFunc
        

        return F

    ##---   ----    ----       SCIPY      ----   ----   ----   ----
    ##---   ----    ----    ----   ----   ----   ----   ----   ----
    import numpy as np    
    res = scipy.optimize.minimize(objective_function, x0, options={'maxiter':maxIter})

    print('')
    print(' Optimization Results Scipy')
    print('')
    print('Best x:',res.x)
    print('func:',res.fun)
    print('nfev:',res.nfev)
    print('nit:',res.nit)
    norm = np.linalg.norm(res.jac)

    print('gradientnorm:', norm)

    ##---   ----    ----    ----   ----   ----   ----   ----   ----
    ##---   ----    ----     Steepest Descent    ----   ----   ----    


    ### Initialising Variables
    print(' ')
    print('----- Steepest Descent -----')

    converge = False
    iter = 0    

    delta_vec = zeros((numDofs, numDofs))
    np.fill_diagonal(delta_vec, delta)   
    
    gamma = []
    gamma.append(0.01)

    gamma_inner = np.zeros(50)
    gamma_inner[0] = 0.01

    evaluated_design_variables = []
    evaluated_results = []
    evaluated_gradientnorm = []
    count_iterates = []

    ### ----------------------- ###

    # function evaluation
    F = objective_function( x0 )   
    df = zeros( len( x0 ))
    for i in range( len( x0 )):
        F_delta = objective_function( np.add( x0, delta_vec[i] ))
        df[i] = -(( F_delta - F )/( delta )) 
    df_norm = np.linalg.norm( df )

    cur_x = x0
    print(' ')
    print('Startwert:',cur_x,'Funktionswert:',F)
    print('Iterationen:')

    while   iter < maxIter and df_norm >= epsilon and not converge:
        i = 0
        gamma.append(0)

    ###### NEWTON - LIKE    
        
        phi1 = ( objective_function( np.add( cur_x, ( gamma[iter] + delta ) * df )) - objective_function( np.add( cur_x, ( gamma[iter] ) * df ))) / delta
        phi2 = ( -2 * objective_function( np.add( cur_x, ( gamma[iter] ) * df )) + objective_function( np.add( cur_x, ( gamma[iter] + delta ) * df )) + objective_function( np.add( cur_x, ( gamma[iter] - delta ) * df ))) / delta**2
        
        q = phi1 / phi2
        gamma_inner[i+1] = (gamma_inner[i] - float( q ))
        
        while   abs((gamma_inner[i+1] ) - gamma_inner[i] ) >= epsilon:
                i += 1
                
                phi1 = ( objective_function( np.add( cur_x, ( gamma_inner[i] + delta ) * df )) - objective_function( np.add( cur_x, ( gamma_inner[i] ) * df ))) / delta
                phi2 = ( -2*objective_function( np.add( cur_x, ( gamma_inner[i] ) * df )) + objective_function( np.add( cur_x, ( gamma_inner[i] + delta ) * df )) + objective_function( np.add( cur_x, ( gamma_inner[i] - delta ) * df ))) / delta**2
                
                q = phi1 / phi2
                
                gamma_inner[i+1] = (gamma_inner[i] - float( q ))
                
                if i>20:
                        print('Newton linesearch failed')
                        break
        
        if gamma_inner[i+1] > 0 :                
                gamma[iter] = gamma_inner[i+1]
        else:
                pass
                
        print('cur_x:',cur_x,'F:',F)
        
        prev_x = cur_x
        prev_F = F
        prev_df = df

        cur_x = np.add( cur_x, gamma[iter] * df )

        if np.linalg.norm( cur_x - prev_x )  <= epsilon:
            converge = True
            print('converge = True')            
        else:
            iter += 1

        F = objective_function(cur_x)

        for i in range( len( x0 )):
            F_delta = objective_function( np.add( cur_x, delta_vec[i] ))
            df[i] = -(( F_delta - F )/( delta ))
        
        df_norm = np.linalg.norm( df )

        evaluated_design_variables.extend( [prev_x] )
        evaluated_results.append( prev_F )
        evaluated_gradientnorm.append( df_norm )
        count_iterates.append( iter )
        
    print(' ')
    print('iterates:',iter,'  best x:',prev_x,'  func:',prev_F, 'gradientnorm:',df_norm, 'Displacement:',model_part.Nodes[1].variables[0]["DISPLACEMENT_Y"])
    print(' ')


    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(count_iterates,evaluated_results,'k')
    plt.ylabel('f(x)',size=11)
    plt.xlabel('iterations',size=11)

    ax = fig.add_subplot(1, 2, 2)
    plt.plot(count_iterates,evaluated_gradientnorm,'k-')
    plt.ylabel('||f '"'"'(x)||',size =11)
    plt.xlabel('iterations',size=11)
