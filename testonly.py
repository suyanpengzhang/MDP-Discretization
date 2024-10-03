# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 09:58:43 2022

@author: suyan
"""
import gurobipy as gp
from gurobipy import GRB
import  numpy as np
import pandas as pd
import pickle
import random
import time
import os

# Replace this with the path to your actual license file
path_to_license_file = '/Users/suyanpengzhang/gurobi.lic'

# Set the environment variable
os.environ['GRB_LICENSE_FILE'] = path_to_license_file


# Record the start time
start_time = time.time()



for limit_site in range(6,7):
    print('#######################################################################')
    try:
    
        # Create a new model
        lm = gp.Model("lm")
    
        # Create variables
        
        pi = lm.addVars(60,vtype=GRB.BINARY, name="actions") 
        x = lm.addVars(3,60,vtype=GRB.CONTINUOUS, name="states")
        # Set objective
        ##
        cost_function = gp.quicksum(x[1, t] for t in range(60))

        lm.setObjective(cost_function,GRB.MINIMIZE)

        # 
        #upper bound on x
        lm.addConstr(x[0,0]==0.9999)
        lm.addConstr(x[1,0]==0.0001)
        lm.addConstr(x[2,0]==0.0)
        for t in range(1,60):
            x[0,t]=x[0,t-1]-1.4*x[0,t-1]*x[1,t-1]
            x[1,t]=x[1,t-1]+1.4*x[0,t-1]*x[1,t-1]-0.49*x[1,t-1]
            x[2,t]=x[2,t-1]+0.49*x[1,t-1]

        # Optimize model
        #lm.setParam('TimeLimit', 10)
        lm.Params.Threads = 18
        lm.Params.OutputFlag = 1
        lm.Params.LogToConsole = 1
        lm.optimize()
        count=0
        print('Obj: %g' % lm.ObjVal)
        
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')

