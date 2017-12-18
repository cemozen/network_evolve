#!/usr/bin/python
import networkx as nx
from networkx.generators.random_graphs import *
import matplotlib.pyplot as plt
from scipy import integrate
#import pylab as p
import numpy as np
import random
#import numpy.random
#import time
import sys
import ast
import odespy
import simplejson as json
import sqlite3
import matplotlib.cm as cm
import peakdetect as pdt

def generate_random_network(nodes, links, seed=None, directed=True):
    return gnm_random_graph(nodes, links, seed=seed, directed=directed)
    
def test_connected(G):
    return nx.is_weakly_connected(G)

def add_regulation_modes(G, num_act, num_rep):
    if num_act + num_rep == G.number_of_edges():
        for u, v, y in G.edges(data=True):
            y['regulation']='R' ; y['aweight']=0 ; y['rweight']=1 ; y['color']='black'; y['arrowhead']='tee'
        for u, v, y in np.random.permutation(G.edges(data=True))[0:num_act]:
            y['regulation']='A' ; y['aweight']=1 ; y['rweight']=0 ; y['color']='red' ; y['arrowhead']='normal'

def init_random_network(num_nodes,num_links,num_act,num_rep,seed=None):
    # returns a weakly connected network of num_act directed activator edges
    # and num_rep directed repressor edges.
    directed=True   # we demand directed networks anyway
    G = generate_random_network(num_nodes,num_links,seed,directed)
    add_regulation_modes(G,num_act,num_rep)
    connected=test_connected(G)
    #counter=0
    while not connected:
       #counter+=1
       #print "trying to make connected:", counter
#       seed +=1
       G = generate_random_network(num_nodes,num_links,seed,directed)
       add_regulation_modes(G,num_act,num_rep)
       connected=test_connected(G)
    return G

def network_info(G,regulated=False):
    print "number of nodes=", G.number_of_nodes()
    print "number of edges=", G.number_of_edges()
    if regulated:
        L=[j['regulation'] for j in [i[2] for i in G.edges(data=True)][:]]
        print "number of repressor edges=", L.count('R')
        print "number of activator edges=", L.count('A')
    print "number of self loops=", G.number_of_selfloops()
    
def graph_to_list(G):
    g_list= [(m[0],m[1],m[2]["regulation"]) for m in G.edges(data=True)]
    return g_list
    
def list_to_graph(edges):
   # creates a network where each edge's nodes and regulation type is specified by custom.
   # each edge has to be passed on using a triplet (i,j,mode) where node i regulates node j
   # by the given mode (repression or activation). The default values belong to the represillator. 
   rep_edge_data={'arrowhead': 'tee','aweight': 0,'color': 'black','regulation': 'R','rweight': 1}
   act_edge_data={'arrowhead': 'normal','aweight': 1,'color': 'red','regulation': 'A','rweight': 0}
   G=nx.DiGraph() 
   for i,j,mode in edges: 
      if mode=="R":
          G.add_edge(i,j,**rep_edge_data)
      else:
          G.add_edge(i,j,**act_edge_data)
   return G

def network_plot(G,file,format='eps',prog='dot'):
    K=nx.to_agraph(G)
    K.draw(file,format=format,prog=prog)
    
def network_plot_from_edges():
    # Given edge data as string, generates the graph and prints its network picture
    try:
       edges=ast.literal_eval(sys.argv[1])
       outfile=sys.argv[2]
#      print edges, type(edges)
    except IndexError:
       print "Usage:",sys.argv[0],"edges outfile"
       sys.exit(1)  
    G=list_to_graph(edges)
    network_plot(G,outfile,format='eps',prog='dot')    
    
def dynamix(G,phi,n,u0,tstart,tstop,dt,method="RK4"):
    A,R = extract_matrices(G)
    u_i=initialize_nodes(G,u0)
    t=np.arange(tstart,tstop,dt)
    if method == "RK4":        # the best
       solver = odespy.RK4(lambda u,t: du_dt(n,u,A,R,phi,t=None))
       solver.set_initial_condition(u_i)
       u, t = solver.solve(t)
    elif method == "ABM3":     # less accurate than RK4 but faster
       solver = odespy.AdamsBashMoulton3(lambda u,t: du_dt(n,u,A,R,phi,t=None))
       solver.set_initial_condition(u_i)
       u, t = solver.solve(t)
    elif method == "odeint":    # not a good one, fails often
        u=u_i
        u=integrate.odeint(lambda u,t: du_dt(n,u,A,R,phi,t=None),u,t)
    else:
        raise SystemExit('ODE method %s is not implemented yet!' % method)
    return u,t
    

def dynamix_plot_from_edges():
    # Given edge data as string, generates the graph and calculates its dynamics then
    # outputs the 0th node dynamics
    try:
       edges=ast.literal_eval(sys.argv[1])
       phi=float(sys.argv[2])     # unitless regulation strength parameter defined as in Kobayashi et.al.
       n=float(sys.argv[3])
#      print edges, type(edges)
    except IndexError:
       print "Usage:",sys.argv[0],"edges phi n"
       sys.exit(1)    
    G=list_to_graph(edges)
    u0=0.1    # initial value for node 0, other nodes will be assigned 0.
    tstart=0. # start time
    tstop=300. # stop time
    dt=0.01 # time interval
    u,t=dynamix(G,phi,n,u0,tstart,tstop,dt,method="RK4")
    for i in range(len(t)):
        sys.stdout.write("%8.5f" % (t[i]))
        for j in u[i]:
           sys.stdout.write("%8.5f" % (j))
        sys.stdout.write("\n")
    print dynamic_checker(4,t,u)
        
def dynamix_3dplot_from_edges():
    # Given edge data as string, generates the graph and calculates its dynamics then
    # outputs the 0th node dynamics
    try:
       edges=ast.literal_eval(sys.argv[1])
       phi=float(sys.argv[2])     # unitless regulation strength parameter defined as in Kobayashi et.al.
       n=float(sys.argv[3])
#      print edges, type(edges)
    except IndexError:
       print "Usage:",sys.argv[0],"edges phi n"
       sys.exit(1)    
    G=list_to_graph(edges)
    u0=0.1    # initial value for node 0, other nodes will be assigned 0.
    tstart=0. # start time
    tstop=300. # stop time
    dt=0.01 # time interval
    u,t=dynamix(G,phi,n,u0,tstart,tstop,dt,method="RK4")
    np.savetxt(sys.stdout,u)


def extract_matrices(G):
    # matrices defined in such a way that element ij means gene i is
    # regulated (activated or repressed) directly by gene j.
    # that's why we return the transpose matrices.
    A = np.array(nx.adjacency_matrix(G,weight='aweight'))
    R = np.array(nx.adjacency_matrix(G,weight='rweight'))
    return A.T, R.T
  
def rewire(G):
    possible_edges=[(i,j) for i in G.nodes() for j in G.nodes()]
    used_edges=G.edges()
    nonused_edges=list(set(possible_edges).difference(set(used_edges)))
    edge_to_remove=random.choice(used_edges)
    d=G.get_edge_data(*edge_to_remove)
    G.remove_edge(*edge_to_remove)
    edge_to_add=random.choice(nonused_edges)
    G.add_edge(*edge_to_add,**d)

def mutate(G):
    # connected to start with.
    rewire(G)
    connected=test_connected(G)
    while not connected:
       rewire(G)
       connected=test_connected(G)

def initialize_nodes(G,u0=0.1):
    u=np.zeros(G.number_of_nodes())
    u[0]=u0    
    # below is intended only for a particular application to Zhang et.al.
    #u=np.array([0.0,0.0,0.05])
    return u

def du_dt(n,u,A,R,phi,t=None):
    # This function determines the dynamics according to Kobayashi et.al, Eur. Phys. J. B 76, 167, (2010). 
    term1=1./(1.+(phi*np.dot(R,u))**n)
    term2=(phi*np.dot(A,u))**n/(1.+(phi*np.dot(A,u))**n)
    for x in np.nditer(term2, flags=['refs_ok'], op_flags=['readwrite']):
        if abs(x)<1e-6:
           x[...]=1.0
    term3=np.multiply(term1,term2)
    return term3-u    

def x_at_y(x1,y1,x2,y2,y=0.):
   # return x for a given value of y, satisfying y=ax+b, a and b 
   # found from given values of (x1,y1), (x2,y2).
   x=(x2-x1)*y-(x2*y1-x1*y2)
   try:
      return x/(y2-y1)
   except ValueError:
      print "Division by zero encountered in x_at_y"
      return None

def event_finder(x,y):
   # return events defined by abscissae for which y=y_ave and y'>0
   events=[]
   h=(min(y)+max(y))/2.
   eps=1.e-04   # define a tolerance to avoid artifical events, this depends on numerical accuracy
   for i in range(len(y)-1):
      if y[i+1]>h+eps/2. and y[i]<h-eps/2.:
         # to find values y=y_ave, we use linear interpolation:
         events.append(x_at_y(x[i],y[i],x[i+1],y[i+1],y=h)) 
#   print len(events)
   return events
      
def dynamic_checker(num_nodes,t,u):
   # for a more accurate calculation for period, 
   # determine events by searching only on the second half of the plot 
   # we care only for the observation node, which is taken as the 0th.
   events=event_finder(t[(len(t)/2):],u.T[0][(len(t)/2):]) 
   print "events:"
   print events
   if len(events) > 1: 
       intervals=np.diff(events)
       print "intervals:"
       print intervals
       T=np.mean(intervals) ; sig=np.std(intervals)
       return (T,sig)

def cost(T0,T,sig):
    return (T-T0)**2/T0**2+(sig/T)**2
    

def ensemble():
    try:
       num_nodes=int(sys.argv[1])
       num_links=int(sys.argv[2])
       num_act=int(sys.argv[3])
       phi=float(sys.argv[4])     # unitless regulation strength parameter defined as in Kobayashi et.al.
       n=float(sys.argv[5])             # cooperativity constant
       ensemble_size=int(sys.argv[6])
    except IndexError:
       print "Usage:",sys.argv[0],"num_nodes, num_links, num_act, phi, n, ensemble_size"
       sys.exit(1)
    
    num_rep=num_links-num_act
    tstart=0.        # start time
    tstop=60.        # stop time
    #  WARNING: Maximum period that can be caught is about (tstop-tstart)/2. 
    dt=0.01          # time interval
    u0=0.1           # initial value for the observer gene concentration
    count=0
    dyn_count=0

    t=np.arange(tstart,tstop,dt)
    G=init_random_network(num_nodes,num_links,num_act,num_rep,seed=None)
    sys.stdout.write("starting the random ensemble:\n")
    sys.stdout.flush()
    sys.stdout.write("%8s %8s %7s %7s\n" % ("count", "dynamic", "T", "sig"))
    sys.stdout.flush()
    # Create new or update existing database. set create_new=True if you intend to overwrite the existing database.
#    ensemble_update_db_table(create_new=True)
    #--------------------------------
    # ensure it is all-dynamic, i.e., no gene is expressed in a static steady state.
    for iter in range(ensemble_size):
         dynamic=None
         while not dynamic:
              count+=1
              G=init_random_network(num_nodes,num_links,num_act,num_rep,seed=None)
              A,R = extract_matrices(G)
              u_i=initialize_nodes(G,u0)
                           
              solver = odespy.RK4(lambda u,t: du_dt(n,u,A,R,phi,t=None))
              solver.set_initial_condition(u_i)
              u, t = solver.solve(t)
              dynamic=dynamic_checker(num_nodes,t,u)

         # obtain the period and its variance for the observer gene's (bt default node=0) output for
         # the initial network, then get its score:
         dyn_count+=1
         T,sig=dynamic
         sys.stdout.write("%8d %8d %7.4f %7.4f\n" % (count, dyn_count, T, sig))
         sys.stdout.write("%s %s\n" %  ("Network:", graph_to_list(G)))
         sys.stdout.flush()
         # fill in the database:
#         ensemble_fill_db_table(count, dyn_count, T, sig, G)
    sys.stdout.write("Ensemble completed!\n")
    sys.stdout.flush()



def accept_reject(eps,eps_t,mu):
    accepted=False
    diff=eps_t-eps
    if diff<0.:
       accepted=True
    else:
       p=np.exp(-diff/eps/mu)
       eta=np.random.rand()
       if eta<=p:
           accepted=True
    return accepted

def sweep(G,T,sig,eps,n,phi,T0,mu,t,u0,num_nodes,num_act,num_rep,nsweep,nmutation,naccept):
   # start with a mutated trial copy of G, keep mutating it until it is an all-dynamic network
   max_mutation=10000   # maximum number of mutations to prevent infinite loops
   nmutation=0          # number of mutations per sweep
   nsweep+=1            # counter for sweep rounds
   G_t=G.copy()
   # ensure we deal with all-dynamic networks
   dynamic=None
   while not dynamic:
      mutate(G_t)
      nmutation+=1      # counter for mutations
      A_t,R_t = extract_matrices(G_t)      
      u_t = initialize_nodes(G_t,u0)
      solver = odespy.RK4(lambda u,t: du_dt(n,u,A_t,R_t,phi,t=None))
      solver.set_initial_condition(u_t)
      u, t = solver.solve(t)
      dynamic=dynamic_checker(num_nodes,t,u)
      if nmutation > max_mutation:   # if still not found a dynamic network, it is regarded as rejected
         nmutation=-nmutation    # negative integer is a catch condition in the calling function
         return nsweep,nmutation,naccept,G,T,sig,eps

   # obtain the period and its variance for the observer gene's (bt default node=0) output for
   # the trial network, then get its score:
   T_t,sig_t=dynamic      
   eps_t=cost(T0,T_t,sig_t)
   # accept or reject the trial network:
   accepted=accept_reject(eps,eps_t,mu)
   # if the network is accepted, update the current network and its score to the trial one and its score. 
   if accepted:
       G=G_t
       T=T_t
       sig=sig_t
       eps=eps_t
       naccept+=1
   return nsweep,nmutation,naccept,G,T,sig,eps

def metropolis(num_nodes, num_links, num_act, phi, n, trial_threshold, T0):
      num_rep=num_links-num_act      
      mu=0.2
      delT=0.05
      delsig=0.001
      tstart=0.        # start time
      tstop=100.     # stop time
      dt=0.01          # time interval
      u0=0.1       # initial value for the observer gene concentration
      t=np.arange(tstart,tstop,dt)
      # start with a weakly connected network 
      G=init_random_network(num_nodes,num_links,num_act,num_rep)
      # ensure it is all-dynamic, i.e., no gene is expressed in a static steady state.
      dynamic=None
      while not dynamic:
          G = init_random_network(num_nodes,num_links,num_act,num_rep)
          A,R = extract_matrices(G)
          u_i = initialize_nodes(G,u0)
          solver = odespy.RK4(lambda u,t: du_dt(n,u,A,R,phi,t=None))
          solver.set_initial_condition(u_i)
          u, t = solver.solve(t)
          dynamic=dynamic_checker(num_nodes,t,u)
      # obtain the period and its variance for the observer gene's (bt default node=0) output for
      # the initial network, then get its score:
      T,sig=dynamic      
      eps=cost(T0,T,sig)
      convergence=abs(T-T0)       

      nsweep=0
      nmutation=0
      naccept=0
      acc_ratio=0.0
           
      successful=False      
      while nsweep < trial_threshold :
          nsweep,nmutation,naccept,G,T,sig,eps=sweep(G,T,sig,eps,n,phi,T0,mu,t,u0,num_nodes,num_act,num_rep,nsweep,nmutation,naccept)
          acc_ratio=float(naccept)/float(nsweep)       
          if nmutation < 0:  # catch condition for exceeding maximum allowed mutation number
             return (successful,)
          data=graph_to_list(G)
#          sys.stdout.write("%8d %12d %9.4f %7.4f %7.4f %11.4f %7.4f\n %s\n " % (nsweep, nmutation, acc_ratio, eps, T, abs(T-T0), sig, data)) 
#          sys.stdout.flush()                          
          if abs(T-T0)/T0<delT and sig<delsig :
             successful=True # loop terminated successfully 
             break
      if successful:
          return (successful,nsweep,nmutation,acc_ratio,eps,T,sig,data)
      else:
          return (successful,)

def metropolis_ensemble():
    try:
       num_nodes=int(sys.argv[1])
       num_links=int(sys.argv[2])
       num_act=int(sys.argv[3])
       phi=float(sys.argv[4])     # unitless regulation strength parameter defined as in Kobayashi et.al.
       n=float(sys.argv[5])             # cooperativity constant
       trial_threshold=int(sys.argv[6])
       T0=float(sys.argv[7])
       ensemble_size=int(sys.argv[8])
    except IndexError:
       print "Usage:",sys.argv[0],"num_nodes, num_links, num_act, phi, n, trial_threshold, T0, size_of_ensemble "
       sys.exit(1)
#    sys.stdout.write("starting the evolutionary game for phi:%d, n:%d, T0:%d\n\n" % (phi, n ,T0))
    sys.stdout.write("%9s  %9s  %9s  %7s %7s %7s %7s %7s\n" % ("ens_no", "nsweep", "nmutat", "acc_rat", "eps", "T","conv","sig"))
    sys.stdout.flush()
    cnt=0  # counter for succesful members of the ensemble achieving the target.
    for ensemble_no in range(1,ensemble_size+1):
        mtrp=metropolis(num_nodes, num_links, num_act, phi, n, trial_threshold, T0)
        if mtrp[0]:   # if metropolis succceeds for the given circuit 
            nsweep,nmutation,acc_ratio,eps,T,sig,data=mtrp[1:]
            sys.stdout.write("%9d  %9d  %9d  %7.4f %7.4f %7.4f %7.4f %7.4f\n" % (ensemble_no, nsweep,nmutation,acc_ratio,eps,T,abs(T-T0)/T0, sig))
            sys.stdout.write("%s %s\n" %  ("Network:", data))
            sys.stdout.flush() 
            cnt+=1 
        else:
            sys.stdout.write("Ensemble member %9d fails to yield the target period!\n"  % (ensemble_no))
            sys.stdout.flush()              
    success_rate=(float(cnt)/ensemble_size)
    sys.stdout.write("%30s  %7.4f\n"  % ("Target success rate is:", success_rate))
    sys.stdout.flush()
               
    # log plot of "target_periods vs % evolution success"
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.set_yscale('log')
#    ax.set_xscale('log')
#    ax.set_ylabel('% evolution success')
#    ax.set_xlabel('target periods')
#    ax.plot(target_periods, success_rates, color="blue", linewidth="2")
#    plt.show()



def define_db(database_file):
      global db
      db = sqlite3.connect(database_file)
      global cursor
      cursor = db.cursor()

def ensemble_update_db_table(create_new=False):
      if create_new:
           # create new by ignoring an existing table
           cursor.execute('''DROP TABLE IF EXISTS ensemble_table''')
      cursor.execute('''CREATE TABLE ensemble_table
                 (count INT, dyn_count INT, T REAL, sig REAL, data_str TEXT)''')

def ensemble_fill_db_table(count, dyn_count, T, sig, G):
      # edge data (list of tuples)           
      data=graph_to_list(G)
      # convert "data" to string type
      data_str = json.dumps(data)
      # insert many rows of data
      cursor.execute('''INSERT INTO ensemble_table VALUES(?,?,?,?,?)''',(count, dyn_count, T, sig, data_str))
      # commit the change
      db.commit()
     

    

def ensemble_query_db_table(dyn_count=None, T=None, sig=None, tol=None):
      # tol above is to be supplied by the user for giving a tolerance value for a real-valued search item
      found=False
      if dyn_count !=None:
          found=True
          cursor.execute('''SELECT data_str FROM ensemble_table WHERE dyn_count = ? ''', (dyn_count,))
          (call,) = cursor.fetchone()
          print call
          data_list = json.loads(call)
          F=list_to_graph(edges=data_list)
          plotname='graph_from_row'+str(dyn_count)+".eps"
          network_plot(F,plotname,format='eps',prog='dot')
          print "Network plot has been saved to:", plotname
          graph_dynamics(F,T0=10.0,phi=117.65,n=3.)
          
#          ensemble_graph(F,dyn_count,T0=10.0,phi=117.65,n=3.)          
          
      elif T != None:
          cursor.execute('''SELECT dyn_count, T FROM ensemble_table''')
          T_list = cursor.fetchall()
          for (db_dyn_count, db_T) in T_list:
              if abs(T-db_T) < tol:
                  found=True
                  cursor.execute('''SELECT data_str FROM ensemble_table WHERE dyn_count = ? ''', (db_dyn_count,))
                  print "db_dyn_count=", db_dyn_count
                  (call,) = cursor.fetchone()
                  print call
                  data_list = json.loads(call)
                  F=list_to_graph(edges=data_list)    
                  plotname='graph_from_row'+str(db_dyn_count)+".eps"
                  network_plot(F,plotname,format='eps',prog='dot')
                  print "Network plot has been saved to:", plotname
#                  ensemble_graph(F,dyn_count,T0=10.0,phi=117.65,n=3.)          

      elif sig != None:    
          cursor.execute('''SELECT dyn_count, sig FROM ensemble_table''')
          sig_list = cursor.fetchall()
          for (db_dyn_count, db_sig) in sig_list:
              if abs(sig-db_sig) < tol:
                  found=True
                  cursor.execute('''SELECT data_str FROM ensemble_table WHERE dyn_count = ? ''', (db_dyn_count,))
                  print "db_dyn_count=", db_dyn_count
                  (call,) = cursor.fetchone()
                  print call
                  data_list = json.loads(call)
                  F=list_to_graph(edges=data_list)    
                  plotname='graph_from_row'+str(db_dyn_count)+".eps"
                  network_plot(F,plotname,format='eps',prog='dot')
                  print "Network plot has been saved to:", plotname
#                  ensemble_graph(F,dyn_count,T0=10.0,phi=117.65,n=3.)          

      if not found:
          print "no match has been found!"
      db.close()



def metropolis_update_db_table(create_new=False):
      if create_new:
           # create new by ignoring an existing table
           cursor.execute('''DROP TABLE IF EXISTS metropolis_table''')
      cursor.execute('''CREATE TABLE metropolis_table
                 (nsweep INT, nmutation INT, eps REAL, T REAL, sig REAL, convergence REAL, data_str TEXT)''')

def metropolis_fill_db_table(nsweep, nmutation, eps, T, convergence, sig, G):
      # edge data (list of tuples)           
      data=graph_to_list(G)
      # convert "data" to string type
      data_str = json.dumps(data)
      # insert many rows of data
      cursor.execute('''INSERT INTO metropolis_table VALUES(?,?,?,?,?,?,?)''',(nsweep, nmutation, eps, T, convergence, sig, data_str))
      # commit the change
      db.commit()

def metropolis_query_db_table(nsweep=None, eps=None, T=None, sig=None, tol=None):
      # tol above is to be supplied by the user for giving a tolerance value for a real-valued search item
      found=False
      if nsweep !=None:
          found=True
          cursor.execute('''SELECT data_str FROM metropolis_table WHERE nsweep = ? ''', (nsweep,))
          (call,) = cursor.fetchone()
          print call
          data_list = json.loads(call)
          F=list_to_graph(edges=data_list)
          plotname='graph_from_row'+str(nsweep)+".eps"
          network_plot(F,plotname,format='eps',prog='dot')
          print "Network plot has been saved to:", plotname
          
      elif eps != None:  
          cursor.execute('''SELECT nsweep, eps FROM metropolis_table''')
          eps_list = cursor.fetchall() # returns a tuple list
          for (db_nsweep, db_eps) in eps_list:
              if abs(eps-db_eps) < tol:
                  found=True
                  cursor.execute('''SELECT data_str FROM metropolis_table WHERE nsweep = ? ''', (db_nsweep,))
                  print "db_nsweep=", db_nsweep
                  (call,) = cursor.fetchone()
                  print call
                  data_list = json.loads(call)
                  F=list_to_graph(edges=data_list)    
                  plotname='graph_from_row'+str(db_nsweep)+".eps"
                  network_plot(F,plotname,format='eps',prog='dot')
                  print "Network plot has been saved to:", plotname

      elif T != None:
          cursor.execute('''SELECT nsweep, T FROM metropolis_table''')
          T_list = cursor.fetchall()
          for (db_nsweep, db_T) in T_list:
              if abs(T-db_T) < tol:
                  found=True
                  cursor.execute('''SELECT data_str FROM metropolis_table WHERE nsweep = ? ''', (db_nsweep,))
                  print "db_nsweep=", db_nsweep
                  (call,) = cursor.fetchone()
                  print call
                  data_list = json.loads(call)
                  F=list_to_graph(edges=data_list)    
                  plotname='graph_from_row'+str(db_nsweep)+".eps"
                  network_plot(F,plotname,format='eps',prog='dot')
                  print "Network plot has been saved to:", plotname

      elif sig != None:    
          cursor.execute('''SELECT nsweep, sig FROM metropolis_table''')
          sig_list = cursor.fetchall()
          for (db_nsweep, db_sig) in sig_list:
              if abs(sig-db_sig) < tol:
                  found=True
                  cursor.execute('''SELECT data_str FROM metropolis_table WHERE nsweep = ? ''', (db_nsweep,))
                  print "db_nsweep=", db_nsweep
                  (call,) = cursor.fetchone()
                  print call
                  data_list = json.loads(call)
                  F=list_to_graph(edges=data_list)    
                  plotname='graph_from_row'+str(db_nsweep)+".eps"
                  network_plot(F,plotname,format='eps',prog='dot')
                  print "Network plot has been saved to:", plotname

      if not found:
          print "no match has been found!"
      db.close()
      
def histogram_from_db(database_name):
    
      db = sqlite3.connect(database_name)  
      cursor = db.cursor()
      
      period_call=[]
      cursor.execute('''SELECT T FROM ensemble_table''')
      period_call=cursor.fetchall()
      period_list = [i[0] for i in period_call]
      print period_list
          
      bin_size=1; lower=0. ; upper=40.
      num_bins = (upper-lower)/bin_size    
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.set_ylabel('# of dynamic networks')
      ax.set_xlabel('period')
      ax.hist(period_list, bins=num_bins,range=(lower,upper))
      
      cursor.execute('''SELECT count, dyn_count FROM ensemble_table 
      WHERE dyn_count IN (SELECT max(dyn_count) FROM ensemble_table) ''')
      label_call = cursor.fetchall()
      dyn_count = [(i[0]) for i in label_call]
      print dyn_count

      ax.text(0.5, 1,"%i objects are chosen from training set" % dyn_count[0], fontsize=12) # the label position is problematic!
      plt.show()

def extract_all(database_name, ensemble=None, metropolis=None):    
      db = sqlite3.connect(database_name)  
      cursor = db.cursor()      
# query for "ensemble" routine:
      if ensemble != None:
          cursor.execute('''SELECT count, dyn_count, T, sig, data_str FROM ensemble_table''')      
          print "outputs (count, dyn_count, T, sig, data_str) of the randomly generated networks: "
          for row in cursor:
              print('{0}, {1}, {2}, {3}, {4}'.format(row[0], row[1], row[2], row[3], row[4]))
# query for "metropolis" routine:
      elif metropolis != None:
          cursor.execute('''SELECT nsweep, nmutation, eps, T, convergence, sig, data_str FROM output_table''')      
          print "outputs (nsweep, nmutation, eps, T, abs(T-T0), sig, data_str) of the evolutionary game: "
          for row in cursor:
              print('{0}, {1}, {2}, {3}, {4}, {5}, {6}'.format(row[0], row[1], row[2], row[3], row[4], row[5], row[6]))

def ensemble_graph(G,dyn_count,T0=10.0,phi=117.65,n=3.):
    A,R = extract_matrices(G)
    u_i=initialize_nodes(G,u0=0.1)
    tstart=0.        # start time
    tstop=10.0*T0     # stop time
    dt=0.01          # time interval
    t=np.arange(tstart,tstop,dt)
    solver = odespy.RK4(lambda u,t: du_dt(n,u,A,R,phi,t=None))
    solver.set_initial_condition(u_i)
    u, t = solver.solve(t)
    plt.figure()
    plt.plot(t,u.T[0])
    
    for i in range(G.number_of_nodes()): 
        filename='graph_'+ str('%04d' % dyn_count) + '.txt'
        filename='./repository/'+filename
        array_file = open(filename, "w")
        array_file.write("t:\n"+" ".join(str(elem) for elem in t)+"\n")
        for j in range(len(u.T)):
            array_file.write("\nu.T["+ str(j)+"]:\n")
            array_file.write(" ".join(str(elem) for elem in u.T[j])+"\n")
        array_file.close()



def bifurcate(edges):
    G=list_to_graph(edges)
    K_vals=np.linspace(0.1,0.115,50)
    n=3.
    u0=0.1
    tstart=0.
    tstop=350.
    dt=0.01
    cutoff=2000
    for i, K in enumerate(K_vals):
       phi=1./K
       u,t=dynamix(G,phi,n,u0,tstart,tstop,dt,method="RK4")
       x=u.T[0]
       x=x[-cutoff:] ; t=t[-cutoff:]
       mycolor=cm.spectral(i/float(len(K_vals)),1)
       _max, _min = pdt.peakdetect(x, t, 1, 0.0001) 
       try:
          tmax,umax=zip(*_max) 
          order=len(umax)
          for peak in umax:
             sys.stdout.write("%10.6f %10.6f %9d\n" % (K,peak,order))
             sys.stdout.flush() 
       except ValueError:
          pass


if __name__ == '__main__':   

#    ensemble()
#    metropolis_ensemble()
#    network_plot_from_edges()
    dynamix_plot_from_edges()
     # Zhang network:
#    bifurcate([(1,0,'R'),(0,2,'R'),(2,1,'R'),(0,1,'A'),(2,0,'A'),(2,2,'A')])    
#    dynamix_3dplot_from_edges()
     # network 13
#    bifurcate([(0, 1, 'A'), (0, 3, 'A'), (1, 0, 'A'), (2, 0, 'R'), (2, 1, 'R'), (2, 3, 'R'), (3, 1, 'R'), (3, 2, 'A')])
#     # network 18:
#     bifurcate([(0, 1, 'A'), (0, 2, 'A'), (0, 3, 'R'), (1, 0, 'R'), (1, 2, 'R'), (2, 0, 'A'), (2, 3, 'A'), (3, 1, 'R')])
#*********** OUTDATED STUFF ****************************
#    define_db("ensemble_5node_10link_0act.db")
#    metropolis(T0=20)
#    metropolis_ensemble(12.0,10)
#    ensemble_query_db_table(dyn_count=366)
#    graph_dynamics(G,T0=10.0,phi=117.65,n=3.)
#    metropolis_plot(period_range=10,num_of_trials=20)     
#    ensemble_query_db_table(T=10.5,tol=0.001)
#    metropolis_query_db_table(eps=1.12,tol=0.01)
#    extract_all("genera_10act.db",ensemble)
#    histogram_from_db("genera_10act.db")
