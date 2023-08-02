import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

import scipy.sparse
import scipy.sparse.linalg


def model_channels(W):
    # define structure
    d = {'Xb': [1, 1, 1, 1, 0, 0, 1, 1, 2, 0, 0, 1, 1, 2, 1], 
     'Yb': [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4], 
     'Xe': [1, 0, 1, 2, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 1], 
     'Ye': [1, 2, 2, 2, 2, 3, 2, 3, 3, 3, 4, 3, 4, 4, 5],
     'Width': [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
     'Name': ['OA','AB','AC','AD','BC','BE','CD','CF','DG','EF','EH','FG','FH','GH', 'HI']}
    df = pd.DataFrame(data=d)

    #Support geometry functions
    def section_len(df):
        return ((df.Xe-df.Xb)**2+(df.Ye-df.Yb)**2)**0.5

    def section_angle(df):
        return np.arctan((df.Ye-df.Yb)/(df.Xe-df.Xb))
    
    # Segmentation
    dx = 0.1
    L = section_len(df)
    alpha = section_angle(df)
    X_len = df.Xe-df.Xb 
    Y_len = df.Ye-df.Yb
    ids = 0

    DFN_Mat2 = np.array([[],[],[],[],[],[],[]]).T
    for i in range(0,len(df)):
        n_seg = np.floor(L[i]/dx)
        d_seg = L/n_seg

        for j in range(0,int(n_seg)):
            newline00 = int(ids)
            newline0 = df.Xb[i] + j*X_len[i]/n_seg
            newline1 = df.Yb[i] + j*Y_len[i]/n_seg
            newline2 = df.Xb[i] + (j+1)*X_len[i]/n_seg
            newline3 = df.Yb[i] + (j+1)*Y_len[i]/n_seg
            newline4 = df.Width[i]
            newline5 = df.Name[i]
            
            DFN_Mat2 = np.vstack([DFN_Mat2,[newline00,newline0,newline1,newline2,newline3,newline4,newline5]])
            ids = ids+1
            

    sf = pd.DataFrame(DFN_Mat2, columns = ['id','Xb','Yb','Xe','Ye','Width','Name'])
    sf.iloc[:,1:6] = sf.iloc[:,1:6].astype('float')

    sf.Xb[87] =  sf.Xb[87]+0.05*W

    #Connectivity matrix
    eps = 0.001
    Connect_Mat = np.ones([len(sf),6])*(-999) 


    for i in range(0,len(sf)):
            D55 = abs(sf.Xb[:]-sf.Xb[i])
            D66 = abs(sf.Yb[:]-sf.Yb[i])
            D77 = abs(sf.Xe[:]-sf.Xe[i])
            D88 = abs(sf.Ye[:]-sf.Ye[i])
            D57 = abs(sf.Xe[:]-sf.Xb[i])
            D68 = abs(sf.Ye[:]-sf.Yb[i])
            D75 = abs(sf.Xb[:]-sf.Xe[i])
            D86 = abs(sf.Yb[:]-sf.Ye[i])
            
            index = np.where(((D55<=eps) & (D66<=eps)) | ((D77<=eps) & (D88<=eps)) | ((D57<= eps) & (D68<=eps)) | ((D75<=eps) & (D86<=eps)))
            k = 0
            for j in range(0,min(5,len(index[0]))):
                if (index[0][j]!=i):
                    Connect_Mat[i,k] = index[0][j]
                    k = k+1


    #Physics part
    #Boundary conditions (pressure for now)
    BC_coord = np.array(((-1,1,0),(1,1,5)))

    P_Inj = 1     # Injection pressure [Pa]
    P_Prod = 0    # Production pressure [Pa]
    P_ini = 0    # Initial pressure in the fractures [Pa]


    # +++++++++++++++ these need to be adjusted:
    Dt_P = 100      # Pressure time step for steady state [sec]
    DP_time = 10000
    n_time_P = int(DP_time/Dt_P) # Maximum number of pressure time steps

    num_Frac = len(sf)

    AP = np.zeros([num_Frac,num_Frac])  # Init pressure conductivity Matrix(LHS)
    BP = np.zeros(num_Frac)# Init pressure storativity Vector (RHS)



    #Solve for steady-state pressure
    tp = -1
    #timeP = np.linspace(Dt_P,DP_time,n_time_P)

    #P_new = [0]*num_Frac
    BP = np.zeros(num_Frac)
    APl = np.copy(AP)   #Local variable

    #Hydraulic parameters
    Gamma_P = np.zeros(num_Frac)
    Gamma_P_tot = np.zeros(num_Frac)
    Etha_P = np.zeros(num_Frac)

    Ls = section_len(sf)
    #A = 1 #aperture originally
    V = 1 #viscosity originally
    C = 1 # fluid compressibility originally


    for i in range(0,num_Frac):
            Gamma_P[i] = Dt_P*sf.Width[i]**3/(6*V*Ls[i])   # Transmissivity 1
            Etha_P[i] = sf.Width[i]*C*Ls[i]                 # Storativity 1


    for i in range(0,num_Frac):
        Gamma_P_tot[i] = Gamma_P[i]
        for j in range(0,6):
            cn = int(Connect_Mat[i,j])
            if (cn >= 0):
                Gamma_P_tot[i] = Gamma_P_tot[i] + Gamma_P[cn]   # if and g1g2 based on index
                AP[i,cn] = Gamma_P[i]*Gamma_P[cn]
                AP[i,i] = AP[i,i]-Gamma_P[i]*Gamma_P[cn]
        AP[i,i] = AP[i,i]/Gamma_P_tot[i] - Etha_P[i]
        for j in range(0,6):
            cn = int(Connect_Mat[i,j])
            if (cn >= 0):
                AP[i,cn] = AP[i,cn]/Gamma_P_tot[i]

    #SOLVER
     # Solve Finite Difference pressure
    P_old = np.ones(num_Frac)
    P_new = np.zeros([n_time_P,num_Frac])    

    tp = -1
    while (tp < n_time_P-1):
        tp = tp + 1

        #Apply initial conditions
        if (tp==0):
            P_old = P_old*P_ini     # Initial pressure
        else:
            P_old[0:num_Frac] = P_new[tp-1,:]

    #            for i in range(0,num_Frac):
    #                BP[i] = -Etha_P[i]*P_old[i]
        BP = -Etha_P*P_old

        # Apply boundary conditions
        APl[0,:] = 0                      #ASK REZA - does not do a lot
        APl[0,0] = 1
        BP[0] = P_Inj
        
        APl[len(sf)-1,:] = 0
        APl[len(sf)-1,len(sf)-1] = 1
        BP[len(sf)-1] = P_Prod
        
        # Injection points
        #APl[int(Inj_id[idinj]),:] = 0
        #APl[int(Inj_id[idinj]),int(Inj_id[idinj])] = 1
        #BP[int(Inj_id[idinj])] = P_Inj
        
        #Production points
    #        for i in range(0,np.shape(Prod_id)[0]):
    #            APl[int(Prod_id[i]),:] = 0
    #            APl[int(Prod_id[i]),int(Prod_id[i])] = 1
    #            BP[int(Prod_id[i])] = P_Prod

        # Solver for pressure
        #AP2 =  #SPARSE MATRIX NEEDED for performance
        AP2 = scipy.sparse.csc_matrix(APl)

        P_new[tp,:] = scipy.linalg.solve(AP,BP.T) #standard solution
        #P_new[tp,:] = scipy.sparse.linalg.spsolve(AP2,BP.T,permc_spec='MMD_ATA') #Fastest but can be singular
        #P_new[tp,:] = scipy.linalg.lstsq(AP,BP.T)[0] #numerical approx
        
    tss = tp

    P_new = P_new[0:tss+1,:]

    #Estimate velocity
    num_up = np.zeros(num_Frac)
    V = np.zeros(num_Frac)
    for i in range(0,num_Frac):
        P_Max = P_new[tss,i]
        num_up[i] = i
        for j in range(0,6):
            cn = int(Connect_Mat[i,j])
            if (cn >= 0 and P_new[tp,cn] > P_Max):
                P_Max = P_new[tp,cn]
                num_up[i] = cn

            V[i] = 1**2/(12*1*2)*(P_Max-P_new[tp,i])


    return sf, P_new, V


def plot_stuff(sf,P_new,V):
    fig = plt.figure()#figsize = [12,6])
    plt.subplot(2,2,1)
    #plt.gca().set_aspect('equal')
    plt.title('Channel structure')
    for i in range(0,len(sf)):
        plt.plot([sf.Xb[i],sf.Xe[i]],[sf.Yb[i],sf.Ye[i]],lw=sf.Width[i]*2)
        #plt.text(sf.Xb[i],sf.Yb[i],str(sf.id[i]))
        
        
    plt.subplot(2,2,2)
    #plt.gca().set_aspect('equal')
    plt.title('P convergence')
    plt.plot(P_new[:,23])

    plt.subplot(2,2,3)
    #plt.gca().set_aspect('equal')
    plt.title('Head distribution')
    time = np.shape(P_new)[0]-1
    P_norm = (P_new[time,:]-min(P_new[time,:]))/(max(P_new[time,:])-min(P_new[time,:]))
    cmap = matplotlib.cm.get_cmap('jet')

    for i in range(0,len(sf)):
        plot = plt.plot([sf.Xb[i],sf.Xe[i]],[sf.Yb[i],sf.Ye[i]],color=cmap(P_norm[i]))
        
        
    plt.subplot(2,2,4)
    #plt.gca().set_aspect('equal')
    plt.title('Discharge')
    V2 = V*sf.Width*1
    V2_norm = (V2-min(V2))/(max(V2)-min(V2))
    cmap = matplotlib.cm.get_cmap('jet')

    for i in range(0,len(sf)):
        plt.plot([sf.Xb[i],sf.Xe[i]],[sf.Yb[i],sf.Ye[i]],color=cmap(V2_norm[i]))
    plt.tight_layout()
    return fig


sf, P_new, V = model_channels(1)
plot_stuff(sf, P_new, V)