from sympy import *
import numpy as np 
import pandas  as pd
import pickle
import os

##############################↓↓↓↓↓↓↓↓↓↓    Model building    ↓↓↓↓↓↓↓↓###############################
class E_coli():
    def __init__(self):      #initial value setting            
        self.icd=False
        self.x0=0.02#biomass
        self.x=list()
        self.cg=list()
        self.V=list()
        self.ce=list()

        self.Ksgx=0.02#Monod glucose to biomass saturation constant
        self.Ksex=0#Monod AKG to biomass saturation constant
        self.Kie=0.5#Inhibition coefficient of AKG
        self.Keu=0#Generation coefficient of AKG
        self.km=0.02#Rate constant for endogenous metabolism
        self.Yxg=0.35#Biomass yield coefficient on glucose
        self.Yxe=3#Biomass yield coefficient on AKG
        self.Yeg=0.5#AKG yield coefficient on glucose
        self.umg=1#Maximum specific growth rate
        

    def gK(self,Ks,c):
        if c>=0:
            rs=(c/(Ks+c))
        else:
            rs=0
        return rs
    def gKi(self,Ki,c):
        if Ki==0:
            rs=1
        else:
            if c>=0:
                rs=(Ki/(Ki+c))
            else:
                rs=1
        return rs
    def gKp(self,Kp,c):
        if Kp!=0:
            if c>=0:
                rs=(c/(Kp+c))
            else:
                rs=0
        else:
            rs=1
        return rs
class Multi_c(E_coli):
    def __init__(self):
        super().__init__()
    
    
    def grow_dt(self,b1=None):
        self.ug=self.umg*self.gK(self.Ksgx,b1.cg0)*self.gKp(self.Ksex,b1.ce0)*self.gKi(self.Kie,b1.ce0)

        self.ve=self.Keu*self.x0*self.ug#AKG generation rate g/h
        if self.icd:
            try:
                self.ve_c=self.x0*self.ug/self.Yxe#AKG consuming rate g/h
            except:
                self.ve_c=0
        else:
            self.ve_c=0

        self.ux0=self.ug-self.km#specific growth rate
        self.dx=self.x0*self.ux0
        self.de=self.ve-self.ve_c#AKG net formation rate g/h
        self.dg=-self.x0*self.ug/self.Yxg-self.ve/self.Yeg#Glc net formation rate （negative）

        
        self.x0+=self.dx*b1.dt#biomass increased during dt


        if self.x0<=0:
            self.x0=0

        self.x.append(self.x0)

class Strain_s(Multi_c):#AKG dornor strain (△sucA)
    def __init__(self):
        super().__init__()
        self.Keu=1.33
        self.umg=0.3
        self.Ksex=0
        self.Yxe=0
class Strain_i(Multi_c):#AKG acceptor strain (△icd)
    def __init__(self):
        super().__init__()
        self.icd=True
        self.Keu=0
        self.umg=0.66
        self.Ksex=0.02
        self.Yxg=0.3
        self.Yxe=3

class Bio_r():#bioreactor
    def __init__(self):
       
        self.dt=0.001#The time interval
        self.t0=0
        self.cg0=20#Glc Conc.
        self.ce0=0#AKG Conc.
        self.V0=1#volume
        self.Fr=0#Initial feeding rate g/l/h (Glc)
        self.Sr=0.6#Glc Conc. in feeding medium g/ml

        self.cg=list()
        self.V=list()
        self.ce=list()
        self.tl=list()
        self.Fr1=5#Auto Feeding Rate
    
    def feed_dt(self):
        self.F=self.Fr/self.Sr#Feeding rate ml/h
        self.V1=self.V0+self.F/1000*self.dt#Volume increased caused by feeding 
        self.cg0=(self.F*self.Sr*self.dt+self.cg0*self.V0)/self.V1#Glc Conc. updated g/l
        self.ce0=(self.ce0*self.V0)/self.V1#AKG Conc. updated g/L
        self.V0=self.V1 #Volume updated

        if self.Fr==0: #Auto feeding start
            if self.cg0<0.001:
                self.Fr=self.Fr1
            else:
                pass
    def sum_dt(self):
        self.cg.append(self.cg0)
        self.V.append(self.V0)
        self.ce.append(self.ce0)
        self.tl.append(self.t0)
    
    def grow_t(self,t,fr=5,strain1=None,strain2=None):
        self.t4=np.arange(0,t,self.dt)
        
        for i in self.t4:
            self.t0+=self.dt
            self.feed_dt()
            strain1.grow_dt(b1=self)
            strain2.grow_dt(b1=self)
            self.cg0+=(strain1.dg+strain2.dg)*self.dt/self.V0
            self.ce0+=(strain1.de+strain2.de)*self.dt/self.V0            
            self.sum_dt()

    def to_np(self,strain1=None,strain2=None):
        self.npx1=np.array(strain1.x)
        self.npx2=np.array(strain2.x)

        print(np.size(self.tl))

        self.npcg=np.array(self.cg)
        self.npce=np.array(self.ce)
        self.npx1ratio=self.npx1/(self.npx1+self.npx2)
        
##############################↓↓↓↓↓↓↓↓↓↓    Different um2    ↓↓↓↓↓↓↓↓############################### 
    
def topdu(name):#Different um2
    pdce=pd.DataFrame()
    pdx1=pd.DataFrame()
    pdx2=pd.DataFrame()
    pdcg=pd.DataFrame()

    pdx1ratio=pd.DataFrame()  

    ulist=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
    
    n=0
    for i in ulist:
        b1=Bio_r()
        b1.dt=0.01
        b1.Fr1=5

        e1=Strain_s()
        e2=Strain_i()
        e1.x0=0.2
        e2.x0=0.2

        e1.umg=0.5
        e2.umg=i 
        e2.Yxg=0.5   
        e1.Yxg=e2.Yxg*0.3

        e2.Kie=e1.Kie=0
        e1.km=0.03
        e2.km=0.03

        b1.ce0=0
        b1.cg0=20
        
        b1.grow_t(t=100,strain1=e1,strain2=e2)
        b1.to_np(strain1=e1,strain2=e2)

        pdx1[i]=pd.Series(b1.npx1,index=b1.tl)
        pdx2[i]=pd.Series(b1.npx2,index=b1.tl)
        pdcg[i]=pd.Series(b1.npcg,index=b1.tl)
        pdce[i]=pd.Series(b1.npce,index=b1.tl)

        pdx1ratio[i]=pd.Series(b1.npx1ratio,index=b1.tl)

        print(str(i)+' Finished')
        n+=1

    file=[pdce,pdx1,pdx2,pdcg,pdx1ratio,ulist]
    
    
    with open(name,'wb') as usr_file:     
        pickle.dump(file,usr_file) 
    print('total: '+str(n)+' Finished')
    
    pdx1ratio.to_csv('./unidirectional/ulist/x1ratio.csv')
    pdX1X2=3*pdx1+3*pdx2
    pdX1X2.to_csv('./unidirectional/ulist/ODsum.csv')
    pdce.to_csv('./unidirectional/ulist/cAKG.csv')
    print('To csv Finished')
    


##############################↓↓↓↓↓↓↓↓↓↓    Different initial popultaion ratio    ↓↓↓↓↓↓↓↓############################### 

def topd_initial_ratio(name):
    pdce=pd.DataFrame()
    pdx1=pd.DataFrame()
    pdx2=pd.DataFrame()
    pdcg=pd.DataFrame()

    pdx1ratio=pd.DataFrame()

    

    ratio=[(0,1),(1,9),(1,1),(9,1),(1,0)]
    n=0
    for i in ratio:
        b1=Bio_r()
        b1.dt=0.001

        e1=Strain_s()
        e2=Strain_i()
        e1.x0=0.1/(i[0]+i[1])*i[0]
        e2.x0=0.1/(i[0]+i[1])*i[1]

        e2.Yxg=0.5   
        e1.Yxg=0.1666
        e1.Keu=1.33

        e2.Kie=e1.Kie=0.8
        e1.km=e2.km=0.02    

        b1.ce0=0
        b1.cg0=20
        b1.Fr1=5

        b1.grow_t(t=100,strain1=e1,strain2=e2)
        b1.to_np(strain1=e1,strain2=e2)

        pdx1[i]=pd.Series(b1.npx1,index=b1.tl)
        pdx2[i]=pd.Series(b1.npx2,index=b1.tl)
        pdcg[i]=pd.Series(b1.npcg,index=b1.tl)
        pdce[i]=pd.Series(b1.npce,index=b1.tl)
        pdx1ratio[i]=pd.Series(b1.npx1ratio,index=b1.tl)
        print(str(i)+' Finished')
        n+=1

    file=[pdce,pdx1,pdx2,pdcg,pdx1ratio,ratio]
    
    
    with open(name,'wb') as usr_file:     
        pickle.dump(file,usr_file) 
    print('total: '+str(n)+' Finished')
    
    pdx1ratio.to_csv('./unidirectional/initial population ratio/x1ratio.csv')
    pdX1X2=3*pdx1+3*pdx2
    pdX1X2.to_csv('./unidirectional/initial population ratio/ODsum.csv')
    pdce.to_csv('./unidirectional/initial population ratio/cAKG.csv')
    print('To csv Finished')
  
##############################↓↓↓↓↓↓↓↓↓↓ Choose one to run ↓↓↓↓↓↓↓↓###############################

try:
    os.makedirs('./unidirectional/')
    os.makedirs('./unidirectional/initial population ratio')
    os.makedirs('./unidirectional/ulist')
except:
    pass

topdu('unidirectional.pickle')
topd_initial_ratio('unidirectional.pickle')



