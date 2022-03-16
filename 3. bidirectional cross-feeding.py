from sympy import *
import numpy as np 
import pandas  as pd
import pickle
import os

##############################↓↓↓↓↓↓↓↓↓↓    Model building    ↓↓↓↓↓↓↓↓###############################

class E_coli():
    def __init__(self):
        self.x0=0.02#biomass DCW g/L
        self.cg0=20#Glc Conc.
        self.V0=1#Volume
        self.x=list()
        self.cg=list()
        self.V=list()
        self.Ksgx=2#Monod glucose to biomass saturation constant
        self.Ksh1x=0.0#Monod H1 to biomass saturation constant
        self.Ksh2x=0.0#Monod H2 to biomass saturation constant
        self.Kh1u=0#Generation coefficient of H1
        self.Kh2u=0#Generation coefficient of H2
        self.km=0.02#Rate constant for endogenous metabolism
        self.Yxg=0.35#Biomass yield coefficient on glucose
        self.Yxh1=0#Biomass yield coefficient on H1
        self.Yxh2=0#Biomass yield coefficient on H2
        self.Yh1g=0.1#H1 yield coefficient on glucose
        self.Yh2g=0.1#H2 yield coefficient on glucose
        self.umg=1#Maximum specific growth rate
       

    def gK(self,Ks,c):
        if c>=0:
            rs=(c/(Ks+c))
        else:
            rs=0
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
        self.ug=self.umg*self.gK(self.Ksgx,b1.cg0)*self.gKp(self.Ksh1x,b1.ch1_0)*self.gKp(self.Ksh2x,b1.ch2_0)
        self.vh1=self.Kh1u*self.x0*self.ug#H1 generation rate g/h
        self.vh2=self.Kh2u*self.x0*self.ug#H2 generation rate g/h
        try:
            self.vh1_c=self.x0*self.ug/self.Yxh1#H1 consuming rate g/h
        except:
            self.vh1_c=0

        try:
            self.vh2_c=self.x0*self.ug/self.Yxh2#H2 consuming rate g/h
        except:
            self.vh2_c=0


        self.ux0=self.ug-self.km#specific growth rate
        self.dx=self.x0*self.ux0
        self.dh1=self.vh1-self.vh1_c#H1 net formation rate g/h
        self.dh2=self.vh2-self.vh2_c#H2 net formation rate g/h
        self.dg=-self.x0*self.ug/self.Yxg-self.vh1/self.Yh1g-self.vh2/self.Yh2g   #Glc net formation rate （negative）
        self.x0+=self.dx*b1.dt#biomass increased during dt
        b1.ch1_0+=self.dh1*b1.dt/b1.V0
        b1.ch2_0+=self.dh2*b1.dt/b1.V0

        if self.x0<=0:
            self.x0=0

        self.x.append(self.x0)
 
class Strain_h1(Multi_c): #H1 dornor strain
    def __init__(self):
        super().__init__()
        self.Kh1u=0.09
        self.Ksh2x=0.01
        self.Yxh2=1

class Strain_h2(Multi_c): #H2 dornor strain
    def __init__(self):
        super().__init__()
        self.Kh2u=1.0
        self.Ksh1x=0.01
        self.Yxh1=10

class Bio_r():#bioreactor
    def __init__(self):
                    

        self.dt=0.01
        self.t0=0
        self.cg0=20#initial Glc Conc.
        self.ch1_0=0
        self.ch2_0=0

        self.V0=1
        self.Fr=0#initial feeding rate g/l/h
        self.Fr1=5#Auto feeding rate g/l/h
        self.Sr=0.6#Glc Conc. in feeding medium g/ml

        self.cg=list()
        self.V=list()
        self.ch1=list()
        self.ch2=list()
        self.tl=list()
        self.Frl=list()
    
    def feed_dt(self):
        self.F=self.Fr/self.Sr#Feeding rate ml/h
        self.V1=self.V0+self.F/1000*self.dt#Volume update
        self.cg0=(self.F*self.Sr*self.dt+self.cg0*self.V0)/self.V1
        self.ch1_0=(self.ch1_0*self.V0)/self.V1
        self.ch2_0=(self.ch2_0*self.V0)/self.V1
        self.V0=self.V1 
        
        if self.Fr==0: #Auto feeding
            if self.cg0<0.001:
                self.Fr=self.Fr1
            else:
                pass
        

    def sum_dt(self):

        self.cg.append(self.cg0)
        self.V.append(self.V0)
        self.ch1.append(self.ch1_0)
        self.ch2.append(self.ch2_0)
        self.tl.append(self.t0)
        self.Frl.append(self.Fr)
    
    def grow_t(self,t,strain1=None,strain2=None):
        self.t4=np.arange(0,t,self.dt)
        
        for i in self.t4:
            self.t0+=self.dt
            self.feed_dt()
            strain1.grow_dt(b1=self)
            strain2.grow_dt(b1=self)
            self.cg0+=(strain1.dg+strain2.dg)*self.dt/self.V0
            self.sum_dt()
    def to_np(self,strain1=None,strain2=None):
        self.npx1=np.array(strain1.x)
        self.npx2=np.array(strain2.x)
        self.npcg=np.array(self.cg)
        self.npch1=np.array(self.ch1)
        self.npch2=np.array(self.ch2)
        self.npx1ratio=self.npx1/(self.npx1+self.npx2)
        print(np.size(self.tl))

##############################↓↓↓↓↓↓↓↓↓↓    Different um2    ↓↓↓↓↓↓↓↓###############################        

def topdu(name):#Different um2
    pdch1=pd.DataFrame()
    pdch2=pd.DataFrame()
    pdx1=pd.DataFrame()
    pdx2=pd.DataFrame()
    pdcg=pd.DataFrame()
    pdx1ratio=pd.DataFrame()
    
    ulist=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
    print(ulist)
    n=0
    for i in ulist:
        b1=Bio_r()
        b1.Fr1=5
        b1.dt=0.01

        e1=Strain_h1()
        e2=Strain_h2()
        e1.x0=e2.x0=0.2#inital biomass

        
        e2.Yxg=e1.Yxg=0.5
        e1.km=e2.km=0.05

        e1.Kh1u=0.5
        e1.Yh1g=1
        e2.Yxh1=3
        
        e2.Kh2u=0.5
        e2.Yh2g=1
        e1.Yxh2=3
        
        e1.Ksh2x=0.05
        e2.Ksh1x=0.05
        e1.umg=0.9
        e2.umg=i        

        b1.ch1_0=0.15
        b1.ch2_0=0.15
        b1.cg0=20
        
        b1.grow_t(t=100,strain1=e1,strain2=e2)
        b1.to_np(strain1=e1,strain2=e2)

        pdx1[i]=pd.Series(b1.npx1,index=b1.tl)
        pdx2[i]=pd.Series(b1.npx2,index=b1.tl)
        pdcg[i]=pd.Series(b1.npcg,index=b1.tl)
        pdch1[i]=pd.Series(b1.npch1,index=b1.tl)
        pdch2[i]=pd.Series(b1.npch2,index=b1.tl)
        pdx1ratio[i]=pd.Series(b1.npx1ratio,index=b1.tl)
        print(str(i)+' Finished')
        n+=1

    file=[pdx1,pdx2,pdcg,pdch1,pdch2,ulist,pdx1ratio]  
    
    with open(name,'wb') as usr_file:     
        pickle.dump(file,usr_file) 
    print('total: '+str(n)+' Finished')

topdu('bidirectional.pickle')

##############################↓↓↓↓↓↓↓↓↓↓Export the file to CSV format↓↓↓↓↓↓↓↓###############################
try:
    os.makedirs('./bidirectional/')
except:
    pass

with open('bidirectional.pickle','rb') as usr_file:    
    pd1=pickle.load(usr_file)
pdx1,pdx2,pdcg,pdch1,pdch2,ulist,pdx1ratio=pd1


pdx1ratio.to_csv('./bidirectional/x1ratio.csv')
pdX1X2=3*pdx1+3*pdx2
pdX1X2.to_csv('./bidirectional/ODsum.csv')
pdch1.to_csv('./bidirectional/cH1.csv')
pdch2.to_csv('./bidirectional/cH2.csv')
print('To csv Finished')

