'''
  File : Ex_FE_B2_01_Beam.py
  This program is written for static analysis of beam.
'''
import numpy as np
import numpy.linalg as la
import math


class Object(object):
    def __init__(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value

    def update(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value
                      
class Node(Object):
    def __init__(self,**kwargs):
        self.nid=1
        self.snid=''
        self.x=0
        self.y=0
        self.z=0
        Object.__init__(self,**kwargs)

class Element(Object):
    def __init__(self,**kwargs):
        Object.__init__(self,**kwargs)

class Beam(Element):
    def __init__(self,**kwargs):
        self.modelling='consistent'
        Element.__init__(self,**kwargs)
                
        self.mdof=4        
        self.ke=np.zeros((4,4),float)
        self.ok=np.zeros((4,4),float)
        self.om=np.zeros((4,4),float)
        
        self.ot=np.zeros((4,4),float)
        self.jdofv=np.zeros((4),int)

    def computeMatrices(self,NODE):
        
        self.c=np.zeros((self.mdof,self.mdof),float)
       
        for i in range(self.mdof):
            for j in range(self.mdof):
                self.c[i,j]=0.0

        for i in range(4):
            for j in range(4):
                self.ok[i,j]=0.0
                self.ot[i,j]=0.0
            self.ot[i,i]=1.0
           
        n1=self.nodes[0]
        n2=self.nodes[1]
        xx=NODE[n2].x-NODE[n1].x
        yy=NODE[n2].y-NODE[n1].y
        L=math.sqrt((xx*xx+yy*yy))
        self.L=L
        
        Iz=self.Iz
        E=self.E

        print(n1,n2,L,E,Iz)

        EI_by_L3 =  E * Iz / (L*L*L)
        EI_by_L2 =  E * Iz / (L*L)
        EI_by_L =  E * Iz / L

        self.ok[0,0]= 12.0*EI_by_L3
        self.ok[0,1]= 6.0*EI_by_L2
        self.ok[0,2]= -12.0*EI_by_L3
        self.ok[0,3]= 6.0*EI_by_L2

        self.ok[1,0]= 6.0*EI_by_L2
        self.ok[1,1]= 4.0*EI_by_L
        self.ok[1,2]= -6.0*EI_by_L2
        self.ok[1,3]= 2.0*EI_by_L

        self.ok[2,0]= -12.0*EI_by_L3
        self.ok[2,1]= -6.0*EI_by_L2
        self.ok[2,2]= 12.0*EI_by_L3
        self.ok[2,3]= -6.0*EI_by_L2

        self.ok[3,0]= 6.0*EI_by_L2
        self.ok[3,1]= 2.0*EI_by_L
        self.ok[3,2]= -6.0*EI_by_L2
        self.ok[3,3]= 4.0*EI_by_L
        
        if self. modelling == 'consistent':
            self.om[0,0]=156.0
            self.om[0,1]=22.0*L
            self.om[0,2]=54.0
            self.om[0,3]=-13.0*L
            self.om[1,0]=22.0*L
            self.om[1,1]=4.0*L*L
            self.om[1,2]=13.0*L
            self.om[1,3]=-3.0*L*L
            self.om[2,0]=54.0
            self.om[2,1]=13.0*L
            self.om[2,2]=156.0
            self.om[2,3]=-22.0*L
            self.om[3,0]=-13.0*L
            self.om[3,1]=-3.0*L*L
            self.om[3,2]=-22.0*L
            self.om[3,3]=4.0*L*L

            for i in range(4):
                for j in range(4):
                    self.om[i,j]=self.om[i,j]*self.m_bar*L/420.0
        else:
            self.om[0,0]=self.m_bar*L/2.0        
            self.om[1,1]=self.m_bar*L/2.0        
            self.om[2,2]=self.m_bar*L/2.0        
            self.om[3,3]=self.m_bar*L/2.0

        self.jdofv[0]=NODE[n1].dof[1]
        self.jdofv[1]=NODE[n1].dof[5]
        self.jdofv[2]=NODE[n2].dof[1]
        self.jdofv[3]=NODE[n2].dof[5]
        print(self.jdofv)
           
        del self.c

    def computeElementForce(self):
            self.femfl=[0.0,0.0,0.0,0.0]
            self.femfg=[0.0,0.0,0.0,0.0]
            L=self.L
            if 'w' in vars(self):
                w=self.w
                self.femfl[0]=self.femfl[0]+0.5*w*L
                self.femfl[1]=self.femfl[1]+w*L*L/12.0
                self.femfl[2]=self.femfl[2]+0.5*w*L
                self.femfl[3]=self.femfl[3]-w*L*L/12.0
            if 'P' in vars(self):
                P=self.P
                a=self.a
                self.femfl[0]=self.femfl[0]+P*(L-a)/L
                self.femfl[1]=self.femfl[1]+P*(L-a)*(L-a)*a/(L*L)
                self.femfl[2]=self.femfl[2]+P*a/L
                self.femfl[3]=self.femfl[3]-P*a*a*(L-a)/(L*L)
            for i in range(4):
                for j in range(4):
                    self.femfg[i]=self.femfg[i]+self.ot[i,j]*self.femfl[j]
            

    def assembleStiffness(self,gk):
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                for j in range(self.mdof):
                    jj=self.jdofv[j]
                    if jj > 0:
                        temp=gk[ii-1,jj-1]+self.ok[i,j] 
                        gk[ii-1,jj-1]=temp
                        
    def assembleLoad(self,gp):
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                temp=gp[ii-1][0]+self.femfg[i] 
                gp[ii-1][0]=temp        
                        

    def computeMemberForces(self,NODE):
            self.computeMatrices(NODE)
            mfg=[0,0,0,0]
            mfl=[0,0,0,0]
            disp=[0,0,0,0]            
            n1=self.nodes[0]
            n2=self.nodes[1]
            dispN1=NODE[n1].disp
            dispN2=NODE[n2].disp
            disp[0]=dispN1[1]
            disp[1]=dispN1[5]
            disp[2]=dispN2[1]
            disp[3]=dispN2[5]

            for i in range(4):
                mfg[i]=0.0
                for j in range(4):
                    mfg[i]=mfg[i]+self.ok[i,j]*disp[j]
                    
            for i in range(4):
                mfl[i]=0.0
                for j in range(4):
                    mfl[i]=mfl[i]+self.ot[i,j]*mfg[j]

            if 'femfl' in vars(self):
               femfl=self.femfl
                        
               for i in range(4):
                   mfl[i]=mfl[i]-femfl[i]                    

            self.mfl=mfl
            print('mfl:',mfl)
            
class Structure(Object):
    def __init__(self,**kwargs):
        self.numnode=0
        self.numelem=0
        self.DEFAULT_DOF=dict(idx=0,idy=0,idz=0,
                              irx=0,iry=0,irz=0)

        self.dof_index=dict(idx=0,idy=1,idz=2,
                            irx=3,iry=4,irz=5)

        self.force_index=dict(Fx=0,Fy=1,Fz=2,
                              Mx=3,My=4,Mz=5)
        Object.__init__(self,**kwargs)

    def setDefaultDOF(self,**kwargs):
        for key, value in kwargs.items():
            self.DEFAULT_DOF[key] = value
            
    def node(self,**kwargs):
        if 'nid' in kwargs:
            nid=kwargs['nid']
        if 'NODE_LIST' not in vars(self):
            self.NODE_LIST=list()
        if 'NODE' not in vars(self):
            self.NODE=dict()
        if nid not in self.NODE_LIST:
            self.NODE_LIST.append(nid)
            self.NODE[nid]=Node(**kwargs)
            self.numnode=self.numnode+1
        else:
            self.NODE[nid].update(**kwargs)
            
    def element(self,**kwargs):
        if 'eid' in kwargs:
            eid=kwargs['eid']
        if 'ELEM_LIST' not in vars(self):
            self.ELEM_LIST=list()
        if 'ELEM' not in vars(self):
            self.ELEM=dict()
            
        if eid not in self.ELEM_LIST:
            self.ELEM_LIST.append(eid)
            if 'etype' in kwargs:
                self.etype=kwargs['etype']
            if self.etype == 'Beam':
                self.ELEM[eid]=Beam(**kwargs)
            else:
                self.ELEM[eid]=Element(**kwargs)
            self.numelem=self.numelem+1
        else:
            self.ELEM[eid].update(**kwargs)
            
    def printResults(self):
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            idx=0
            idy=1
            irz=1
            if 'idx' in vars(N):
                idx=N.idx
            if 'idy' in vars(N):
                idy=N.idy
            if 'irz' in vars(N):
                irz=N.irz
                
            Fx=0.0
            Fy=0.0
            Mz=0.0
            if 'Fx' in vars(N):
                Fx=N.Fx
            if 'Fy' in vars(N):
                Fy=N.Fy
            if 'Mz' in vars(N):
                Mz=N.Mz
                
            print('nid:',N.nid,'snid:',N.snid,'x:',N.x,'y:',N.y,
                  'idx=',idx,'idy=',idy,'irz:',irz,'Fx=',Fx,'Fy=',Fy,'Mz=',Mz,'Dx:',N.disp[0],'Dy:',N.disp[1],'Rz:',N.disp[5])

        for eid in self.ELEM_LIST:
            E=self.ELEM[eid]
            print('eid:',E.eid,' seid:',E.seid,
                  ' Nodes:(%d,%d)'%(E.nodes[0],E.nodes[1]),
                  ' Iz :',E.Iz,'E : %.4e'%E.E,' m_bar : %.4e'%E.m_bar,
                  'MFL:',E.mfl[0],E.mfl[1],E.mfl[2],E.mfl[3])

    def decideActiveDOF(self):
        self.nndof=0
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            self.dof=[0,0,0,0,0,0]
            for jd in ['idx', 'idy', 'idz',
                       'irx', 'iry', 'irz']:
                if self.DEFAULT_DOF[jd] == 1:
                    self.nndof=self.nndof+1
                    index=self.dof_index[jd]
                    self.dof[index]=self.nndof
                self.NODE[nid].update(dof=self.dof)

    def computeSystemMatrices(self,**kwargs):
            
        self.gk=np.zeros((self.nndof,self.nndof),float)
        self.gm=np.zeros((self.nndof,self.nndof),float)
        
        self.gp=np.zeros((self.nndof,1),float)

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMatrices(self.NODE)
            
            self.ELEM[eid].assembleStiffness(self.gk)
            
            self.ELEM[eid].computeElementForce()
            self.ELEM[eid].assembleLoad(self.gp)

        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            for fd in ['Fx', 'Fy', 'Fz',
                       'Mx', 'My', 'Mz']:
                index=self.force_index[fd]
                ii=N.dof[index]
                if fd in vars(N):
                    self.gp[ii-1][0]=vars(N)[fd]

    def setSupportCondition(self):
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            for jd in ['idx', 'idy', 'idz',
                       'irx', 'iry', 'irz']:
                if jd in vars(N):
                    idj=vars(N)[jd]
                    index=self.dof_index[jd]
                    ii=N.dof[index]
                    if idj == 0:
                        for j in range(self.nndof):
                            self.gk[ii-1][j]=0.0
                            self.gk[j][ii-1]=0.0
                           
                        self.gk[ii-1][ii-1]=1.0
                        
                        self.gp[ii-1][0]=0.0

    def solveStatic(self):
        print(self.gk)
        print(self.gp)
        self.disp=np.linalg.solve(self.gk,self.gp)
        print('disp:',self.disp)
        for nid in self.NODE_LIST:
            disp=[0,0,0,0,0,0]
            dof=self.NODE[nid].dof
            for i in range(6):
                ii=dof[i]
                if ii > 0:
                   disp[i]=self.disp[ii-1][0]
            print(dof,disp)
            self.NODE[nid].update(disp=disp)        

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMemberForces(self.NODE)
            

if __name__ == "__main__":

    def example_101():
        P=10000.0
        L=3.0
        A=150.0e-06
        Iz=1000.0e-09
        E=2.1e11
        rho=7850.0
        deltaP=P*L*L*L/(3.0*E*Iz)
        m_bar=rho*A

        pstr=Structure(title='Cantilever Beam with concentrated load at the free end')

        pstr.setDefaultDOF(idy=1,irz=1)

        pstr.node(nid=1,snid='A')
        pstr.node(nid=2,snid='B',x=1.5)
        pstr.node(nid=3,snid='C',x=3.0)

        pstr.node(nid=1,idy=0,irz=0)

        pstr.node(nid=3,Fy=-10000.0)
##
        pstr.element(eid=1,seid='AB',etype='Beam',nodes=(1,2))
        pstr.element(eid=2,seid='BC',nodes=(2,3))
 
        pstr.element(eid=1,Iz=1000.0e-09,E=2.1e11,m_bar=m_bar)
        pstr.element(eid=2,Iz=1000.0e-09,E=2.1e11,m_bar=m_bar)
##
        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
##
        pstr.printResults()
##
        print('Displacement due to concentrated load:',deltaP)
        
    def example_102():
        w=1000.0
        L=3.0
        A=150.0e-06
        Iz=1000.0e-09
        rho=7850.0
        E=2.1e11
        m_bar=rho*A
        
        deltaw=w*L*L*L*L/(8.0*E*Iz)    

        pstr=Structure(title='Cantilever Beam with UDL')

        pstr.setDefaultDOF(idy=1,irz=1)

        pstr.node(nid=1,snid='A')
        pstr.node(nid=2,snid='B',x=1.5)
        pstr.node(nid=3,snid='C',x=3.0)

        pstr.node(nid=1,idy=0,irz=0)
##
        pstr.element(eid=1,seid='AB',etype='Beam',nodes=(1,2))
        pstr.element(eid=2,seid='BC',nodes=(2,3))
 
        pstr.element(eid=1,Iz=1000.0e-09,E=2.1e11,m_bar=m_bar,w=-1000.0)
        pstr.element(eid=2,Iz=1000.0e-09,E=2.1e11,m_bar=m_bar,w=-1000.0)
##
        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
##
        pstr.printResults()
##

        print('Displacement due to UDL:',deltaw)

    def example_201():
        alpha=['A','B','C','D','E','F','G','H','I','J','K','L','M',
               'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        L=3.0
        A=150.0e-06
        Iz=1000.0e-09
        rho=7850.0 
        E=2.1e11
        P=10000.0
        deltaP=P*L*L*L/(48.0*E*Iz)
        m_bar=rho*A
        
        nnodes=11
        nelems=nnodes-1
        delx=L/(nnodes-1)
        nP=int((nnodes+1)/2)
        
        pstr=Structure(title='Simply Supported Beam with concentrated load at the mid point')

        pstr.setDefaultDOF(idy=1,irz=1)

        for i in range(nnodes):
            x=i*delx
            pstr.node(nid=i+1,snid=alpha[i],x=x,y=0.0)

        pstr.node(nid=nP,Fy=-P)
        pstr.node(nid=1,idy=0)
        pstr.node(nid=nnodes,idy=0)
        
        for i in range(nelems):
             pstr.element(eid=i+1,seid=alpha[i]+alpha[i+1],etype='Beam',nodes=(i+1,i+2),Iz=Iz,E=E,m_bar=m_bar)

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
##
        pstr.printResults()
##
        print('Displacement due to concentrated load at midspan:',deltaP)
        
    def example_202():
        alpha=['A','B','C','D','E','F','G','H','I','J','K','L','M',
               'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        L=3.0
        A=150.0e-06
        Iz=1000.0e-09
        rho=7850.0 
        E=2.1e11
        w=1000.0
        deltaw=5.0*w*L*L*L*L/(384.0*E*Iz)
        m_bar=rho*A
        
        nnodes=11
        nelems=nnodes-1
        delx=L/(nnodes-1)
        nP=int((nnodes+1)/2)
        
        pstr=Structure(title='Simply Supported Beam with concentrated load at the mid point')

        pstr.setDefaultDOF(idy=1,irz=1)

        for i in range(nnodes):
            x=i*delx
            pstr.node(nid=i+1,snid=alpha[i],x=x,y=0.0)

        pstr.node(nid=1,idy=0)
        pstr.node(nid=nnodes,idy=0)
        
        for i in range(nelems):
             pstr.element(eid=i+1,seid=alpha[i]+alpha[i+1],etype='Beam',nodes=(i+1,i+2),Iz=Iz,E=E,m_bar=m_bar,w=-w)

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
##
        pstr.printResults()
##
        print('Displacement due to UDL:',deltaw)
        
    example_101()
#    example_102()
#    example_201()
#    example_202()       
