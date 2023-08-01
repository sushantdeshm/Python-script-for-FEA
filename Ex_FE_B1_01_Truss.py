from __future__ import print_function

import numpy as np
import math
class Node(object):
    def __init__(self,**kwargs):
        self.nid=1
        self.snid=''
        self.x=0
        self.y=0
        self.z=0
        for key, value in kwargs.items():
            vars(self)[key]=value

    def update(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value
##
class Element(object):
    def __init__(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value

    def update(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value
##
class TrussXY(Element):
    def __init__(self,**kwargs):
        Element.__init__(self,**kwargs)
                
        self.mdof=4        
        self.ke=np.zeros((4,4),float)
        self.ok=np.zeros((4,4),float)
        
        self.ot=np.zeros((4,4),float)
        self.jdofv=np.zeros((4),int)

##                         
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
##                
        n1=self.nodes[0]
        n2=self.nodes[1]
        xx=NODE[n2].x-NODE[n1].x
        yy=NODE[n2].y-NODE[n1].y
        L=math.sqrt((xx*xx+yy*yy))
        C = xx / L
        S = yy / L
        self.ot[0,0]=C
        self.ot[0,1]=S
        self.ot[1,0]=-S
        self.ot[1,1]=C
        self.ot[2,2]=C
        self.ot[2,3]=S
        self.ot[3,2]=-S
        self.ot[3,3]=C
##
        A=self.A
        E=self.E

        AE_by_L = A * E / L

        self.ke[0,0]= AE_by_L
        self.ke[0,2]=-AE_by_L
    
        self.ke[2,0]=-AE_by_L
        self.ke[2,2]= AE_by_L
##
        for i in range(4):
            for j in range(4):
                self.c[i,j]=0.0
                for k in range(4):
                    temp=self.ot[k,i]*self.ke[k,j]                    
                    self.c[i,j]=self.c[i,j]+temp

        for i in range(4):
            for j in range(4):
                self.ok[i,j]=0.0
                for k in range(4):
                    temp=self.c[i,k]*self.ot[k,j]
                    self.ok[i,j]=self.ok[i,j]+temp
##
        self.jdofv[0]=NODE[n1].dof[0]
        self.jdofv[1]=NODE[n1].dof[1]
        self.jdofv[2]=NODE[n2].dof[0]
        self.jdofv[3]=NODE[n2].dof[1]                    
           
        del self.c
##
    def assembleStiffness(self,gk):
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                for j in range(self.mdof):
                    jj=self.jdofv[j]
                    if jj > 0:
                        temp=gk[ii-1,jj-1]+self.ok[i,j] 
                        gk[ii-1,jj-1]=temp        
##
    def computeMemberForces(self,NODE):
            self.computeMatrices(NODE)
            mfg=[0,0,0,0]
            mfl=[0,0,0,0]
            disp=[0,0,0,0]            
            n1=self.nodes[0]
            n2=self.nodes[1]
            dispN1=NODE[n1].disp
            dispN2=NODE[n2].disp
            disp[0]=dispN1[0]
            disp[1]=dispN1[1]
            disp[2]=dispN2[0]
            disp[3]=dispN2[1]
##            
            for i in range(4):
                mfg[i]=0.0
                for j in range(4):
                    mfg[i]=mfg[i]+self.ok[i,j]*disp[j]
                    
            for i in range(4):
                mfl[i]=0.0
                for j in range(4):
                    mfl[i]=mfl[i]+self.ot[i,j]*mfg[j]
##                        
            if 'femfl' in vars(self):
               femfl=self.femfl
                        
               for i in range(4):
                   mfl[i]=mfl[i]-femfl[i]                    

            self.mfl=mfl
            
##
class Structure(object):
    def __init__(self,**kwargs):
        self.numnode=0
        self.numelem=0
        self.DEFAULT_DOF=dict(idx=0,idy=0,idz=0,
                              irx=0,iry=0,irz=0)

        self.dof_index=dict(idx=0,idy=1,idz=2,
                            irx=3,iry=4,irz=5)

        self.force_index=dict(Fx=0,Fy=1,Fz=2,
                              Mx=3,My=4,Mz=5)
        for key,value in kwargs.items():
            vars(self)[key]=value
##
    def setDefaultDOF(self,**kwargs):
        for key, value in kwargs.items():
            self.DEFAULT_DOF[key] = value
##               
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
##
    def element(self,**kwargs):
        if 'eid' in kwargs:
            eid=kwargs['eid']
        if 'ELEM_LIST' not in vars(self):
            self.ELEM_LIST=list()
        if 'ELEM' not in vars(self):
            self.ELEM=dict()
##            
        if eid not in self.ELEM_LIST:
            self.ELEM_LIST.append(eid)
            if 'etype' in kwargs:
                self.etype=kwargs['etype']
            if self.etype == 'TrussXY':
                self.ELEM[eid]=TrussXY(**kwargs)
            else:
                self.ELEM[eid]=Element(**kwargs)
            self.numelem=self.numelem+1
        else:
            self.ELEM[eid].update(**kwargs)
##
    def printNodes(self):
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            print('nid:',N.nid,'snid:',N.snid,
                  'x:',N.x,'y:',N.y,'z:',N.z)
##
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            idx=1
            idy=1
            supportnode=0
            if 'idx' in vars(N):
                idx=N.idx
                supportnode=supportnode+1
            if 'idy' in vars(N):
                idy=N.idy
                supportnode=supportnode+1
            if supportnode > 0:
                print('nid:',nid,'idx=',idx,'idy=',idy)
##
        for nid in self.NODE_LIST:
            N=self.NODE[nid]

            Fx=0.0
            Fy=0.0
            loadednode=0
            if 'Fx' in vars(N):
                Fx=N.Fx
                loadednode=loadednode+1
            if 'Fy' in vars(N):
                Fy=N.Fy
                loadednode=loadednode+1

            if loadednode > 0:
               print('nid:',nid,'Fx=',Fx,'Fy=',Fy)
##
    def printElements(self):
        for eid in self.ELEM_LIST:
            E=self.ELEM[eid]
            print('eid:',E.eid,'seid:',E.seid,
                  'Nodes:(%d,%d)'%(E.nodes[0],E.nodes[1]))
##
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
##
    def computeSystemMatrices(self,**kwargs):
            
        self.gk=np.zeros((self.nndof,self.nndof),float)
        self.gp=np.zeros((self.nndof,1),float)

        print('initial gk:',self.gk)

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMatrices(self.NODE)
            self.ELEM[eid].assembleStiffness(self.gk)
            print('after assembly of elem:',eid)
            print('gk:',self.gk)
##
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            for fd in ['Fx', 'Fy', 'Fz',
                       'Mx', 'My', 'Mz']:
                index=self.force_index[fd]
                ii=N.dof[index]
                if fd in vars(N):
                    self.gp[ii-1][0]=vars(N)[fd]

        print('gp:',self.gp)
##
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
##            
    def solveStatic(self):
        self.disp=np.linalg.solve(self.gk,self.gp)

        for nid in self.NODE_LIST:
            disp=[0,0,0,0,0,0]
            dof=self.NODE[nid].dof
            for i in range(6):
                ii=dof[i]
                if ii > 0:
                   disp[i]=self.disp[ii-1][0]
            self.NODE[nid].update(disp=disp)        

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMemberForces(self.NODE)
##            

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

class PlotStructure(object):
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            vars(self)[key]=value
##
    def show(self,**kwargs):
        for key,value in kwargs.items():
            vars(self)[key]=value
        plt.axis((-2.0,6.0,-2.0,6.0))
        ax=plt.gca()
        plt.axis('off')
        for eid in self.strobj.ELEM_LIST:
            n1=self.strobj.ELEM[eid].nodes[0]
            n2=self.strobj.ELEM[eid].nodes[1]
            N1=self.strobj.NODE[n1]
            N2=self.strobj.NODE[n2]
            l=mlines.Line2D([N1.x,N2.x],[N1.y,N2.y])
            ax.add_line(l)
            
        plt.show()
##
if __name__ == "__main__":

    def example_01():
    
        pstr=Structure(title='Three Bar Truss')

        pstr.setDefaultDOF(idx=1,idy=1)

        pstr.node(nid=1,snid='A')
        pstr.node(nid=2,snid='B',x=3.0)
        pstr.node(nid=3,snid='C',x=1.0,y=3.0)

        pstr.node(nid=1,idx=0,idy=0)
        pstr.node(nid=2,idy=0)

        pstr.node(nid=3,Fy=-10000.0)
##
        pstr.element(eid=1,seid='AB',etype='TrussXY',nodes=(1,2))
        pstr.element(eid=2,seid='AC',nodes=(1,3))
        pstr.element(eid=3,seid='BC',nodes=(2,3))

        pstr.element(eid=1,A=150.0e-06,E=2.1e11)
        pstr.element(eid=2,A=150.0e-06,E=2.1e11)
        pstr.element(eid=3,A=150.0e-06,E=2.1e11)
##

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()

##

        pstr.printNodes()
        pstr.printElements()

        for nid in pstr.NODE_LIST:
            N=pstr.NODE[nid]
            print(N.nid,N.disp)

        for eid in pstr.ELEM_LIST:
            E=pstr.ELEM[eid]
            print(E.eid,E.mfl)

##
        plotstr=PlotStructure(strobj=pstr)
        plotstr.show()
        
    example_01()
        
