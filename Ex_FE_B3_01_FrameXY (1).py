'''
  File : Ex_B3_01_FrameXY.py
  This program is written for static and dynamic analysis of frames in XY plane.
  
  Author : Prof. Nirjhar Dhang
           Department of Civil Engineering
           Indian Institute of Technology Kharagpur
           Kharagpur-721302
           West Bengal
           India
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

    def assign(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value

class Unit(Object):
    def __init__(self,**kwargs):
        self.F='N'
        self.M='Nm'
        self.E='N/mm^2'
        self.rho='kg/m^3'
        self.mbar='kg/m'
        self.w='N/m'
        self.P='N'
        self.Ix='m^4'
        self.Iy='m^4'
        self.G='N/mm^2'
        Object.__init__(self,**kwargs)
        
## Node class : for storing node related data            
class Node(Object):
    def __init__(self,**kwargs):
        self.nid=1
        self.snid=''
        self.x=0
        self.y=0
        self.z=0
        Object.__init__(self,**kwargs)
        
## Element class : for storing node related data          
class Element(Object):
    def __init__(self,**kwargs):
        Object.__init__(self,**kwargs)

class FrameXY(Element):
    def __init__(self,**kwargs):
        self.modelling='consistent'
        Element.__init__(self,**kwargs)
                
        self.mdof=6        
        self.ke=np.zeros((6,6),float)
        self.ok=np.zeros((6,6),float)
        self.om=np.zeros((6,6),float)
        
        self.ot=np.zeros((6,6),float)
        self.jdofv=np.zeros((6),int)
        
    def getLength(self,NODE):
        
        n1=self.nodes[0]
        n2=self.nodes[1]
        xx=NODE[n2].x-NODE[n1].x
        yy=NODE[n2].y-NODE[n1].y
        L=math.sqrt((xx*xx+yy*yy))

        return L, xx, yy
    
    def computeMatrices(self,NODE):
        
        self.cm=np.zeros((self.mdof,self.mdof),float)
       
        for i in range(self.mdof):
            for j in range(self.mdof):
                self.cm[i,j]=0.0

        for i in range(self.mdof):
            for j in range(self.mdof):
                self.ok[i,j]=0.0
                self.ot[i,j]=0.0
            self.ot[i,i]=1.0
           
        n1=self.nodes[0]
        n2=self.nodes[1]
        xx=NODE[n2].x-NODE[n1].x
        yy=NODE[n2].y-NODE[n1].y
        L=math.sqrt((xx*xx+yy*yy))
        self.L=L
        C = xx / L
        S = yy / L
        self.ot[0,0]=C
        self.ot[0,1]=S
        self.ot[1,0]=-S
        self.ot[1,1]=C
        self.ot[2,2]=1.0
        self.ot[3,3]=C
        self.ot[3,4]=S
        self.ot[4,3]=-S
        self.ot[4,4]=C
        self.ot[5,5]=1.0

        A=self.A
        Iz=self.Iz
        E=self.E

        print(n1,n2,L,A,E,Iz)

        AE_by_L = A * E / L

        self.ke[0,0]= AE_by_L
        self.ke[0,3]=-AE_by_L
    
        self.ke[3,0]=-AE_by_L
        self.ke[3,3]= AE_by_L

        EI_by_L3 =  E * Iz / (L*L*L)
        EI_by_L2 =  E * Iz / (L*L)
        EI_by_L =  E * Iz / L

        self.ke[1,1]= 12.0*EI_by_L3
        self.ke[1,2]= 6.0*EI_by_L2
        self.ke[1,4]= -12.0*EI_by_L3
        self.ke[1,5]= 6.0*EI_by_L2

        self.ke[2,1]= 6.0*EI_by_L2
        self.ke[2,2]= 4.0*EI_by_L
        self.ke[2,4]= -6.0*EI_by_L2
        self.ke[2,5]= 2.0*EI_by_L

        self.ke[4,1]= -12.0*EI_by_L3
        self.ke[4,2]= -6.0*EI_by_L2
        self.ke[4,4]= 12.0*EI_by_L3
        self.ke[4,5]= -6.0*EI_by_L2

        self.ke[5,1]= 6.0*EI_by_L2
        self.ke[5,2]= 2.0*EI_by_L
        self.ke[5,4]= -6.0*EI_by_L2
        self.ke[5,5]= 4.0*EI_by_L

        for i in range(self.mdof):
            for j in range(self.mdof):
                self.cm[i,j]=0.0
                for k in range(self.mdof):
                    temp=self.ot[k,i]*self.ke[k,j]                    
                    self.cm[i,j]=self.cm[i,j]+temp

        for i in range(self.mdof):
            for j in range(self.mdof):
                self.ok[i,j]=0.0
                for k in range(self.mdof):
                    temp=self.cm[i,k]*self.ot[k,j]
                    self.ok[i,j]=self.ok[i,j]+temp

        if self.modelling == 'consistent':
            mbar=self.mbar
            self.om[0,0]=140.0
            self.om[0,3]=70.0
            self.om[3,0]=70.0
            self.om[3,3]=140.0
            
            self.om[1,1]=156.0
            self.om[1,2]=22.0*L
            self.om[1,4]=54.0
            self.om[1,5]=-13.0*L
            self.om[2,1]=22.0*L
            self.om[2,2]=4.0*L*L
            self.om[2,4]=13.0*L
            self.om[2,5]=-3.0*L*L
            self.om[4,1]=54.0
            self.om[4,2]=13.0*L
            self.om[4,4]=156.0
            self.om[4,5]=-22.0*L
            self.om[5,1]=-13.0*L
            self.om[5,2]=-3.0*L*L
            self.om[5,4]=-22.0*L
            self.om[5,5]=4.0*L*L
            for i in range(6):
                for j in range(6):
                    self.om[i,j]=self.om[i,j]*mbar*L/420.0
            
        else:   
            alpha=1.0
            self.om[0,0]=mbar*L/2.0
            self.om[1,1]=mbar*L/2.0
            self.om[2,2]=alpha*mbar*L/2.0
            self.om[3,3]=mbar*L/2.0
            self.om[4,4]=mbar*L/2.0
            self.om[5,5]=alpha*mbar*L/2.0
            
        for i in range(self.mdof):
            for j in range(self.mdof):
                self.cm[i,j]=0.0
                for k in range(self.mdof):
                    self.cm[i,j]=self.cm[i,j]+self.ot[k,i]*self.om[k,j]

        for i in range(self.mdof):
            for j in range(self.mdof):
                self.om[i,j]=0.0
                for k in range(self.mdof):
                    self.om[i,j]=self.om[i,j]+self.cm[i,k]*self.ot[k,j]
                    
                    
        self.jdofv[0]=NODE[n1].dof[0]
        self.jdofv[1]=NODE[n1].dof[1]
        self.jdofv[2]=NODE[n1].dof[5]
        
        self.jdofv[3]=NODE[n2].dof[0]
        self.jdofv[4]=NODE[n2].dof[1]
        self.jdofv[5]=NODE[n2].dof[5]
        print(self.jdofv)
           
        del self.cm

    def computeElementForce(self):
            self.femfl=[0.0,0.0,0.0,0.0,0.0,0.0]
            self.femfg=[0.0,0.0,0.0,0.0,0.0,0.0]
            L=self.L
            if 'w' in vars(self):
                w=self.w
                self.femfl[1]=self.femfl[1]+0.5*w*L
                self.femfl[2]=self.femfl[2]+w*L*L/12.0
                self.femfl[4]=self.femfl[4]+0.5*w*L
                self.femfl[5]=self.femfl[5]-w*L*L/12.0
            if 'P' in vars(self):
                P=self.P
                a=self.a
                self.femfl[1]=self.femfl[1]+P*(L-a)/L
                self.femfl[2]=self.femfl[2]+P*(L-a)*(L-a)*a/(L*L)
                self.femfl[4]=self.femfl[4]+P*a/L
                self.femfl[5]=self.femfl[5]-P*a*a*(L-a)/(L*L)
            for i in range(self.mdof):
                for j in range(self.mdof):
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

    def assembleMassMatrix(self,gm):
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                for j in range(self.mdof):
                    jj=self.jdofv[j]
                    if jj > 0:
                        temp=gm[ii-1,jj-1]+self.om[i,j] 
                        gm[ii-1,jj-1]=temp

    def assembleLoad(self,gp):
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                temp=gp[ii-1][0]+self.femfg[i] 
                gp[ii-1][0]=temp        
                        
   
    
    def computeMemberForces(self,NODE):
            self.computeMatrices(NODE)
            mfg=[0,0,0,0,0,0]
            mfl=[0,0,0,0,0,0]
            disp=[0,0,0,0,0,0]            
            n1=self.nodes[0]
            n2=self.nodes[1]
            dispN1=NODE[n1].disp
            dispN2=NODE[n2].disp
            disp[0]=dispN1[0]
            disp[1]=dispN1[1]
            disp[2]=dispN1[5]
            disp[3]=dispN2[0]
            disp[4]=dispN2[1]
            disp[5]=dispN2[5]

            for i in range(self.mdof):
                mfg[i]=0.0
                for j in range(self.mdof):
                    mfg[i]=mfg[i]+self.ok[i,j]*disp[j]
                    
            for i in range(self.mdof):
                mfl[i]=0.0
                for j in range(self.mdof):
                    mfl[i]=mfl[i]+self.ot[i,j]*mfg[j]

            if 'femfl' in vars(self):
               femfl=self.femfl
                        
               for i in range(self.mdof):
                   mfl[i]=mfl[i]-femfl[i]                    

            self.mfl=mfl
            print('mfl:',mfl)

            
class Structure(Object):
    def __init__(self,**kwargs):
        self.modelling='consistent'
        self.assembly='compact'
        self.numnode=0
        self.numelem=0
        self.U=Unit(L='m',F='N',M='Nm',E='N/mm^2',G='N/mm^2',rho='kg/m^3',mbar='kg/m',A='m^2',Iz='m^4',w='N/m',P='N',disp='m',rot='rad')
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
        if 'SNODE' not in vars(self):
            self.SNODE=dict()
            
        if nid not in self.NODE_LIST:
            self.NODE_LIST.append(nid)
            self.NODE[nid]=Node(**kwargs)
            if 'snid' in vars(self.NODE[nid]):
               snid=self.NODE[nid].snid
               self.SNODE[snid]=self.NODE[nid]
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
        if 'SELEM' not in vars(self):
            self.SELEM=dict()
            
        if eid not in self.ELEM_LIST:
            self.ELEM_LIST.append(eid)
            if 'etype' in kwargs:
                self.etype=kwargs['etype']
            if self.etype == 'FrameXY':
                self.ELEM[eid]=FrameXY(**kwargs)
            else:
                self.ELEM[eid]=Element(**kwargs)
            if 'seid' in vars(self.ELEM[eid]):
               seid=self.ELEM[eid].seid
               self.SELEM[seid]=self.ELEM[eid]
                
            self.numelem=self.numelem+1
        else:
            self.ELEM[eid].update(**kwargs)

    def getLength(self,eid):
        E=self.ELEM[eid]
        return E.getLength(self.NODE)           

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
                  ' Nodes:(%d,%d)'%(E.nodes[0],E.nodes[1]), ' Area : %.4e'%E.A, ' Iz :',E.Iz,'E : %.4e'%E.E,
                  'MFL:',E.mfl[0],E.mfl[1],E.mfl[2],E.mfl[3],E.mfl[4],E.mfl[5])

    def decideActiveDOF(self):
        self.ndof=0
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            self.dof=[0,0,0,0,0,0]
            for jd in ['idx', 'idy', 'idz',
                       'irx', 'iry', 'irz']:
                if self.DEFAULT_DOF[jd] == 1:
                    if self.assembly == 'compact':
                        if jd in vars(N):
                            if vars(N)[jd] == 1:
                                self.ndof=self.ndof+1
                                index=self.dof_index[jd]
                                self.dof[index]=self.ndof
                        else:
                            self.ndof=self.ndof+1
                            index=self.dof_index[jd]
                            self.dof[index]=self.ndof
                    else:
                        self.ndof=self.ndof+1
                        index=self.dof_index[jd]
                        self.dof[index]=self.ndof
                        
                self.NODE[nid].update(dof=self.dof)

    def computeSystemMatrices(self,**kwargs):
            
        self.gk=np.zeros((self.ndof,self.ndof),float)
        self.gm=np.zeros((self.ndof,self.ndof),float)
        
        self.gp=np.zeros((self.ndof,1),float)

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMatrices(self.NODE)
            self.ELEM[eid].assembleStiffness(self.gk)
            self.ELEM[eid].assembleMassMatrix(self.gm)
            
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
        if self.assembly == 'full':
            for nid in self.NODE_LIST:
                N=self.NODE[nid]
                for jd in ['idx', 'idy', 'idz',
                           'irx', 'iry', 'irz']:
                    if jd in vars(N):
                        idj=vars(N)[jd]
                        index=self.dof_index[jd]
                        ii=N.dof[index]
                        if idj == 0:
                            for j in range(self.ndof):
                                self.gk[ii-1][j]=0.0
                                self.gk[j][ii-1]=0.0
                                self.gm[ii-1][j]=0.0
                                self.gm[j][ii-1]=0.0
                                
                            self.gk[ii-1][ii-1]=1.0
                            self.gm[ii-1][ii-1]=1.0
                            
                            self.gp[ii-1][0]=0.0

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
            
    def solveEigen(self):
         print ('Eigen value analysis started')
         self.Lgm=la.cholesky(self.gm)
         self.Lgmi=la.inv(self.Lgm)
         self.LiK=np.dot(self.Lgmi,self.gk)
         self.K_bar=np.dot(self.LiK,self.Lgmi.T)
         self.EigenValues,self.EigenVectorsL=la.eig(self.K_bar)
         print ('Eigen value analysis over')

         
#         Ev=np.dot(la.inv(Lgm.T),self.EigenVectors)
         self.EigenVectors=np.dot(self.Lgmi.T,self.EigenVectorsL)
         self.EigenVectorsMass=np.dot(self.Lgmi.T,self.EigenVectors)
         
         '''
         self.gmi=la.inv(self.gm)
         self.K_bar=np.dot(self.gmi,self.gk)
         self.EigenValues,self.EigenVectors=la.eig(self.K_bar)
         '''

         self.omega=np.zeros(self.ndof,float)
         self.frequency=np.zeros(self.ndof,float)
         self.timeperiod=np.zeros(self.ndof,float)

         for i in range(self.ndof):
              if self.EigenValues[i] > 0.0:
                  self.omega[i]=math.sqrt(self.EigenValues[i])
              else:
                  self.omega[i]='Neg'
                  print( 'Eigen values negative')
                  self.errormessage.append('Error : Eigen values negative')
              self.frequency[i]=self.omega[i]/(2.0*3.14159)
              if self.omega[i] > 0.0:
                  self.timeperiod[i]=2.0*3.14159/self.omega[i]
              else:
                  self.timeperiod[i]='Neg'
                  self.errormessage.append('Warning : Omega [%d]=%.4e'%(i,self.omega[i]))

    def saveAsDisplacement(self,evec):
        for nid in self.NODE_LIST:
            disp=[0,0,0,0,0,0]
            dof=self.NODE[nid].dof
            for i in range(6):
                ii=dof[i]
                if ii > 0:
                   disp[i]=evec[ii-1]
            self.NODE[nid].update(evdisp=disp)
            
    def getEigenValue(self,**kwargs):
        attribs={}
        for key,value in kwargs.items():
            attribs[key]=value
        if 'index' in attribs:
            index=attribs['index']
        else:
            index=1
        sortindex=self.EigenValues.argsort()
        return self.EigenValues[sortindex[index-1]]
    
    def getEigenVector(self,**kwargs):
        attribs={}
        for key,value in kwargs.items():
            attribs[key]=value
        if 'index' in attribs:
            index=attribs['index']
        else:
            index=1            
        sortindex=self.EigenValues.argsort()

        return self.EigenVectors[:,sortindex[index-1]]

    
    def csvStaticNodalTable(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value
       
        if 'fcsvnode' not in vars(self):
            if 'fileid' in vars(self):
                self.csvnodefilename=self.fileid+'_node.csv'
            self.fcsvnode=open(self.csvnodefilename,'w')
        
        self.fcsvnode.write("%s\n"%self.title)
        self.fcsvnode.write("Table-1 : Nodal Data and Results\n")
        self.fcsvnode.write("recno,nid,snid,x,y,z,idx,idy,idz,irx,iry,irz,lcase,Fx,Fy,Fz,Mx,My,Mz,Dx,Dy,Dz,Rx,Ry,Rz\n")
        self.fcsvnode.write(",,,%s,%s,%s,,,,,,,,"%(self.U.L,self.U.L,self.U.L)+
                            "%s,%s,%s,%s,%s,%s,"%(self.U.F,self.U.F,self.U.F,self.U.M,self.U.M,self.U.M)+
                            "%s,%s,%s,%s,%s,%s\n"%(self.U.disp,self.U.disp,self.U.disp,self.U.rot,self.U.rot,self.U.rot))
                                                        
        recno=0 
        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            Fx=Fy=Fz=Mx=My=Mz=0.0
            Dx=Dy=Dz=Rx=Ry=Rz=0.0
            idx=idy=irz=1
            idz=irx=iry=0
            Dx=N.disp[0]
            Dy=N.disp[1]
            Rz=N.disp[5]
            supportnode=0
            if 'idx' in vars(N):
                idx=N.idx
                supportnode=supportnode+1
            if 'idy' in vars(N):
                idy=N.idy
                supportnode=supportnode+1
            if 'irz' in vars(N):
                irz=N.irz
                supportnode=supportnode+1
                
            loadednode=0
            if 'Fx' in vars(N):
                Fx=N.Fx
                loadednode=loadednode+1
            if 'Fy' in vars(N):
                Fy=N.Fy
                loadednode=loadednode+1
            if 'Mz' in vars(N):
                Mz=N.Mz
                loadednode=loadednode+1


            recno=recno+1
            self.fcsvnode.write("%d,%d,%s,%.3f,%.3f,%.3f,%d,%d,%d,%d,%d,%d,"%(recno,N.nid,N.snid,N.x,N.y,N.z,idx,idy,idz,irx,iry,irz)+
                                "%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"%(self.lcase,Fx,Fy,Fz,Mx,My,Mz)+
                                "%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n"%(Dx,Dy,Dz,Rx,Ry,Rz))

        self.fcsvnode.close()
        
                        
    def csvStaticElementTable(self,**kwargs):
        for key, value in kwargs.items():
            vars(self)[key]=value
        
        if 'fcsvelem' not in vars(self):
            if 'fileid' in vars(self):
                self.csvelemfilename=self.fileid+'_elem.csv'
            self.fcsvelem=open(self.csvelemfilename,'w')
        
        self.fcsvelem.write("%s\n"%self.title)
        
        self.fcsvelem.write("Table-2 : Element Data and Results\n")
        self.fcsvelem.write("recno,eid,seid,E,G,Nu,rho,A,Ix,Iy,Iz,mbar,lcase,w,P,a,nidpos,nid,Fx,Fy,Fz,Mx,My,Mz\n")
        self.fcsvelem.write(",,,%s,%s,,%s,"%(self.U.E,self.U.G,self.U.rho)+
                            "%s,%s,%s,%s,%s,"%(self.U.A,self.U.Ix,self.U.Iy,self.U.Iz,self.U.mbar)+
                            ",%s,%s,%s,"%(self.U.w,self.U.F,self.U.L)+
                            ",,%s,%s,%s,%s,%s,%s\n"%(self.U.F,self.U.F,self.U.F,self.U.M,self.U.M,self.U.M))
        recno=0
        for eid in self.ELEM_LIST:
            E=self.ELEM[eid]
            G=Nu=Ix=Iy=w=P=a=0.0
            Fx1=Fy1=Fz1=Mx1=My1=Mz1=0.0
            Fx2=Fy2=Fz2=Mx2=My2=Mz2=0.0
            if 'G' in vars(E):
                G=E.G
            if 'Nu' in vars(E):
                Nu=E.Nu
            if 'Ix' in vars(E):
                Ix=E.Ix
            if 'Iy' in vars(E):
                Iy=E.Iy
            if 'w' in vars(E):
                w=E.w
            if 'P' in vars(E):
                P=E.P
            if 'a' in vars(E):
                a=E.a
            Fx1=E.mfl[0]
            Fy1=E.mfl[1]
            Mz1=E.mfl[2]
            Fx2=E.mfl[3]
            Fy2=E.mfl[4]
            Mz2=E.mfl[5]
            nidpos1=1
            nidpos2=2
            recno=recno+1
            self.fcsvelem.write("%d,%d,%s,%.3e,%.3e,%.3f,%.1f,"%(recno,E.eid,E.seid,E.E,G,Nu,E.rho)+
                                "%.4e,%.4e,%.4e,%.4e,%4e,"%(E.A,Ix,Iy,E.Iz,E.mbar)+
                                "%d,%.3f,%.3f,%.3f,"%(self.lcase,w,P,a)+
                                "%d,%d,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n"%(nidpos1,E.nodes[0],Fx1,Fy1,Fz1,Mx1,My1,Mz1))
            recno=recno+1
            self.fcsvelem.write("%d,%d,%s,%.3e,%.3e,%.3f,%.1f,"%(recno,E.eid,E.seid,E.E,G,Nu,E.rho)+
                                "%.4e,%.4e,%.4e,%.4e,%4e,"%(E.A,Ix,Iy,E.Iz,E.mbar)+
                                "%d,%.3f,%.3f,%.3f,"%(self.lcase,w,P,a)+
                                "%d,%d,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n"%(nidpos2,E.nodes[1],Fx2,Fy2,Fz2,Mx2,My2,Mz2))
            
        self.fcsvelem.close()          


if __name__ == "__main__":

    def example_01():
        rho=7850.0
        A=150.0e-06
        Iz=10000.0e-09
        E=2.1e11
        w=-1000.0
        mbar=A*rho
    
        pstr=Structure(title='Cantilever Beam')

        pstr.setDefaultDOF(idx=1,idy=1,irz=1)

        pstr.node(nid=1,snid='A')
        pstr.node(nid=2,snid='B',x=1.5)
        pstr.node(nid=3,snid='C',x=3.0)

        pstr.node(nid=1,idx=0,idy=0,irz=0)

#        pstr.node(nid=3,Fy=-10000.0)
##
        pstr.element(eid=1,seid='AB',etype='FrameXY',nodes=(1,2))
        pstr.element(eid=2,seid='BC',nodes=(2,3))
 
        pstr.element(eid=1,A=150.0e-06,Iz=1000.0e-09,E=2.1e11,rho=rho,mbar=mbar,w=-1000.0)
        pstr.element(eid=2,A=150.0e-06,Iz=1000.0e-09,E=2.1e11,rho=rho,mbar=mbar,w=-1000.0)
##

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

#        pstr.setSupportCondition()

        pstr.solveStatic()
        pstr.solveEigen()
        print(pstr.EigenValues)

##

        pstr.printResults()

##
        w=1000.0
        P=10000.0
        L=3.0
        Iz=1000.0e-09
        E=2.1e11
        delta=P*L*L*L/(3.0*E*Iz)
        deltaw=w*L*L*L*L/(8.0*E*Iz)
        print('Displacement:',delta,deltaw)


    def example_02():
        rho=7850.0
        A=150.0e-06
        Iz=10000.0e-09
        E=2.1e11
        w=-1000.0
        P=-10000.0
        mbar=A*rho
        L=3.0
        nelem=10
        lx=L/nelem
    
        pstr=Structure(title='Cantilever Beam')

        pstr.setDefaultDOF(idx=1,idy=1,irz=1)

        for i in range(nelem+1):
            pstr.node(nid=i+1,snid=chr(65+i),x=i*lx)

        pstr.node(nid=1,idx=0,idy=0,irz=0)

#        pstr.node(nid=3,Fy=-10000.0)
##
        for i in range(nelem):
            pstr.element(eid=i+1,seid=chr(65+i)+chr(65+i+1),etype='FrameXY',nodes=(i+1,i+2),E=E,A=A,Iz=Iz,rho=rho,mbar=mbar,w=w)
##

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

#        pstr.setSupportCondition()

        pstr.solveStatic()
        pstr.solveEigen()
        print(pstr.EigenValues)

##

        pstr.printResults()

##
        delta=P*L*L*L/(3.0*E*Iz)
        deltaw=w*L*L*L*L/(8.0*E*Iz)
        print('Displacement:',delta,deltaw)
        
    def example_03():
        rho=2500.0 # kg/m^3
        bc=300.0 # mm
        Dc=300.0 # mm
        bb=250.0 # mm
        Db=350.0 # mm
        fck=25.0 # N/mm^2
        E=5000.0*math.sqrt(fck)*1.0e06 # N/m^2
        Ac=bc*Dc*1.0e-06
        Izc=bc*Dc*Dc*Dc*1.0e-12/12.0
        Ab=bb*Db*1.0e-06
        Izb=bb*Db*Db*Db*1.0e-12/12.0
        w=-35000.0 # N/m
        mbarc=Ac*rho
        mbarb=Ab*rho
        L=5.0
        H=3.0 # Floor height
    
        pstr=Structure(title='Portal frame',assembly='compact',lcase=1)

        pstr.setDefaultDOF(idx=1,idy=1,irz=1)

        pstr.node(nid=1,snid='A0')
        pstr.node(nid=2,snid='C0',x=L)
        pstr.node(nid=3,snid='A1',x=0.0,y=H)
        pstr.node(nid=4,snid='B1',x=0.5*L,y=H)
        pstr.node(nid=5,snid='C1',x=L,y=H)

        pstr.node(nid=1,idx=0,idy=0,irz=0)
        pstr.node(nid=2,idx=0,idy=0,irz=0)

#        pstr.node(nid=3,Fy=-10000.0)
##
        pstr.element(eid=1,seid='A0-A1',etype='FrameXY',nodes=(1,3),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=2,seid='C0-C1',nodes=(2,5),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=3,seid='A1-B1',nodes=(3,4),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb)
        pstr.element(eid=4,seid='A1-C1',nodes=(4,5),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb)

##

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
        pstr.solveEigen()
        print(pstr.EigenValues)

##

        pstr.printResults()


    def example_04():
        rho=2500.0 # kg/m^3
        bc=300.0 # mm
        Dc=300.0 # mm
        bb=250.0 # mm
        Db=350.0 # mm
        fck=25.0 # N/mm^2
        E=5000.0*math.sqrt(fck)*1.0e06 # N/m^2
        Ac=bc*Dc*1.0e-06
        Izc=bc*Dc*Dc*Dc*1.0e-12/12.0
        Ab=bb*Db*1.0e-06
        Izb=bb*Db*Db*Db*1.0e-12/12.0
        w=-35000.0 # N/m
        mbarc=Ac*rho
        mbarb=Ab*rho
        L=5.0
        H0=1.2 # Plinth height
        H1=3.0 # First floor height
        H2=3.0 # Second floor height
        H3=3.0 # Third floor height
        H4=3.0 # Fourth floor height
    
        pstr=Structure(title='Four storeyed single bay building frame',assembly='compact',lcase=1)

        pstr.setDefaultDOF(idx=1,idy=1,irz=1)

        pstr.node(nid=1,snid='A0')
        pstr.node(nid=2,snid='C0',x=L)
        pstr.node(nid=3,snid='AP',x=0.0,y=H0)
        pstr.node(nid=4,snid='BP',x=0.5*L,y=H0)
        pstr.node(nid=5,snid='CP',x=L,y=H0)
        pstr.node(nid=6,snid='A1',x=0.0,y=H0+H1)
        pstr.node(nid=7,snid='B1',x=0.5*L,y=H0+H1)
        pstr.node(nid=8,snid='C1',x=L,y=H0+H1)
        pstr.node(nid=9,snid='A2',x=0.0,y=H0+H1+H2)
        pstr.node(nid=10,snid='B2',x=0.5*L,y=H0+H1+H2)
        pstr.node(nid=11,snid='C2',x=L,y=H0+H1+H2)
        pstr.node(nid=12,snid='A3',x=0.0,y=H0+H1+H2+H3)
        pstr.node(nid=13,snid='B3',x=0.5*L,y=H0+H1+H2+H3)
        pstr.node(nid=14,snid='C3',x=L,y=H0+H1+H2+H3)
        pstr.node(nid=15,snid='A4',x=0.0,y=H0+H1+H2+H3+H4)
        pstr.node(nid=16,snid='B4',x=0.5*L,y=H0+H1+H2+H3+H4)
        pstr.node(nid=17,snid='C4',x=L,y=H0+H1+H2+H3+H4)
        

        pstr.node(nid=1,idx=0,idy=0,irz=0)
        pstr.node(nid=2,idx=0,idy=0,irz=0)

#        pstr.node(nid=3,Fy=-10000.0)
##
        pstr.element(eid=1,seid='A0-AP',etype='FrameXY',nodes=(1,3),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=2,seid='AP-A1',nodes=(3,6),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=3,seid='A1-A2',nodes=(6,9),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=4,seid='A2-A3',nodes=(9,12),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=5,seid='A3-A4',nodes=(12,15),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        
        pstr.element(eid=6,seid='C0-CP',nodes=(2,5),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=7,seid='CP-C1',nodes=(5,8),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=8,seid='C1-C2',nodes=(8,11),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=9,seid='C2-C3',nodes=(11,14),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)
        pstr.element(eid=10,seid='C3-C4',nodes=(14,17),A=Ac,Iz=Izc,E=E,rho=rho,mbar=mbarc)

        pstr.element(eid=11,seid='AP-BP',nodes=(3,4),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        pstr.element(eid=12,seid='BP-CP',nodes=(4,5),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        
        pstr.element(eid=13,seid='A1-B1',nodes=(6,7),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        pstr.element(eid=14,seid='B1-C1',nodes=(7,8),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)

        pstr.element(eid=15,seid='A2-B2',nodes=(9,10),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        pstr.element(eid=16,seid='B2-C2',nodes=(10,11),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        
        pstr.element(eid=17,seid='A3-B3',nodes=(12,13),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        pstr.element(eid=18,seid='B3-C3',nodes=(13,14),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        
        pstr.element(eid=19,seid='A4-B4',nodes=(15,16),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)
        pstr.element(eid=20,seid='B4-C4',nodes=(16,17),A=Ab,Iz=Izb,E=E,rho=rho,mbar=mbarb,w=w)

##

        pstr.decideActiveDOF()

        pstr.computeSystemMatrices()

        pstr.setSupportCondition()

        pstr.solveStatic()
        pstr.solveEigen()
        print(pstr.EigenValues)

##

        pstr.printResults()
        pstr.csvStaticNodalTable(fileid='four_storeyed_single_bay')
        pstr.csvStaticElementTable(fileid='four_storeyed_single_bay')

        for i in range(1,6):
            eigenvalue=pstr.getEigenValue(index=i)
            omega=math.sqrt(eigenvalue)            
            frequency=omega/(2.0*3.14159)
            timeperiod=2.0*3.14159/omega            
            print('Eigen value: %2d,%.3f,%.3f,%.3f,%.3f'%(i,eigenvalue,omega,frequency,timeperiod))
    example_04()
        
