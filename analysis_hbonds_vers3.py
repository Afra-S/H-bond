
from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math


numb=int(input("number of traj  "));
typetop=input("Enter type of trajectory [pdb/xtc] ");

if(typetop=='pdb'):
    name=input("Enter name trajctory: ");
    traj = md.load(name+'.pdb')
elif(typetop=='xtc'):
    name=input("Enter name trajctory: ");
    topname=input("Enter name topology with extension (ex start.pdb): ");

ref = md.load(topname)
topref= ref.topology

Listname=[]
for ii,rr in enumerate(topref.residues):
    Listname.append(rr.name)


#print(nsize)

h,w = 501000, 250;
tothbonds=np.zeros((h, w))
tothbondsWC=np.zeros((h, w))
tothbondsNWC=np.zeros((h, w))
listrhbonds=np.zeros((h, w))
reshbonds=np.zeros((h, w))
reshbondsWC=np.zeros((h, w))
reshbondsNWC=np.zeros((h, w))

###For WC HB
ListA=['N6','N1']
ListC=['N4','N3','O2']
ListU=['O4','N3']
ListG=['N1','N2','O6']

DicAtoms = {'A': ListA, 'C': ListC,'G': ListG,'U': ListU}
ntframes=0
iloc=-1
for i in range(0,numb):
    traj = md.load(name+'_'+str(i)+'.xtc',top=topname)
    print(traj)
    topology = traj.topology
    nres=traj.n_residues
    print(nres)
    print(topology)
                    
    #hbonds = md.baker_hubbard(traj, periodic=False,freq=0.05)
    hbonds=md.wernet_nilsson(traj)

    ntframesloc=traj.n_frames

    ntframes=ntframes+ntframesloc
    
    for t in range(0,ntframesloc):
        iloc=iloc+1
        sizear=len(hbonds[t])
        #    print(hbonds[i])
        listres=np.full((250, 250), True)
        listresWC=np.full((250, 250), True)
        listresNWC=np.full((250, 250), True)
    
    
    
        for j in range(0,sizear):
            #print(i,j,hbonds[i][j][0],hbonds[i][j][2],traj.topology.atom(hbonds[i][j][2]).residue.index)
            hb1=traj.topology.atom(hbonds[t][j][0])
            hb2=traj.topology.atom(hbonds[t][j][2])
            indx=hb1.residue.index
            indy=hb2.residue.index
            namehb1=hb1.name
            namehb2=hb2.name
            nameres1=hb1.residue.name
            nameres2=hb2.residue.name
            lim1=len(DicAtoms[nameres1])
            lim2=len(DicAtoms[nameres2])
            logwc1=False
            logwc2=False
            for ll in range(0,lim1):
                if(DicAtoms[nameres1][ll]==namehb1):
                    logwc1=True
                    #print(indx,DicAtoms[nameres1][ll],namehb1)
                    break
            for ll in range(0,lim2):
                if(DicAtoms[nameres2][ll]==namehb2):
                    logwc2=True
                    #print(indx,DicAtoms[nameres1][ll],namehb1)
                    break
            if(logwc1==True and logwc2==True):
                tothbondsWC[iloc][indx]=tothbondsWC[iloc][indx]+1
                tothbondsWC[iloc][indy]=tothbondsWC[iloc][indy]+1
                if(listresWC[indx,indy]==True and listresWC[indy,indx]==True):
                    #           print(listres[indx,indy])
                    reshbondsWC[iloc][indx]=reshbondsWC[iloc][indx]+1
                    reshbondsWC[iloc][indy]=reshbondsWC[iloc][indy]+1
                    listresWC[indx,indy]=False
                    listresWC[indy,indx]=False
            else:
                tothbondsNWC[iloc][indx]=tothbondsNWC[iloc][indx]+1
                tothbondsNWC[iloc][indy]=tothbondsNWC[iloc][indy]+1
                if(listresNWC[indx,indy]==True and listresNWC[indy,indx]==True):
                    #           print(listres[indx,indy])
                    reshbondsNWC[iloc][indx]=reshbondsNWC[iloc][indx]+1
                    reshbondsNWC[iloc][indy]=reshbondsNWC[iloc][indy]+1
                    listresNWC[indx,indy]=False
                    listresNWC[indy,indx]=False
            tothbonds[iloc][indx]=tothbonds[iloc][indx]+1
            tothbonds[iloc][indy]=tothbonds[iloc][indy]+1
            listrhbonds[iloc][indy]=1
            listrhbonds[iloc][indx]=1
            #        print(i,indx,tothbonds[i][indx],indy,traj.topology.atom(hbonds[i][j][0]),traj.topology.atom(hbonds[i][j][2]))
            if(listres[indx,indy]==True and listres[indy,indx]==True):
                #           print(listres[indx,indy])
                reshbonds[iloc][indx]=reshbonds[iloc][indx]+1
                reshbonds[iloc][indy]=reshbonds[iloc][indy]+1
                listres[indx,indy]=False
                listres[indy,indx]=False
                #          print(listres[indx,indy])
            

aver=np.zeros((nres,2))
aver1=np.zeros((nres,2))
aver2=np.zeros((nres,2))
averWC=np.zeros((nres,2))
averNWC=np.zeros((nres,2))
averWCr=np.zeros((nres,2))
averNWCr=np.zeros((nres,2))

file1a=open('hbonds_tot_aver_traj.dat','w')
file1b=open('hbonds_tot_WC_aver_traj.dat','w')
file1c=open('hbonds_tot_NWC_aver_traj.dat','w')

file2a=open('hbonds_res_aver_traj.dat','w')
file2b=open('hbonds_res_WC_aver_traj.dat','w')
file2c=open('hbonds_res_NWC_aver_traj.dat','w')
#filetest=open('test','w')

for i in range(0,nres):
    #print(listrhbonds[0:ntframes+1,i:i+1])

    anumb='%5d' % (i+1)
    
    
    aver[i,0]=np.mean(listrhbonds[0:ntframes,i:i+1])
    aver[i,1]=np.std(listrhbonds[0:ntframes,i:i+1])


    
    aver2[i,0]=np.mean(tothbonds[0:ntframes,i:i+1])
    aver2[i,1]=np.std(tothbonds[0:ntframes,i:i+1])

    aaver='%8.3f' % (aver2[i,0])
    afluc='%8.3f' % (aver2[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    if(i==0):
        #print(tothbonds[0:ntframes+1,i:i+1])
        np.savetxt('test.dat',tothbonds[0:ntframes,i:i+1],fmt='%8.2f')
    file1a.write(stringa+'\n')
    
    averWC[i,0]=np.mean(tothbondsWC[0:ntframes,i:i+1])
    averWC[i,1]=np.std(tothbondsWC[0:ntframes,i:i+1])

    aaver='%8.3f' % (averWC[i,0])
    afluc='%8.3f' % (averWC[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    file1b.write(stringa+'\n')
    
    averNWC[i,0]=np.mean(tothbondsNWC[0:ntframes,i:i+1])
    averNWC[i,1]=np.std(tothbondsNWC[0:ntframes,i:i+1])

    aaver='%8.3f' % (averNWC[i,0])
    afluc='%8.3f' % (averNWC[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    file1c.write(stringa+'\n')
    
    aver1[i,0]=np.mean(reshbonds[0:ntframes,i:i+1])
    aver1[i,1]=np.std(reshbonds[0:ntframes,i:i+1])
    aaver='%8.3f' % (aver1[i,0])
    afluc='%8.3f' % (aver1[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    file2a.write(stringa+'\n')

    averWCr[i,0]=np.mean(reshbondsWC[0:ntframes,i:i+1])
    averWCr[i,1]=np.std(reshbondsWC[0:ntframes,i:i+1])    

    aaver='%8.3f' % (averWCr[i,0])
    afluc='%8.3f' % (averWCr[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    file2b.write(stringa+'\n')
    
    
    averNWCr[i,0]=np.mean(reshbondsNWC[0:ntframes,i:i+1])
    averNWCr[i,1]=np.std(reshbondsNWC[0:ntframes,i:i+1])
    
    aaver='%8.3f' % (averNWCr[i,0])
    afluc='%8.3f' % (averNWCr[i,1])
    stringa=anumb+' '+Listname[i]+'  '+aaver+' '+afluc
    file2c.write(stringa+'\n')



np.savetxt('totWC_traj.dat',tothbondsWC[0:ntframes,0:nres],fmt='%5d')
np.savetxt('resWC_traj.dat',reshbondsWC[0:ntframes,0:nres],fmt='%5d')  
    
#np.savetxt('hbonds_formation_traj.dat',listrhbonds[0:ntframes+1,0:nres],fmt='%5d')
#np.savetxt('hbonds_formation_aver_traj.dat',aver,fmt='%5.2f')

#np.savetxt('hbonds_res_aver_traj.dat',aver1,fmt='%5.2f')
#np.savetxt('hbonds_res_WC_aver_traj.dat',averWCr,fmt='%5.2f')
#np.savetxt('hbonds_res_NWC_aver_traj.dat',averNWCr,fmt='%5.2f')

#np.savetxt('hbonds_tot_aver_traj.dat',aver2,fmt='%5.2f')
#np.savetxt('hbonds_tot_WC_aver_traj.dat',averWC,fmt='%5.2f')
#np.savetxt('hbonds_tot_NWC_aver_traj.dat',averNWC,fmt='%5.2f')


