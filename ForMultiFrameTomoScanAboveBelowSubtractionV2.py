import h5py
import cv2
import os
from os import listdir
from os.path import isfile, join
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import ion
from scipy.signal.signaltools import wiener
import EllipseForKaz
from scipy import math
from skimage.restoration import unwrap_phase



        
      
def myRec(obj,continueLoop,pathTot,dataFolder):  
    ### recursive function to look for the data database
    temp=None
    i=1
    tempPath=''
    for name, value in obj.items():
        if continueLoop:
            #check if the object is a group
            if isinstance(obj[name], h5py.Group):
                tempPath='/'+name
                if len(obj[name])>0:
                    continueLoop,temp,tempPath= myRec(obj[name],continueLoop,tempPath,dataFolder)
                else:
                    continue
            else:
                test=obj[name]
                temp1='/'+dataFolder
                if temp1 in test.name:
                    continueLoop=False
                    tempPath=pathTot+'/'+name
                    return continueLoop,test.name,tempPath
            i=i+1
        if (i-1)>len(obj.items()):
            tempPath=''
    pathTot=pathTot+tempPath
    return continueLoop,temp, pathTot

def findBestFlat(imgRef,pippo):
    p1=0
    p2=10
    goodp2=0
    goodStd=10
    for i in range(-30,30):
        p2=i
        print 'Im looking for the best flat'
        M = np.float32([[1,0,-p1],[0,1,-p2]])  
        rows,cols = imgRef.shape  
        print rows,cols, 'rows'
        test = cv2.warpAffine(imgRef.astype(np.float64),M,(cols,rows))
        #plt.figure(100)
        #plt.imshow(test[rows-10:rows,:])
        #plt.show()
        bho = np.ones([abs(p2),cols])# * imgRef[rows-p2-1,:]
        #test[rows-p2:rows,:]=bho
        print np.std(np.std(imgRef)), 'std1'
        test=test/avg1*avg2
        print np.std(np.std(test,0)), 'std2'
        
        if p2>0:
            bho = np.ones([abs(p2),cols]) * imgRef[rows-1,:]
            test[rows-p2:rows,:]=bho    
        else:
            bho = np.ones([abs(p2),cols]) * imgRef[0,:]
            test[0:abs(p2),:]=bho
        
        
        #plt.figure(1)
        
        #plt.imshow(test,'Greys')
        #plt.figure(2)
        #plt.imshow(imgRef,'Greys')
        #plt.figure(3)
        uffa=pippo/test
        print np.std(np.std(uffa,0)), 'std3'
        #plt.figure(111)
        #plt.imshow(uffa,'Greys',clim=(0.8,1.2))
        #plt.show()
        myStd=np.std(np.std(uffa,0))
        if myStd<goodStd:
            goodStd=myStd
            goodp2=p2
            print 'new std2',goodp2, goodStd
    bestP2=goodp2      
    mySmallerStep=0.1        
    
    return goodp2

def flatImages(pathAbove,pathFlatAbove,pathFlatAbove,pathSave):
        
    '''
    LOADING ABOVE
    '''
    mypathAbove=h5py.File(pathAbove,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotAbove=''
    contLoop, pathToData, pathTotAbove=myRec(mypathAbove,contLoop,pathTotAbove,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images',np.shape(mypathAbove[str(pathTotAbove)])
    a,b,c=np.shape(mypathAbove[str(pathTotAbove)])
    
    mypathFlatAbove=h5py.File(pathFlatAbove,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotFlatAbove=''
    contLoop, pathToData, pathTotFlatAbove=myRec(mypathFlatAbove,contLoop,pathTotFlatAbove,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images flat above',np.shape(mypathFlatAbove[str(pathTotFlatAbove)])
    #a,b,c=np.shape(mypathFlatAbove[str(pathTotFlatAbove)])
    
    mypathDarkAbove=h5py.File(pathFlatAbove,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotDarkAbove=''
    contLoop, pathToData, pathTotDarkAbove=myRec(mypathDarkAbove,contLoop,pathTotDarkAbove,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images dark above',np.shape(mypathDarkAbove[str(pathTotDarkAbove)])
    
    
    '''
    LOADING BELOW
    '''
    mypathBelow=h5py.File(pathBelow,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotBelow=''
    contLoop, pathToData, pathTotBelow=myRec(mypathBelow,contLoop,pathTotBelow,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images',np.shape(mypathBelow[str(pathTotBelow)])
    a2,b2,c2=np.shape(mypathBelow[str(pathTotBelow)])
    
    mypathFlatBelow=h5py.File(pathFlatBelow,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotFlatBelow=''
    contLoop, pathToData, pathTotFlatBelow=myRec(mypathFlatBelow,contLoop,pathTotFlatBelow,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images flat below',np.shape(mypathFlatBelow[str(pathTotFlatBelow)])
    #a,b,c=np.shape(mypathFlatAbove[str(pathTotFlatAbove)])
    
    mypathDarkBelow=h5py.File(pathDarkBelow,'r') 
    folder='data'
    print 'looking for "',folder, '" in the tree...'
    contLoop=True
    pathTotDarkBelow=''
    contLoop, pathToData, pathTotDarkBelow=myRec(mypathDarkBelow,contLoop,pathTotDarkBelow,folder)
    #npdata=np.array(mypath[str(pathTot)])
    print 'images dark below',np.shape(mypathDarkBelow[str(pathTotDarkBelow)])
    
    #raw_input('wait for any key')
    
    flatAbove=mypathFlatAbove[str(pathTotFlatAbove)][0,:,:]
    darkAbove=mypathDarkAbove[str(pathTotDarkAbove)][0,:,:]
    flatBelow=mypathFlatBelow[str(pathTotFlatBelow)][0,:,:]
    darkBelow=mypathDarkBelow[str(pathTotDarkBelow)][0,:,:]
    
    
    imagesToAnalise=a2
    key=np.zeros(a2)
    name=pathSave
    vortexIm=h5py.File(name,"w")
    gr=vortexIm.create_group('/entry1/tomo_entry/data/')
    dsetImage=vortexIm.create_dataset('/entry1/tomo_entry/data/data', (a2,b2,c2), 'f')
    dsetAngle=vortexIm.create_dataset('/entry1/tomo_entry/data/rotation_angle', (a2,), 'f')
    #dsetAngle2=merlinIm.create_dataset('/entry1/tomo_entry/data/rotation_angle2', (a,), 'f')
    dsetKey=vortexIm.create_dataset('/entry1/tomo_entry/instrument/detector/image_key', (a2,), 'f')
    dsetKey[...]=key 
    #dsetImage=vortexIm.create_dataset('data', (imagesToAnalise,b2,c2), 'f')
    #dsetImage[...]=totImage3
    
    for i in range(imagesToAnalise):
        '''
        ABOVE
        '''
        imAbove=(mypathAbove[str(pathTotAbove)][i,:,:])-darkAbove#/(flatAbove-darkAbove)
        
        flAbove=flatAbove-darkAbove
        
        imgRef=flAbove[1200:1750,200:700]
        avg1=np.average(np.average(imgRef))
        print 'average1',avg1
        pippo=imAbove[1200:1750,200:700]
        
        avg2=np.average(np.average(pippo))
        print 'average2',avg2
        #ret2 = cv2.phaseCorrelate(pippo.astype(np.float64),imgRef.astype(np.float64))
        
        #p1,p2=ret2

        shift=findBestFlat(imgRef,pippo)
        print shift, 'best shift'
        p1=0
        M = np.float32([[1,0,-p1],[0,1,-shift]])  
        rows,cols = flAbove.shape  
        print rows,cols, 'rows'
        flltBest = cv2.warpAffine(flAbove.astype(np.float64),M,(cols,rows))    
        print np.shape(np.ones([abs(shift),cols]))
        #print np.shape(flAbove[rows-shift-1,:])
        bho = np.ones([abs(shift),cols])
        if shift>0:
            bho = np.ones([abs(shift),cols]) * flAbove[rows-1,:]
            flltBest[rows-abs(shift):rows,:]=bho  
        else:
            bho = np.ones([abs(shift),cols]) * flAbove[0,:]
            flltBest[0:abs(shift),:]=bho  
        uffa=(imAbove)/flltBest
        print np.std(np.std(uffa,0)), 'std3'
        #plt.figure(200)
        #plt.imshow(uffa,'Greys',clim=(0,0.2))
        #plt.figure(300)
        #plt.imshow(imAbove/flAbove,'Greys',clim=(0,0.2))
        #plt.show()
        
        
        '''
        BELOW
        '''
        
        
        imBelow=(mypathBelow[str(pathTotBelow)][i,:,:])-darkBelow#/(flatAbove-darkAbove)
        
        flBelow=flatBelow-darkBelow
        
        imgRefBelow=flBelow[1200:1750,200:700]
        avg1=np.average(np.average(imgRefBelow))
        print 'average1',avg1
        pippoBelow=imBelow[1200:1750,200:700]
        
        avg2=np.average(np.average(pippoBelow))
        print 'average2',avg2
        #ret2 = cv2.phaseCorrelate(pippo.astype(np.float64),imgRef.astype(np.float64))
        
        #p1,p2=ret2

        shiftBelow=findBestFlat(imgRefBelow,pippoBelow)
        print shiftBelow, 'best shift'
        p1=0
        M = np.float32([[1,0,-p1],[0,1,-shiftBelow]])  
        rows,cols = flBelow.shape  
        print rows,cols, 'rows'
        flltBestBelow = cv2.warpAffine(flBelow.astype(np.float64),M,(cols,rows))    
        print np.shape(np.ones([abs(shiftBelow),cols]))
        #print np.shape(flAbove[rows-shift-1,:])
        bhoBelow = np.ones([abs(shiftBelow),cols])
        if shiftBelow>0:
            bhoBelow = np.ones([abs(shiftBelow),cols]) * flBelow[rows-1,:]
            flltBestBelow[rows-abs(shiftBelow):rows,:]=bhoBelow  
        else:
            bhoBelow = np.ones([abs(shiftBelow),cols]) * flBelow[0,:]
            flltBestBelow[0:abs(shiftBelow),:]=bhoBelow  
        uffaBelow=(imBelow)/flltBestBelow
        print np.std(np.std(uffaBelow,0)), 'std3'
        #plt.figure(400)
        #plt.imshow(uffaBelow,'Greys',clim=(0,0.2))
        #plt.figure(500)
        #plt.imshow(flBelow,'Greys',clim=(0,0.2))
        #plt.show()
        
        #plt.figure(600)
        #plt.imshow(uffa-uffaBelow,'Greys',clim=(0,0.05))
        #plt.show()
        
        dsetImage[i,:,:]= uffa-uffaBelow
    
   
#########For testing function
if __name__ == "__main__":
    
    '''
    width=2961
    height=1128
    #binning=10
    
    cropping parameters
    
    minRow=322
    minCol=567
    depthProjections=2#1493
    startFormThisLeg=0
    nFull=depthProjections
    numOffset=6
    legs=[0,2,4,1,3,5]
    key=np.zeros(depthProjections)
    angles=np.zeros(depthProjections)
    '''
    
    pathAbove="/dls/i13/data/2017/mt15805-1/processing/91461Processed.hdf"
    pathFlatAbove="/dls/i13/data/2017/mt15805-1/processing/91461Flat.hdf"
    pathDarkAbove="/dls/i13/data/2017/mt15805-1/processing/91461Dark.hdf"
    pathBelow="/dls/i13/data/2017/mt15805-1/processing/91462Processed.hdf"
    pathFlatBelow="/dls/i13/data/2017/mt15805-1/processing/91462Flat.hdf"
    pathDarkBelow="/dls/i13/data/2017/mt15805-1/processing/91462Dark.hdf"

        
        
        
        
        
    vortexIm.close()
    print 'all done!!!'    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    