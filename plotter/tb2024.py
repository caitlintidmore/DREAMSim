import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import ROOT
import pickle
from scipy import interpolate

def getG4data(name,t,df):
        df["event"]={"calibSS":100./t.calibSph,"calibCC":100./t.calibCph}

        ix3dCC=np.array(t.ix3dCC,dtype=np.int32)
        iy3dCC=np.array(t.iy3dCC,dtype=np.int32)
        ts3dCC=np.array(t.it3dCC,dtype=np.int32)
        ph3dCC=np.array(t.ph3dCC,dtype=np.double)
        chid=ix3dCC*1000000+iy3dCC*1000+ts3dCC
        df["simhitC3d"]={"nhits":len(ph3dCC),"chid":chid, "val":ph3dCC}

        ix2dCC=np.array(t.ix2dCC,dtype=np.int32)
        iy2dCC=np.array(t.iy2dCC,dtype=np.int32)
        ts2dCC=np.array(t.it2dCC,dtype=np.int32)
        ph2dCC=np.array(t.ph2dCC,dtype=np.double)
        chid=ix2dCC*1000000+iy2dCC*1000+ts2dCC
        df["simhitC2d"]={"nhits":len(ph2dCC),"chid":chid, "val":ph2dCC}

        ix3dSS=np.array(t.ix3dSS,dtype=np.int32)
        iy3dSS=np.array(t.iy3dSS,dtype=np.int32)
        ts3dSS=np.array(t.it3dSS,dtype=np.int32)
        ph3dSS=np.array(t.ph3dSS,dtype=np.double)
        chid=ix3dSS*1000000+iy3dSS*1000+ts3dSS
        df["simhitS3d"]={"nhits":len(ph3dSS),"chid":chid, "val":ph3dSS}

        ix2dSS=np.array(t.ix2dSS,dtype=np.int32)
        iy2dSS=np.array(t.iy2dSS,dtype=np.int32)
        ts2dSS=np.array(t.it2dSS,dtype=np.int32)
        ph2dSS=np.array(t.ph2dSS,dtype=np.double)
        chid=ix2dSS*1000000+iy2dSS*1000+ts2dSS
        df["simhitS2d"]={"nhits":len(ph2dSS),"chid":chid, "val":ph2dSS}

        return

def getPulseShapeFromPkl(fname):
   with open('pulsedata_1033.pkl','rb') as fp: 
       df = pickle.load(fp)
       # print("df...") print("type(df",type(df))
       print("(def getPulseShapeFromPkl)  fname= ",fname)
       print("Using run: ",df["run"])
       print("number of pulses (events): ",df["events"])
       print("chlist:",df["chlist"])
   return df

def get_splwave(dfPulseShape):
    print("(def get_splwave) creating a set of single photon pulse shape with ")
    print("spline fits to measured pulse shapes in the pulse shape library from pkl file.")
    n_splwave=0
    splwave={}

    for key in dfPulseShape:
        #  note: length of tuple is 2 and those of other keys are then umbe of charecters
        #  tuple:   [0] for event id numberm and [1] for channel number
        if key[1]==9:
            ts=dfPulseShape[key]["ts"]
            x=np.arange(len(ts))*0.2   # 200 ps per timeslice
            tck=interpolate.splrep(x,ts)
            xnew=np.arange(0.0,40.0,0.005)  # 8000 steps, 5ps/step
            ynew=interpolate.splev(xnew,tck,der=0) 
            splwave[n_splwave]=ynew
            n_splwave=n_splwave+1
            if n_splwave==1:  #plot the first one for visual validation
                print("Example of a pulse shape from spline fit.")
                print("40 ns wide,  5 ps step.  Total 8000 steps.")
                ynew2=ynew[ynew>0.0]
                ysum2=np.sum(ynew2)
                ysum=np.sum(ynew)
                ysum3=np.sum(ynew[ynew<0.0])
                print("ysum=",ysum,"  ysum2=",ysum2," ysum3 ",ysum3,"  ysum2+ysum3",ysum2+ysum3)
                plt.plot(x,ts+0.1,label="200ps/step, DRS data, 0.1 shift in y")
                plt.plot(xnew,ynew,label="5ps/step spline curve")
                # plt.xlim(40,60.0)    # full range is [0,200]
                plt.xlabel("ns")
                plt.legend(loc='upper right')
                plt.grid(True)
                plt.show()
    return n_splwave,splwave,xnew

def f_splinewave(xval,xmax,n_splwave,splwave,rng): 
    # pick a waveform randomly
    dt=0.005
    ir=int(rng.random()*n_splwave)
    x=np.arange(0.0,xmax*2.0,dt)  # x array is wider for scanning beyond xmax.
    nmax=int(xmax/dt)
    x1=int((10.0-xval)/dt)
    x2=x1+nmax
    xx=x[:nmax]
    yy=splwave[ir][x1:x2]
    return xx,yy

def makeNew2d(name,ix,iy,it,ph,calib):
    # this creates course 2D, ie. hits in 3mm SiPM area is merged to the area same as 6 mm SiPM.
    val = {}
    for j, ja in enumerate(ph):
        kx=ix[j]
        ky=iy[j]
        kk=kx*100+math.floor(ky/10)
        val[kk]=val.get(kk,0.0)+ph[j]


    ph2=np.array(list(val.values()),dtype=float)
    arg2=np.argsort(-ph2)
    id2=np.array(list(val.keys()),dtype=int)

    ph3=ph2[arg2]
    ix3=np.array(id2[arg2]/1000,dtype=int)*10
    iy3=np.array(id2[arg2]%100,dtype=int)*10
    '''
    print("  ")
    print("New:  new 2d...")
    for j, ja in enumerate(ph3):
        if j>10:
            break
        s=" 2d k {:3d}  (ix,iy)  ({:3d}, {:3d})  ph {:.2f}  en {:.2f} GeV ".format(j,ix3[j],iy3[j],ph3[j],ph3[j]*calib)
        print(s) 
    '''
    return ix3,iy3,ph3

def makeNew3d(name,ix,iy,it,ph,calib):
    val = {}
    for j, ja in enumerate(ph):
        kx=ix[j]
        ky=iy[j]
        kt=math.floor(it[j])
        kk=(kx*100+math.floor(ky/10))*1000+kt
        val[kk]=val.get(kk,0.0)+ph[j]

    ph2=np.array(list(val.values()),dtype=float)
    arg2=np.argsort(-ph2)
    id2=np.array(list(val.keys()),dtype=int)

    ph3=ph2[arg2]
    kkxy=np.array(id2[arg2]/1000,dtype=int)
    ix3=np.array(kkxy/1000,dtype=int)*10
    iy3=np.array(kkxy%100,dtype=int)*10
    it3=np.array(id2[arg2]%1000,dtype=int)
    '''
    print("  ")
    print("NEW:  new 3d...")
    for j, ja in enumerate(arg2):
        if j>10:
            break
        kkxy=math.floor(id2[ja]/1000)
        kx=math.floor(kkxy/1000)
        ky=math.floor(kkxy%100)
        kt=math.floor(id2[ja]%1000)
        s=" 3d k {:3d}  (ix,iy,it)  ({:3d}, {:3d} {:4d} )  ph {:.2f}  en {:.2f} GeV ".format(j,ix3[j],iy3[j],it3[j],ph3[j],ph3[j]*calib)

        print(s) 
    '''
    return ix3,iy3,it3,ph3


def makeDRSwaveTest(df,n_splwave,splwave,rng):
   ph=10.0
   x=4.0
   xmax=10.0
   xx,yy=f_splinewave(x,xmax,n_splwave,splwave,rng)
   print("len(xx)",len(xx)," len(yy)",len(yy))
   plt.plot(xx,yy)
   plt.grid(True)
   plt.show()   

def makeDRSwave(df,n_splwave,splwave,rng):
   ph3d=df["simhitC3d"]["val"]
   #  accumulate hits in a time frame...
   drshit={}
   h1={}
   xmin=0.0
   xmax=6.0    # limit for histrogram
   xmax2=10.0  # limit for getting spline wave
   nbins=30
   deltaT=0.0001  
   for k, ph in enumerate(ph3d):
      ixy=math.floor(df["simhitC3d"]["chid"][k]/1000)
      it =math.floor(df["simhitC3d"]["chid"][k]%1000)
      s=str(ixy) 
      if s not in h1:
         h1[s]=ROOT.TH1D(s,s,nbins,xmin,xmax)
      x=it*0.05-6.2   # the grid in G4 is 50ps step along z-axis
      y=ph
      xx,yy=f_splinewave(x,xmax2,n_splwave,splwave,rng)
      wt=yy*y
      h1[s].FillN(int(len(xx)),xx+deltaT,wt)

   df["drswaveC"]={}
   for key in h1.keys():
      vals=np.zeros(nbins,dtype=float)
      for i in np.arange(1,nbins+1,1):
          ii=int(i)
          vals[ii-1]=h1[key].GetBinContent(ii)
      df["drswaveC"][key]=vals*(5.0/50.)  # 5ps steps in spline wave and 50 ps from G4
      h1[key].Delete()
   
   # print("makeDRShits...")
   # print(dfevt["drshitC"])

def plotPulseShape(dfPulseShape): 
    print("(def plotPulseShape) type(df",type(dfPulseShape))
    print("run: ",dfPulseShape["run"])
    print("number of events: ",dfPulseShape["events"])
    print("chlist:",dfPulseShape["chlist"]) 
    for k in dfPulseShape["chlist"]:
        print("ch ",k,"  ",dfPulseShape["chnames"][k])
    #  loop over all rows...
    ch0=1
    for key in dfPulseShape:
        # note: length of tuple is 2 and those of other keys are then umbe of␣ ↪charecters
        #  tuple:   [0] for event id numberm and [1] for channel number
        if key[1]==ch0:
            # print("key (selected)=",key,"type(key)",type(key)) 
            ts=dfPulseShape[key]["ts"]*(-1.0)
            x=np.arange(len(ts))
            plt.plot(x,ts)
            #
    plt.ylim(-0.4,1.1)
    plt.title("(def plotPulseShape) Run "+str(dfPulseShape["run"])+": "+dfPulseShape["chnames"][ch0]) 
    plt.xlabel("DRS ts, 200ps")
    plt.ylabel("pulse hight (normalized to 1.0)")
    plt.grid(True)
    plt.show()

def test():
   print("print out in tb2024/test")
   return 1

def main():
   ok=test()
   return

if __name__ == "__main__":
    main()
