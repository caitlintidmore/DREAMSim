import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import ROOT
import pickle
from scipy import interpolate

def getG4data(dataset,t,df):
   df["event"]={"dataset":dataset,"calibSS":100./t.calibSph,"calibCC":100./t.calibCph, \
   "beamID":t.beamID,"beamType":t.beamType,"beamE":t.beamE,"beamX":t.beamX, \
   "beamY":t.beamY,"beamY":t.beamZ
   }

   id3dCC=np.array(t.id3dCC,dtype=np.int32)  # xxxyyyttt, xx0, yyY, Y=0 for 6mm, 1,2,3,4 for 3mm
   ph3dCC=np.array(t.ph3dCC,dtype=np.double)
   df["g4hit3dCC"]={"nhits":len(ph3dCC),"id":id3dCC,"val":ph3dCC,"sum":np.sum(ph3dCC)}

   id3dSS=np.array(t.id3dCC,dtype=np.int32)  # xxxyyyzzz
   ph3dSS=np.array(t.ph3dCC,dtype=np.double)
   df["g4hit3dSS"]={"nhits":len(ph3dSS),"id":id3dSS,"val":ph3dSS,"sum":np.sum(ph3dSS)}

   # create 2D data from 3d data...
   val={}
   for j,ph in enumerate(ph3dCC):
      kk=math.floor(id3dCC[j]/1000)
      val[kk]=val.get(kk,0.0)+ph

   ph2=np.array(list(val.values()),dtype=float)
   id2=np.array(list(val.keys()),dtype=int)

   df["g4hit2dCC"]={"nhits":len(ph2),"id":id2,"val":ph2,"sum":np.sum(ph2)}
   #
   val={}
   for j,ph in enumerate(ph3dSS):
      kk=math.floor(id3dSS[j]/1000)
      val[kk]=val.get(kk,0.0)+ph

   ph2=np.array(list(val.values()),dtype=float)
   id2=np.array(list(val.keys()),dtype=int)

   df["g4hit2dSS"]={"nhits":len(ph2),"id":id2,"val":ph2,"sum":np.sum(ph2)}
   
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

def makeDRSwave_v01(df,n_splwave,splwave,rng):
   ph3d=df["g4hit3dCC"]["val"]
   #  accumulate hits in a time frame...
   drshit={}
   h1={}
   xmin=0.0
   xmax=6.0    # limit for histrogram
   xmax2=10.0  # limit for getting spline wave
   nbins=30
   deltaT=0.0001
   for k, ph in enumerate(ph3d):
      ixy=math.floor(df["g4hit3dCC"]["id"][k]/1000)
      it =math.floor(df["g4hit3dCC"]["id"][k]%1000)
      s=str(ixy)
      if s not in h1:
         h1[s]=ROOT.TH1D(s,s,nbins,xmin,xmax)
      x=it*0.05-6.2   # the grid in G4 is 50ps step along z-axis
      y=ph
      xx,yy=pick_one_splinewave(x,xmax2,n_splwave,splwave,rng)
      wt=yy*y
      h1[s].FillN(int(len(xx)),xx+deltaT,wt)

   df0={}
   for key in h1.keys():
      vals=np.zeros(nbins,dtype=float)
      for i in np.arange(1,nbins+1,1):
          ii=int(i)
          vals[ii-1]=h1[key].GetBinContent(ii)
      df0[int(key)]=vals*(5.0/50.)  # 5ps steps in spline wave and 50 ps from G4
      h1[key].Delete()

   return df0

def pick_one_splinewave(xval,xmax,n_splwave,splwave,rng):
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

def plotPulseTrain_v01(df,chList):
   #  loop over channel list...
   for ch in chList:
      if ch in df["drsCC"].keys():
         print("(def plotPulseTrain_v01) ch",ch)
         a=(((df["g4hit3dCC"]["id"])/1000).astype(int))==ch
         # print("type(a)",type(a)," a=",a)
         g4ph=df["g4hit3dCC"]["val"][a]
         g4x=((df["g4hit3dCC"]["id"][a]%1000)-133)*0.05
         plt.hist(g4x,bins=120,weights=g4ph)
         # print("g4ph",g4ph)
         # print("g4x",g4x)

         drs=df["drsCC"][ch]*(5.0/100.0)
         xdrs=np.arange(0.0,len(drs),1.0)*0.2   # 200 ps sampling)
         plt.plot(xdrs,drs)
         phsum=np.sum(g4ph)
         s=" {}  ch {:5d}  g4 {:6.2f} photons".format(df["event"]["dataset"],ch,phsum)
         plt.title(s)
         plt.xlabel("ns") 
         plt.ylabel("g4(photons),drs(a.u.)")
         plt.grid(True)
         plt.show()

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
