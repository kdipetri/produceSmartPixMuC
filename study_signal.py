import pyLCIO
import ROOT
import glob
import os
import json
from math import *
import numpy as np
import csv

from plothelper import *

# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
ROOT.gROOT.SetBatch()

# setup plotter
plt = PlotHelper()

plot = False

# Set up some options, constants
max_events = 1000 # Set to -1 to run over all events
max_npart = 100000 
Bfield = 3.57 # T for legacy

npart = 0
nevts = 0

sensorAngles = np.arange(-np.pi,np.pi+2*np.pi/8,np.pi/8)

def getYlocalAndGamma(x,y):
    # Get gamma of the point
    gammaP=np.arctan2(y,x)

    # Get two sensor angles that are closest to gammaP

    diff = np.abs(sensorAngles-gammaP)
    index1 = np.argmin(diff)
    gamma1=sensorAngles[index1]

    diff[index1]=3*np.pi
    index2 = np.argmin(diff)
    gamma2=sensorAngles[index2]

    # Rotate x coordinate of the point by each option for gamma
    x1=x*np.cos(-gamma1)-y*np.sin(-gamma1)
    y1=y*np.cos(-gamma1)+x*np.sin(-gamma1)
    x2=x*np.cos(-gamma2)-y*np.sin(-gamma2)
    y2=y*np.cos(-gamma2)+x*np.sin(-gamma2)

    # Determine which x is closest to expected value
    xTrue=30.16475324197002

    diff1=abs(x1-xTrue)
    diff2=abs(x2-xTrue)
    
    # If both x1 and x2 are really close to the ex
    if diff1 < 0.5 and diff2 < 0.5:
        if y1>8.5 or y1<-4.5:
            index=index2
        else:
            index=index1
            
    elif diff1<diff2:
        index=index1
    else:
        index=index2

    if index==index1:
        yentry=y1
    else:
        yentry=y2
    
    ylocal=yentry+12.5*13/2/1000

    if index==0:
        index=16
    if index==17:
        index=1
    
    gamma=sensorAngles[index]

    return ylocal, gamma

# setup ouptut data
tracks = [] # [["cota", "cotb", "p", "flp", "localx", "localy", "pT", "hittime", "PID"]]

# gather input files 
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
directory_path = "./"
for filename in os.listdir(directory_path):

    if "sim.slcio" not in filename: continue
    # Get the full path to the file
    file_path = os.path.join(directory_path, filename)
    print(file_path)
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file_path)
    for ievt,event in enumerate(reader):
        nevts+=1
        if max_events!=-1 and nevts > max_events: break
        print("Processing event %i."%ievt)

        # Print all the collection names in the event
        collection_names = event.getCollectionNames()

        # Get vertex barrel hits
        vtxBarrelHits = event.getCollection("VertexBarrelCollection")

        hit_id_last = -1
        for hit in vtxBarrelHits:
            position = hit.getPosition()
            x,y,z = position[0],position[1],position[2]
            rxy=(x**2+y**2)**0.5
            t=hit.getTime()

            # Get layer
            encoding = vtxBarrelHits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
            decoder = pyLCIO.UTIL.BitField64(encoding)
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            detector = decoder["system"].value()
            layer = decoder['layer'].value()
            side = decoder["side"].value()

            if layer!=0: continue # first layer only

            # get the particle that caused the hit
            mcp = hit.getMCParticle()
            hit_pdg = mcp.getPDG() if mcp else None
            hit_id = mcp.id() if mcp else None

            if abs(hit_pdg) != 11 and abs(hit_pdg) != 211 and abs(hit_pdg) != 13: continue

            # momentum at production
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())

            #momentum at hit
            hit_p = hit.getMomentum()
            hit_tlv = ROOT.TLorentzVector()
            hit_tlv.SetPxPyPzE( hit_p[0], hit_p[1], hit_p[2], mcp.getEnergy())

            prodx,prody,prodz=mcp.getVertex()[0],mcp.getVertex()[1],mcp.getVertex()[2]
            endx,endy,endz=mcp.getEndpoint()[0],mcp.getEndpoint()[1],mcp.getEndpoint()[2]
            prodrxy = (prodx**2 + prody**2)**0.5
            endrxy = (endx**2 + endy**2)**0.5
            
            if plot:
                plt.plot1D("sig_hit_mcp_p"  ,";mcp p [GeV];hits" , mcp_tlv.P(), 100, 0, 100)
                plt.plot1D("sig_hit_mcp_pt"  ,";mcp pt [GeV];hits" , mcp_tlv.Pt(), 100, 0, 100)
                plt.plot1D("sig_hit_mcp_eta" ,";mcp eta;hits" , mcp_tlv.Eta(), 100, -3.2, 3.2)
                plt.plot1D("sig_hit_mcp_theta" ,";mcp theta;hits" , mcp_tlv.Theta(), 100, 0, 3.2)
                plt.plot1D("sig_hit_mcp_phi" ,";mcp phi;hits" , mcp_tlv.Phi(), 100, -3.2, 3.2)
                plt.plot1D("sig_hit_p"  ,";incident p [GeV];hits" , hit_tlv.P(), 100, 0, 100)
                plt.plot1D("sig_hit_pt"  ,";incident pt [GeV];hits" , hit_tlv.Pt(), 100, 0, 100)
                plt.plot1D("sig_hit_eta" ,";incident eta;hits" , hit_tlv.Eta(), 100, -3.2, 3.2)
                plt.plot1D("sig_hit_theta" ,";incident theta;hits" , hit_tlv.Theta(), 100, 0,3.2)
                plt.plot1D("sig_hit_phi" ,";incident phi;hits" , hit_tlv.Phi(), 100, -3.2, 3.2)
                plt.plot1D("sig_hit_mcp_prodrxy" ,";mcp prod rxy [mm];hits" , prodrxy, 100, 0,150)
                plt.plot1D("sig_hit_mcp_prodz"   ,";mcp prod z [mm];hits" , prodz, 100, -1000,1000)
                plt.plot1D("sig_hit_mcp_endrxy"  ,";mcp end rxy [mm];hits" , endrxy, 100, 0, 150)
                plt.plot1D("sig_hit_mcp_endz"    ,";mcp end z [mm];hits" , endz, 100, -1000, 1000)

                phi = 1./np.tan(hit_tlv.Phi()) 
                theta = 1./np.tan(hit_tlv.Theta())

                plt.plot1D("sig_hit_phi"    ,";cota;hits" , phi, 100, -15, 15)
                plt.plot1D("sig_hit_theta"    ,";cotb;hits" , theta, 100, -15, 15)
                plt.plot1D("sig_hit_t"    ,";hit time [ns];hits" , t, 100, -1, 10)

            ylocal, gamma0 = getYlocalAndGamma(x,y)
            xlocal=0
            
            # Define unit vector of track at tracker edge with respect to barrel
            theta=hit_tlv.Theta()
            phi=hit_tlv.Phi()
            x=np.sin(theta)*np.cos(phi)
            y=np.sin(theta)*np.sin(phi)
            z=np.cos(theta) 

            # Transform into rotated coordinate system (sensor coordinate system sort of)
            xp=x*np.cos(np.pi/2-gamma0)+y*np.sin(np.pi/2-gamma0)
            yp=-x*np.sin(np.pi/2-gamma0)+y*np.cos(np.pi/2-gamma0)
            
            beta=np.arctan2(yp,xp)
            alpha=np.arctan2(yp,z)

            # Since we are unflipped, we must adjust alpha and beta 
            cotb = 1./np.tan(beta+np.pi)
            cota = 1./np.tan(2*np.pi-alpha)

            p = mcp_tlv.P()
            pt = mcp_tlv.Pt()
            track = [cota, cotb, p, 0, xlocal, ylocal, pt, hit_pdg]
            tracks.append(track)

            print("")
            print("NEW PARTICLE")

            # helpful printout
            print("x,y,z,t={:.1f},{:.1f},{:.1f},{:.3f}".format(x,y,z,t))
            ##print("  det,lay,side={},{},{}".format(detector,layer,side))
            print("  hit e={:.2e}".format(hit.getEDep()) )
            print("  hit dE/dx={:.2e}".format(hit.getdEdx()) )
            print("  prod rxy,x,y,z={:.1f},{:.1f},{:.1f},{:.1f}".format(prodrxy,prodx,prody,prodz))
            print("  end  rxy,x,y,z={:.1f},{:.1f},{:.1f},{:.1f}".format(endrxy, endx, endy, endz))
            print("  pt,theta,phi,e={:.4f},{:.3f},{:.1f},{:.4f}".format(mcp_tlv.Pt(), mcp_tlv.Theta(), mcp_tlv.Phi(),mcp_tlv.E()))
            print("  pt,theta,phi,e={:.4f},{:.3f},{:.1f},{:.4f}".format(hit_tlv.Pt(), hit_tlv.Theta(), hit_tlv.Phi(),hit_tlv.E()))
            print("  pdg={}".format(hit_pdg))
            #print("  id={}".format(hit_id))

            hit_id_last = hit_id

        # if looking at particles themselves 
        #mcpCollection = event.getCollection("MCParticle")
        
        #print("n mcp: ", len(mcpCollection))
        #for mcp in mcpCollection:
        #    npart+=1
        #    if max_npart!=-1 and npart > max_npart : break

        #    # only look at charged particles
        #    if mcp.getCharge() == 0 : continue

        #    mcp_p = mcp.getMomentum()
        #    mcp_tlv = ROOT.TLorentzVector()
        #    mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
        #    pdg = mcp.getPDG()
        #    print("pt,eta,phi,pdg={:.3f},{:.1f},{:.1f},{}".format(mcp_tlv.Pt(), mcp_tlv.Eta(), mcp_tlv.Phi(), pdg))

        #    status = mcp.getGeneratorStatus()
        #    prodx,prody,prodz=mcp.getVertex()[0],mcp.getVertex()[1],mcp.getVertex()[2]
        #    endx,endy,endz=mcp.getEndpoint()[0],mcp.getEndpoint()[1],mcp.getEndpoint()[2]
        #    prodrxy = (prodx**2 + prody**2)**0.5
        #    endrxy = (endx**2 + endy**2)**0.5
        #    print("prod rxy,z={:.1f},{:.1f}".format(prodrxy,prodz))
        #    print("end  rxy,z={:.1f},{:.1f}".format(endrxy ,endz))

        ## only one event
    # only one file

# save histos to file
fout = ROOT.TFile.Open("plots/sig_out.root","RECREATE")
fout.cd()
plt.drawAll()

# Writing to csv file
filename = "sig_tracklist.txt"
float_precision=5
with open(filename, 'w') as file:
    for track in tracks:

        # set flp to an int
        track = list(track)
        #track[3] = int(track[3])

        formatted_sublist = [f"{element:.{float_precision}f}" if isinstance(element, float) else element for element in track]
        line = ' '.join(map(str, formatted_sublist)) + '\n'
        file.write(line)
