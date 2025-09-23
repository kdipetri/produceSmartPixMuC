import numpy as np
import pyLCIO
import ROOT
import glob
import os
import json
from math import *
import csv

from plothelper import *

# setup plotter
plt = PlotHelper()

sensorAngles = np.arange(-np.pi,np.pi+2*np.pi/8,np.pi/8)

def getYlocalAndGamma(x,y):
    # Get exact angle gamma of the hit position
    gammaP=np.arctan2(y,x)

    # Get two sensor angles that are closest to the exact angle
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
    
    ylocal=-round(yentry/25e-3)*25e-3
    # at some point, add limits to possible ROIs

    if index==0:
        index=16
    if index==17:
        index=1
    
    gamma=sensorAngles[index]

    return ylocal, gamma

def make_plots(mcp_tlv, hit_tlv, hit, prodrxy, prodz, endrxy, endz, t):
    plt.plot1D("hit_mcp_e"  ,";mcp e [GeV];hits" , mcp_tlv.E(), 100, 0, 0.2)
    plt.plot1D("hit_mcp_pt"  ,";mcp pt [GeV];hits" , mcp_tlv.Pt(), 100, 0, 0.2)
    plt.plot1D("hit_mcp_eta" ,";mcp eta;hits" , mcp_tlv.Eta(), 100, -3.2, 3.2)
    plt.plot1D("hit_mcp_theta" ,";mcp theta;hits" , mcp_tlv.Theta(), 100, 0, 3.2)
    plt.plot1D("hit_mcp_phi" ,";mcp phi;hits" , mcp_tlv.Phi(), 100, -3.2, 3.2)
    plt.plot1D("hit_pt"  ,";incident pt [GeV];hits" , hit_tlv.Pt(), 100, 0, 0.2)
    plt.plot1D("hit_eta" ,";incident eta;hits" , hit_tlv.Eta(), 100, -3.2, 3.2)
    plt.plot1D("hit_theta" ,";incident theta;hits" , hit_tlv.Theta(), 100, 0,3.2)
    plt.plot1D("hit_phi" ,";incident phi;hits" , hit_tlv.Phi(), 100, -3.2, 3.2)
    plt.plot1D("hit_eDep" ,";incident e deposit [MeV];hits" , hit.getEDep()*1000, 100, 0, 0.5)
    plt.plot1D("hit_mcp_prodrxy" ,";mcp prod rxy [mm];hits" , prodrxy, 100, 0,150)
    plt.plot1D("hit_mcp_prodz"   ,";mcp prod z [mm];hits" , prodz, 100, -1000,1000)
    plt.plot1D("hit_mcp_endrxy"  ,";mcp end rxy [mm];hits" , endrxy, 100, 0, 150)
    plt.plot1D("hit_mcp_endz"    ,";mcp end z [mm];hits" , endz, 100, -1000, 1000)

    plt.plot1D("hit_time"    ,";hit time [ns];hits" , t, 100, -1, 5)
    print("phi particle,hit {:.2f} {:.2f}".format(hit_tlv.Phi(), mcp_tlv.Phi()))
    print("eta particle,hit {:.2f} {:.2f}".format(hit_tlv.Eta(), mcp_tlv.Eta()))

    # double check if any bugs
    phi = hit_tlv.Phi()
    theta = hit_tlv.Theta()

    plt.plot1D("hit_phi"    ,";cota;hits" , phi, 100, -10,10)
    plt.plot1D("hit_theta"    ,";cotb;hits" , theta, 100, -10,10)
    plt.plot1D("hit_t"    ,";t;hits" , t, 100, -1,10)

def getTracks(file_paths, allowedPIDS=None, plot=False, max_events=-1, flp=0, tracklist_folder=None, binsize=500, float_precision=5, overwrite=False):
    
    # ############## SETUP #############################
    # Prevent ROOT from drawing while you're running -- good for slow remote servers
    ROOT.gROOT.SetBatch()

    if allowedPIDS is None:
        raise Exception("Must provide allowedPIDS to get tracks")
    if file_paths is None or len(file_paths)==0:
        raise Exception("Must provide non-empty file_paths to get tracks")
    
    # check if output directory exists
    if not os.path.isdir(tracklist_folder):
        os.makedirs(tracklist_folder)
    #
    #Check this first, so it doens't waste time reading files if the directory is already in use
    #
    elif overwrite:
        # Empty folder
        files = os.listdir(tracklist_folder)
        for f in files:
            if os.path.isfile(f"{tracklist_folder}/{f}"):
                os.remove(f"{tracklist_folder}/{f}")
    else:
        raise Exception(f"Directory {tracklist_folder} already exists. To overwrite, set overwrite=True")
    
    # convert this to writing to a log file
    print(f"Getting tracks from files: {file_paths}\n")
    if max_events != -1:
        print("Reading all events.")
    else:
        print(f"Limiting to max_events = {max_events}\n")
    print(f"Allowed PIDs: {allowedPIDS}\n")

    tracks=[]
    nevts=0
    for file_path in file_paths:
        if not file_path.endswith(".slcio"):
            raise Exception("Input file must be a .slcio file")

        print("Processing file: ", file_path)

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

            for hit in vtxBarrelHits: 
                position = hit.getPosition()
                
                x,y,hit_z = position[0],position[1],position[2]
                #rxy=(x**2+y**2)**0.5
                t=hit.getTime()

                # Get layer
                encoding = vtxBarrelHits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                decoder = pyLCIO.UTIL.BitField64(encoding)
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                #detector = decoder["system"].value()
                layer = decoder['layer'].value()
                #side = decoder["side"].value()

                if layer!=0: continue # first layer only

                # get the particle that caused the hit
                mcp = hit.getMCParticle()
                hit_pdg = mcp.getPDG() if mcp else None
                hit_id = mcp.id() if mcp else None

                if hit_id is None:
                    print("No MCParticle associated with hit, skipping.")
                    continue

                if abs(hit_pdg) not in allowedPIDS: continue

                # momentum at production
                mcp_p = mcp.getMomentum()
                mcp_tlv = ROOT.TLorentzVector()
                mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())

                # momentum at hit
                hit_p = hit.getMomentum()
                hit_tlv = ROOT.TLorentzVector()
                hit_tlv.SetPxPyPzE( hit_p[0], hit_p[1], hit_p[2], mcp.getEnergy())

                prodx,prody,prodz=mcp.getVertex()[0],mcp.getVertex()[1],mcp.getVertex()[2]
                endx,endy,endz=mcp.getEndpoint()[0],mcp.getEndpoint()[1],mcp.getEndpoint()[2]
                prodrxy = (prodx**2 + prody**2)**0.5
                endrxy = (endx**2 + endy**2)**0.5

                if plot: make_plots(mcp_tlv, hit_tlv, hit, prodrxy, prodz, endrxy, endz, t)

                ylocal, gamma0 = getYlocalAndGamma(x,y)
                zglobal = round(hit_z/25e-3)*25e-3 # round to nearest pixel
                
                # Define unit vector of track at tracker edge with respect to barrel
                theta=hit_tlv.Theta()
                phi=hit_tlv.Phi()
                vx_global=np.sin(theta)*np.cos(phi)
                vy_global=np.sin(theta)*np.sin(phi)
                vz_global=np.cos(theta) 

                # Convert vector to sensor frame to calculate alpha
                vy_local=vx_global*np.sin(np.pi/2-gamma0)+vy_global*np.cos(np.pi/2-gamma0)
                vx_local=vz_global
                alpha=np.arctan2(vy_local,vx_local)
                
                # Calculate beta from phi and gamma0
                beta=phi-(gamma0-np.pi/2)

                # If we are unflipped, we must adjust alpha and beta, and flip y-local
                if flp == 0:
                    beta += np.pi
                    alpha = 2*np.pi-alpha
                    ylocal *= -1

                cota = 1./np.tan(alpha)
                cotb = 1./np.tan(beta)

                p = mcp_tlv.P()
                pt = mcp_tlv.Pt()
                track = [cota, cotb, p, flp, ylocal, zglobal, pt, t, hit_pdg]
                tracks.append(track)

                #print("")
                #print("NEW PARTICLE")

                # helpful printout
                #print("x,y,z,t={:.1f},{:.1f},{:.1f},{:.3f}".format(x,y,z,t))
                ##print("  det,lay,side={},{},{}".format(detector,layer,side))
                #print("  prod rxy,z={:.1f},{:.1f}".format(prodrxy,prodz))
                #print("  end  rxy,z={:.1f},{:.1f}".format(endrxy ,endz))
                #print("  pt,theta,phi,e={:.3f},{:.3f},{:.1f},{:.3f}".format(mcp_tlv.Pt(), mcp_tlv.Theta(), mcp_tlv.Phi(),mcp_tlv.E()))
                #print("  pt,theta,phi,e={:.3f},{:.3f},{:.1f},{:.3f}".format(hit_tlv.Pt(), hit_tlv.Theta(), hit_tlv.Phi(),hit_tlv.E()))
                #print("  pdg={}".format(hit_pdg))
                #print("  id={}".format(hit_id))

        writeTracklists(tracks, tracklist_folder, binsize, float_precision)
            

def writeTracklists(tracks, tracklist_folder=None, binsize=500, float_precision=5, overwrite=False):
    if tracklist_folder is None:
        raise Exception("Must provide tracklist_folder to writeTracklists")
    
    print(f"Writing {len(tracks)} tracks to {tracklist_folder} with {binsize} tracks per file\n")# check if directory exists
    
    if "sig" in tracklist_folder:
        tag = "sig"
    else:
        tag = "BIB"

    numFiles = int(np.ceil(len(tracks)/binsize))
    for fileNum in range(numFiles):
        with open(f"{tracklist_folder}/{tag}_tracklist{fileNum}.txt", 'w') as file:
            for track in tracks[fileNum*binsize:(fileNum+1)*binsize]:

                # set flp to an int
                track = list(track)

                formatted_sublist = [f"{element:.{float_precision}f}" if isinstance(element, float) else element for element in track]
                line = ' '.join(map(str, formatted_sublist)) + '\n'
                file.write(line)