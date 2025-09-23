import pyLCIO
import ROOT
import glob
import os
import json
from math import *
import numpy as np
import csv

#from study_signal import sensorAngles, getYlocalAndGamma

from plothelper import *

# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
ROOT.gROOT.SetBatch()

# setup plotter
plt = PlotHelper()

plot = False

# Set up some options, constants
max_events = -1 # Set to -1 to run over all events
max_npart = 100000 
Bfield = 3.57 # T for legacy

npart = 0
nevts = 0

# setup ouptut data
tracks = [['phi', 'theta', 'x', 'y', 'z', 'pt', 'prodx', 'prody', 'prodz', 't', 'pid']] 

# gather input files 
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
print("running")
directory_path = "/cvmfs/muoncollider.cern.ch/datasets/bib/MuColl_v1/example/"
for filename in os.listdir(directory_path):

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

            if abs(hit_pdg) != 11 and abs(hit_pdg) != 13: continue

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
            
            # Define unit vector of track at tracker edge with respect to barrel
            theta=hit_tlv.Theta()
            phi=hit_tlv.Phi()
            pt = mcp_tlv.Pt()

            track = [phi, theta, x, y, z, pt, prodx, prody, prodz, t, hit_pdg]
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
        #break
    # only one file
    #break
print(nevts)
# save histos to file 
"""
fout = ROOT.TFile.Open("plots/out.root","RECREATE")
fout.cd()
plt.drawAll()
"""
# Writing to csv file
filename = "BIB_hits.txt"
float_precision=5
with open(filename, 'w') as file:
    for track in tracks:

        # set flp to an int
        track = list(track)
        #track[3] = int(track[3])

        formatted_sublist = [f"{element:.{float_precision}f}" if isinstance(element, float) else element for element in track]
        line = ' '.join(map(str, formatted_sublist)) + '\n'
        file.write(line)

