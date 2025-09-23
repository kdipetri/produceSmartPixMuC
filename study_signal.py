from utils import getTracks

plot = False

# Set up some options, constants
max_events = -1 # Set to -1 to run over all events
#Bfield = 3.57 # T for legacy

file_paths = "./output_sim.slcio"
getTracks(file_paths, allowedPIDS=[13], plot=plot, max_events=max_events, flp=0, tracklist_folder="~/MuonColliderSim/Tracklists/signal_tracklists", binsize=500, float_precision=5, overwrite=True)