from utils import getTracks

plot = False

# Set up some options, constants
max_events = -1 # Set to -1 to run over all events
#Bfield = 3.57 # T for legacy

npart = 0
nevts = 0

# gather input files 
directory_path = "/cvmfs/public-uc.osgstorage.org/ospool/uc-shared/public/futurecolliders/BIB10TeV/sim_mm_pruned"
file_paths = glob.glob(f"{directory_path}/*.slcio")

tracks=getTracks(file_paths, allowedPIDS=[11,13,22], plot=False, max_events=max_events, flp=0, tracklist_folder="~/MuonColliderSim/Tracklists/BIB_tracklists", binsize=500, float_precision=5)