
# produceSmartPixMuC

Setup muon collider software environment
```
source setup.sh
```


To produce a BIB dataset for input to PixelAV
```
python study_bib.py
```

To produce a signal dataset for input to PixelAV
* First produce a muon gun sample, in lcio format
* Then run GEANT4 simulation
* Then produce text file for PixelAV
```
source particle_gun.sh
source detector_sim.sh
python study_signal.py
```
