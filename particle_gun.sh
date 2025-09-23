# Runs in a few seconds

# It is useful for performance studies to generate events with single particles (muons, pions, electrons, etc) 
# with specific parameters. This is colloquially called the particle gun. 

# A dedicated Python script allows to generate an LCIO file with stable particles, represented by LCIO::MCParticle objects.
# For example, the following command will generate 10 events with 1 electron per event, where the electron has p = 100 GeV,
# polar angle randomly distributed in the range [10deg, 170deg], and azimuthal angle randomly distributed in the range [0, 2pi]:

# 45000 events should get around 80 tracklist files with 500 tracks per file   ??


rm output_gen.slcio
benchmark_dir=/home/karri/mucLLPs/mucoll-benchmarks
python $benchmark_dir/generation/pgun/pgun_lcio.py \
    -s 12345 \
    -e 50000 \
    --pdg 13 -13 \
    --p 1 100 \
    --theta 10 170 \
    --dz 0 0 1.5 \
    --d0 0 0 0.0009 \
    -- output_gen.slcio
ls