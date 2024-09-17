# Check the geometry
echo $MUCOLL_GEO

# remove prev file
rm output_sim.slcio

# 10 events runs in 1-2 minutes
benchmark_dir=/home/karri/mucLLPs/mucoll-benchmarks
ddsim --steeringFile $benchmark_dir/simulation/ilcsoft/steer_baseline.py \
    --inputFile output_gen.slcio \
    --outputFile output_sim.slcio 
ls

