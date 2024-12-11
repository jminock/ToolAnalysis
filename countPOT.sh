#!/bin/bash
# Based on /annie/app/users/dingpf/GridSub_test.sh


# Create a dummy file in the output directory. This is a hack to get jobs
# that fail to end themselves quickly rather than hanging on for a long time
# waiting for output to arrive.
#DUMMY_OUTPUT_FILE=${CONDOR_DIR_OUTPUT}/${JOBSUBJOBID}_dummy_output
#touch ${DUMMY_OUTPUT_FILE}

# Copy datafile in ToolAnalysis

for RUN in {10..99}
do
    sed -i "s/gntp.0.ghep.root/gntp.${RUN}.ghep.root/g" configfiles/CC_MC_RECO_ntuple/LoadGenieEventConfig
    ./Analyse configfiles/CC_MC_RECO_ntuple/ToolChainConfig
    sed -i "s/gntp.${RUN}.ghep.root/gntp.0.ghep.root/g" configfiles/CC_MC_RECO_ntuple/LoadGenieEventConfig
done



# setup software

#echo "Moving the output files to CONDOR OUTPUT:" >> ${DUMMY_OUTPUT_FILE}

#echo "Cleaning up:" >> ${DUMMY_OUTPUT_FILE}
### END ###
