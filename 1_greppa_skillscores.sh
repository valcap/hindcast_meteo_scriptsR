#!/bin/bash

OBS_DATASET='ARCIS'

Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep CORR |  cut -d';' -f2-5 > ../res/tp_CORR_vs_${OBS_DATASET}_years.csv
Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep RMSE |  cut -d';' -f2-5 > ../res/tp_RMSE_vs_${OBS_DATASET}_years.csv
Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep MBIAS |  cut -d';' -f2-5 > ../res/tp_MBIAS_vs_${OBS_DATASET}_years.csv
Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep OBIAS |  cut -d';' -f2-5 > ../res/tp_OBIAS_vs_${OBS_DATASET}_years.csv
Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep ME |  cut -d';' -f2-5 > ../res/tp_ME_vs_${OBS_DATASET}_years.csv
Rscript 1_verification_tp_skillscor_${OBS_DATASET}.r | grep VER_ | grep MAE |  cut -d';' -f2-5 > ../res/tp_MAE_vs_${OBS_DATASET}_years.csv

exit 0;

