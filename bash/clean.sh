# Clean folder of all output

# Intermediate data
rm data/intermediate/*.csv
rm data/intermediate/*.Rdata
rm data/table-construction-data/*.xlsx

# Log files
rm log/*.log
rm log/basic-analysis/*.log
rm log-bash/*.o*
rm log-bash/*.e*

# MATLAB functions
# rm matlab/function/*t_distribution.m
# .mat files. Right now keep it so that we don't have to redo the MLE
# rm matlab/mat/*.mat

# Output
rm output/constants/*.txt
rm output/figures/*.eps
rm output/figures/*.pdf
rm output/tables/*.tex

