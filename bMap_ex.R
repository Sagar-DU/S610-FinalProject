# Packages must be installed and loaded into R
# remove.packages("BayesBrainMap")
install.packages("devtools")
devtools::install_github("mandymejia/fMRItools", "7.0")
devtools::install_github("mandymejia/BayesBrainMap", "2.0")
#install.packages("BayesBrainMap")
#install.packages("ciftiTools")
library(fMRItools)
library(ciftiTools)
library(BayesBrainMap)

# Set CIFTI Workbench path (download from https://www.humanconnectome.org/software/connectome-workbench)
wb_path <- "C:/Program Files/workbench" # path to where you downloaded it
ciftiTools.setOption("wb_path", wb_path)

# Main function is BrainMap
# Most important parameters are BOLD, prior, and TR

# BOLD which is the user's data for Human Connectome Project can be found in slate
# For each subject there are two paths you can use:
# /N/project/hcp_dcwan/<subject_id>/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii
# /N/project/hcp_dcwan/<subject_id>/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii

# For the prior you can download it from here: https://osf.io/k6vx8/overview?view_only=b614888e9aca42999ee75eb2c3e02877
# Example: priors/GICA15/prior_combined_GICA15_GSR.rds

# TR = 0.72 for HCP data

# Example function call:
BOLD <- "./sample-data/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii"
prior <- "./sample-data/prior/prior_combined_GICA15_GSR.rds"
bMap <- BrainMap(
  BOLD = BOLD, # path to BOLD data
  prior = prior, # path to prior
  TR = 0.72
)

# Can run engagements after with bMap
eng <- engagements(
        bMap
        )
