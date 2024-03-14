

# Load the core tglow functions
for (file in list.files(paste0("/softare/teamtrynka/tglow-core/", "tglow-r-core"), full.names=T)) {
  source(file)
}


# Global parameters
suffix     <- "v1"
prefix     <- "tglow_panel_"
path       <- paste0("../../2_feature_extraction/output/", prefix, suffix)


# Sheet with sample metadata in standardized format
sample.sheet <- "https://docs.google.com/spreadsheets/d/1gdv1ZN9kQm7UHSqZ9urset7VoljR3r4xxt_2xYlaSqA/edit#gid=0"

# Sheet with cell and feature filters following standardized format
filter.sheet <- "https://docs.google.com/spreadsheets/d/1dpgoigdCFhnE8zPEuqe-M-L1OIz2S_vUVzaZkDiKU7A/edit#gid=1278275764"