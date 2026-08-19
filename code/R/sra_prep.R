library(tidyverse)

# Set path to the gzipped files 
path <- "data/raw_gzipped_fastqs_symlinks/"

forward_reads <- sort(list.files(path, pattern = "R1.fastq.gz"))
forward_reads

reverse_reads <- sort(list.files(path, pattern = "R2.fastq.gz"))
reverse_reads

# Get clean names
sample_names <- data.frame(bname = forward_reads) %>% # Extract the basename (file name)
  separate_wider_delim(bname, delim = "_", names_sep = "_", too_few = "align_start") %>% # Split up by the underscores
  transmute(sample_name  = paste(bname_3, bname_4, sep = "_")) %>%
  pull()

sample_names

sra_data <- data.frame(sample_name = sample_names,
                       filename = forward_reads,
                       filename2 = reverse_reads) 

data.frame(sample_name = sample_names,
           filename = forward_reads,
           filename2 = reverse_reads) %>%
  write_csv("sra_files.csv")

library(fun.gus)

load_schmidt_metadata(path_to_data_dir = "../SCHMIDT_LAB_DATA_LOG/",
                      name_vector = c("Station_Log", "Deployment_Log", "Freezer_Log", "FP_Process_Details", "FC_Process_Details", "DNA_Log"),
                      project = "AAM",
                      as_list = FALSE)



sample_info_raw <- Freezer_Log %>%
  dplyr::filter(Location_ID == "ERI") %>%
  select(Sample_ID, Station_ID, Deployment_ID, Date_Collected, Notes) %>%
  left_join(Deployment_Log, by = "Deployment_ID") %>% 
  select(Sample_ID, Station_ID = Station_ID.x, Date_Collected, Deployment_Depth_m, Integrated_Depths_m, Temp_at_Deployment, Notes.x, Notes.y) %>%
  mutate(Depth = ifelse(!is.na(Deployment_Depth_m), Deployment_Depth_m, 
                        ifelse(!is.na(Integrated_Depths_m), 5, NA)),
         Notes = paste(Notes.x, Notes.y)) %>%
  select(-Deployment_Depth_m, -Integrated_Depths_m, -Notes.x, -Notes.y) %>% 
  left_join(Station_Log, by = c("Station_ID", "Date_Collected" = "Date")) %>%
  mutate(Notes = paste(Notes.x, Notes.y)) %>%
  rowwise() %>%
  mutate(Notes = list(c(Notes))) %>%
  ungroup() %>%
  select(Sample_ID, Station_ID, Date = Date_Collected, Depth, Station_Name, Water_Surface_Temp_C, Latitude, Longitude, Notes) %>%
  mutate(Month = case_when(str_detect(Date, "202405") ~ "May",
                           str_detect(Date, "202408") ~ "August",
                           str_detect(Date, "202409") ~ "September"),
         Depth_Class = ifelse(Depth == 5, "Integrated Surface", "Bottom"),
         Rep_ID = paste(Month, Station_ID, Depth_Class, sep = "_")) %>%
  left_join(read_tsv("data/05_Metadata_Curation/basin_info.tsv")) %>%
  left_join(DNA_Log, by = c("Sample_ID")) %>%
  filter(!is.na(DNA_ID)) %>% 
  mutate(sample_name = DNA_ID, 
         sample_title = Sample_ID,
         date = ymd(Date),
         Basin = case_when(Basin == "E" ~ "Eastern",
                           Basin == "C" ~ "Central",
                           Basin == "W" ~ "Western",
                           TRUE ~ NA)) %>%
  select(sample_name, sample_title, Station_ID, date, depth = Depth, EPA_Station_Name = Station_Name, Latitude, Longitude, Month, Depth_Class, Basin)

head(sample_info_raw)




bio_metadata <- sample_info_raw %>%
  mutate(collection_date = as.character(date),
         organism = "lake metagenome",
         env_broad_scale = "freshwater lake biome [ENVO:01000252]",
         env_local_scale = ifelse(Depth_Class == "Integrated_Surface", "epilimnion [ENVO:00002131]", "hypolimnion [ENVO:00002130]"),
         env_medium = "lake water [ENVO:04000007]",
         geo_loc_name = "USA: Lake Erie",
         lat_lon = paste(Latitude, "N", abs(as.numeric(Longitude)), "W"),
         size_frac_up = "0.22 micrometer") %>%
  select(-Latitude, -Longitude)

sra_data %>%
  full_join(bio_metadata) %>%
  write_csv("full_sra_data.csv")


