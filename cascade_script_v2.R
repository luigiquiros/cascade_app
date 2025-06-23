# ===========================================================================================================================================================================================================
# Title:            CASCADE: Contextualizing untargeted Annotation with Semi-quantitative Charged Aerosol Detection for pertinent characterization of natural Extracts.
# Author:           Luis Quiros, Ph.D.
# Organization:     Phytochimie et Produits Naturels Bioactifs, ISPSO, UNIGE, SWITZELAND
# Contact:          luis.guerrero@unige.ch
# Date Created:     03/12/2024
# Version           1.0          
# License:          CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
# Description:      Analyzes UHPLC-HRMS/MS data with semiquantitative information.
# ===========================================================================================================================================================================================================

#===================================== REQUIREMENTS ===========================================================================================================================================================

#   A file (.mzML) containing DDA MS data with an additional detector (PDA, ELSD, CAD)
#   In case you don’t know how to obtain it, see: wiki/How-to-create-a-compliant-mzML-file
#   A file (.csv) containing features, as obtained by mzmine comprehensive export
#   A file (.tsv) containing annotations, as obtained by TIMA
# ===========================================================================================================================================================================================================

#===================================== INSTALLATION ===========================================================================================================================================================

#install.packages(
#  "cascade",
#  repos = c(
#    "https://adafede.r-universe.dev",
#    "https://bioc.r-universe.dev",
#    "https://cloud.r-project.org"
#  )
#)
# ===========================================================================================================================================================================================================

#===================================== CASCADE STARTS HERE===================================================================================================================================================
#load packagea
 library(cascade)
 library(dplyr)   

# Check documentation for the package
#  ?cascade

#check all the functions available
# help(package = "cascade") #shows all the functions available in the package

# =================================== PATHS ==================================================================================================================================================================

# Define the main folder of your project and file names
 
 data_path <- "C:/Users/quirosgu/Documents/GitHub/cascade/data"              #"path/to/your/data"  # Define the path to your data
 filename_pos <- "calendula_pos.mzML"          #your indexed mzml in pos mode 
 filename_neg <- "calendula_neg.mzML"          #your indexed mzml in neg mode 
 feature_table_pos <- "calendula_pos_quant_full.csv"                #your MZmine features table in POS mod               #your MZmine features table in pos mode 
 feature_table_neg <- "calendula_neg_quant_full.csv"                #your MZmine features table in neg mode
 annotations_results_pos <- "calendula_pos_tima_results.tsv"      #your TimaR annotatios results table in pos mode
 annotations_results_neg <- "calendula_neg_tima_results.tsv"      #your TimaR annotatios results table in neg mode
 
# ===========================================================================================================================================================================================================
 # Do not modify this section 
 data_in_path       <- file.path(data_path, "in/")
 file_negative      <- file.path(data_in_path, filename_neg)
 file_positive      <- file.path(data_in_path, filename_pos)
 features_pos       <- file.path(data_in_path, feature_table_pos)
 features_neg       <- file.path(data_in_path, feature_table_neg)
 annotations_pos    <- file.path(data_in_path, annotations_results_pos)
 annotations_neg    <- file.path(data_in_path, annotations_results_neg)
 
 #create automatically a path to export data
 export_path <- file.path(data_path, "results")  # Save results inside the 'cascade' folder
 if (!dir.exists(export_path)) {
   dir.create(export_path)  # Create the directory if it doesn't exist
 }
# ================================= Validate the file paths =================================================================================================================================================
  if (!file.exists(file_negative)) {
   stop("Negative mzML file does not exist in the specified data_path.")
 }
 if (!file.exists(file_positive)) {
   stop("Positive mzML file does not exist in the specified data_path.")
 }
 if (!file.exists(features_pos)) {
   stop("Features table in pos file does not exist in the specified data_path.")
 }
 if (!file.exists(features_neg)) {
   stop("Features table in neg file does not exist in the specified data_path.")
 }
 if (!file.exists(annotations_pos)) {
   stop("Annotation results in pos file does not exist in the specified data_path.")
 }
 if (!file.exists(annotations_neg)) {
   stop("Annotation results in neg file does not exist in the specified data_path.")
 }
# ===========================================================================================================================================================================================================
 
# ================================= VARIABLES ===============================================================================================================================================================
 
# 0.1. define the shifts between the CAD and PDA traces in reference to the MS trace
 cad_shift <- -0.13
 pda_shift <- -0.04
 ms_shift  <- 0.0

 # 0.1. define the rt window you want to consider
 time_min <- 0.1
 time_max <- 25

# ===========================================================================================================================================================================================================

# ================================= ALIGNMENT ===============================================================================================================================================================
 
# 1. Apply the check_chromatograms_alignment function to verify the data
 
 result_plot <- check_chromatograms_alignment(
   file_negative = file_negative,
   file_positive = file_positive,
   time_min = time_min,             # Adjust as needed, chromatographic run time
   time_max = time_max,           # Adjust as needed, chromatographic run time
   cad_shift = cad_shift,           # CAD shift
   pda_shift = pda_shift,           # PDA shift
   fourier_components = 0.01,  # Fourier components
   frequency = 1,              # frequency
   resample = 1,               # Resampling factor
   chromatograms = c("bpi_pos", "cad_pos", "pda_pos", "bpi_neg"),  # Chromatograms to plot
   type = "baselined",         # Choose "baselined" or "improved"
   normalize_intensity = TRUE, # Normalize intensity?
   normalize_time = FALSE      # Normalize time?
 )

 # The function will return a plot
 print("Alignment check completed. Resulting plot:")
 print(result_plot)

 #NOTE: modify the variable accordingly to achieve the best alignment according to your experiment data
# ===========================================================================================================================================================================================================
 
# ================================= INTEGRATIONS ============================================================================================================================================================
 
 # 2. Apply the checK_peaks_integrations
 
 #2.1. cad trace
 
 result_peaks_cad <- check_peaks_integration(
   file = file_positive,        # Positive mzML file
   features = features_pos,     # Features table
   detector = "cad",            # Specify detector (e.g., "cad, bpi or pda")
   chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms
   min_area = 0.005,            # Minimum area threshold for integration
   min_intensity = 1e5,         # Minimum intensity threshold
   shift = cad_shift,                # Shift (CAD shift value)
   show_example = FALSE,        # Use your data, not example data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,              # Minimum retention time
   time_max = time_max,            # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1                 # Resampling factor
 )
 # Display the result
 result_peaks_cad
 
 #2.2. pda trace

 result_peaks_pda <- check_peaks_integration(
   file = file_positive,        # Positive mzML file
   features = features_pos,     # Features table
   detector = "pda",            # Specify detector (e.g., "cad")
   chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms
   min_area = 0.005,            # Minimum area threshold for integration
   min_intensity = 1e5,         # Minimum intensity threshold
   shift = pda_shift,           # Shift (CAD shift value)
   show_example = FALSE,        # Use your data, not example data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,         # Minimum retention time
   time_max = time_max,         # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1                 # Resampling factor
 )
 # Display the result
 result_peaks_pda
 
 #2.3. pos trace
 
 result_peaks_pos <- check_peaks_integration(
   file = file_positive,        # Positive mzML file
   features = features_pos,     # Features table
   detector = "bpi",            # Specify detector (e.g., "cad")
   chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms
   min_area = 0.005,            # Minimum area threshold for integration
   min_intensity = 1e5,         # Minimum intensity threshold
   shift = ms_shift,            # Shift (CAD shift value)
   show_example = FALSE,        # Use your data, not example data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,         # Minimum retention time
   time_max = time_max,         # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1                 # Resampling factor
 )
 # Display the result
 result_peaks_pos
 
 #2.4. neg trace
 
 result_peaks_neg <- check_peaks_integration(
   file = file_negative,        # Positive mzML file
   features = features_neg,     # Features table
   detector = "bpi",            # Specify detector (e.g., "cad")
   chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms
   min_area = 0.005,            # Minimum area threshold for integration
   min_intensity = 1e4,         # Minimum intensity threshold
   shift = ms_shift,            # Shift (CAD shift value)
   show_example = FALSE,        # Use your data, not example data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,         # Minimum retention time
   time_max = time_max,         # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1                 # Resampling factor
 )
 # Display the result
 result_peaks_neg
# ===========================================================================================================================================================================================================

# ================================= MS & CAD LINKS CONSTRUCTION =============================================================================================================================================
 
 # 3. Construct links between CAD, MS and annotations: Determination of major and minor metabolites
 
 # 3.1. compare peaks cad and pos
 
 compare_peaks_cad_pos <- process_compare_peaks(
   file = file_positive,        # Positive mzML file
   features = features_pos,     # Features table
   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
   detector = "cad",            # Specify detector
   type = "baselined",          # Chromatogram type
   min_area = 0.005,            # Minimum area threshold
   min_intensity = 1e5,         # Minimum intensity threshold
   shift = cad_shift,                # CAD shift value
   show_example = FALSE,        # Use your data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,              # Minimum retention time
   time_max = time_max,            # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1,                # Resampling factor
   export_dir = export_path     # Export directory
 )
 
 #Check if the files were exported
 list.files(export_path, full.names = TRUE)

 # 3.2. compare peaks cad and neg
 
 compare_peaks_cad_neg <- process_compare_peaks(
   file = file_negative,        # Positive mzML file
   features = features_neg,     # Features table
   detector = "cad",            # Specify detector
   type = "baselined",          # Chromatogram type
   min_area = 0.005,            # Minimum area threshold
   min_intensity = 1e4,         # Minimum intensity threshold
   shift = cad_shift,                # CAD shift value
   show_example = FALSE,        # Use your data
   fourier_components = 0.01,   # Fourier components for smoothing
   time_min = time_min,              # Minimum retention time
   time_max = time_max,            # Maximum retention time
   frequency = 1,               # Frequency
   resample = 1,                # Resampling factor
   export_dir = export_path     # Export directory
 )
 
 #Check if the files were exported
 list.files(export_path, full.names = TRUE)
# ===========================================================================================================================================================================================================
 
# ================================= PSEUDO CHROMATOGRAMS CONSTRUCTION =======================================================================================================================================
 
 # 4. Construct Pseudo chromatograms
 
# ===========================================================================================================================================================================================================
 # Do not modify this section 
 
 # Define the base name by removing the extension from filename_pos adn filename_neg
 base_name_pos <- sub("\\.mzML$", "", filename_pos)  # Remove '.mzML'
 base_name_neg <- sub("\\.mzML$", "", filename_neg)  # Remove '.mzML'
 # Construct the paths for the informed and not informed features files
 pos_features_informed_path <- file.path(export_path, paste0(base_name_pos, "_featuresInformed_cad.tsv"))
 pos_features_not_informed_path <- file.path(export_path, paste0(base_name_pos, "_featuresNotInformed_cad.tsv"))
 neg_features_informed_path <- file.path(export_path, paste0(base_name_neg, "_featuresInformed_cad.tsv"))
 neg_features_not_informed_path <- file.path(export_path, paste0(base_name_neg, "_featuresNotInformed_cad.tsv"))
# ===========================================================================================================================================================================================================
 
 #4.1. Run generate_pseudochromatograms in pos mode
 
 #4.1.1. combine the data and create the plots
 
 pseudochromatograms_result_pos <- generate_pseudochromatograms(
   annotations = annotations_pos,                # Annotations file
   features_informed = pos_features_informed_path,    # Features informed file
   features_not_informed = pos_features_not_informed_path, # Features not informed file
   file = file_positive,                              # mzML file
   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
   detector = "cad",                              # Detector type
   show_example = FALSE,                          # Use your data, not example data
   min_confidence = 0.4, #0.4                          # Minimum confidence for annotations
   min_similarity_prefilter = 0.6, #0.6               # Pre-filter similarity threshold
   min_similarity_filter = 0.8,  #0.8                 # Filter similarity threshold
   mode = "pos",                                  # Ionization mode ("pos" or "neg")
   organism = "Calendula officinalis",                # Organism name
   fourier_components = 0.01,                     # Fourier components for smoothing
   frequency = 2,                                 # Frequency for resampling
   resample = 2,                                  # Resampling factor
   shift = cad_shift,                                  # Detector shift
   time_min = time_min,                                # Minimum retention time
   time_max = time_max                                # Maximum retention time
 )
 
 histograms_pos_taxo_maj <- pseudochromatograms_result_pos$plots_1$histograms_taxo_maj
 histograms_pos_taxo_min <- pseudochromatograms_result_pos$plots_1$histograms_taxo_min
 histograms_pos_unique_conf_maj <- pseudochromatograms_result_pos$plots_1$histograms_unique_conf_maj
 histograms_pos_unique_conf_min <- pseudochromatograms_result_pos$plots_1$histograms_unique_conf_min
 treemap_pos_peaks_maj <- pseudochromatograms_result_pos$treemaps$peaks_maj
 treemap_pos_peaks_min <- pseudochromatograms_result_pos$treemaps$peaks_min
 treemap_pos_special <- pseudochromatograms_result_pos$treemaps$special
 
 #4.1.2. access the plots by name display the plots
 
 #Taxo (major peaks)
 histograms_pos_taxo_maj  
 
 #Taxo (minor peaks)
 histograms_pos_taxo_min
 
 #Confident unique annotations (major)
 histograms_pos_unique_conf_maj
 
 #Confident unique annotations (major)
 histograms_pos_unique_conf_min
 
 #Treemap semi-quantitative (major)
 treemap_pos_peaks_maj
 
 #Treemap semi-quantitative (minor)
 treemap_pos_peaks_min
 
 #Treemap both (major and minor)
 treemap_pos_special
 
 
 #4.3. Run generate_pseudochromatograms in neg mode
 
 pseudochromatograms_result_neg <- generate_pseudochromatograms(
   annotations = annotations_neg,                # Annotations file
   features_informed = neg_features_informed_path,    # Features informed file
   features_not_informed = neg_features_not_informed_path, # Features not informed file
   file = file_negative,                              # mzML file
   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
   detector = "cad",                              # Detector type
   show_example = FALSE,                          # Use your data, not example data
   min_confidence = 0.4, #0.4                          # Minimum confidence for annotations
   min_similarity_prefilter = 0.6, #0.6               # Pre-filter similarity threshold
   min_similarity_filter = 0.8,  #0.8                 # Filter similarity threshold
   mode = "pos",                                  # Ionization mode ("pos" or "neg")
   organism = "Calendula officinalis",                # Organism name
   fourier_components = 0.01,                     # Fourier components for smoothing
   frequency = 1,                                 # Frequency for resampling
   resample = 1,                                  # Resampling factor
   shift = cad_shift,                                  # Detector shift
   time_min = time_min,                                # Minimum retention time
   time_max = time_max                                # Maximum retention time
 )
 
 histograms_neg_taxo_maj <- pseudochromatograms_result_neg$plots_1$histograms_taxo_maj
 histograms_neg_taxo_min <- pseudochromatograms_result_neg$plots_1$histograms_taxo_min
 histograms_neg_unique_conf_maj <- pseudochromatograms_result_neg$plots_1$histograms_unique_conf_maj
 histograms_neg_unique_conf_min <- pseudochromatograms_result_neg$plots_1$histograms_unique_conf_min
 treemap_neg_peaks_maj <- pseudochromatograms_result_neg$treemaps$peaks_maj
 treemap_neg_peaks_min <- pseudochromatograms_result_neg$treemaps$peaks_min
 treemap_neg_special <- pseudochromatograms_result_neg$treemaps$special
 
 #4.1.2. access the plots by name display the plots
 
 #Taxo (major peaks)
 histograms_neg_taxo_maj  
 
 #Taxo (minor peaks)
 histograms_neg_taxo_min
 
 #Confident unique annotations (major)
 histograms_neg_unique_conf_maj
 
 #Confident unique annotations (major)
 histograms_neg_unique_conf_min
 
 #Treemap semi-quantitative (major)
 treemap_neg_peaks_maj
 
 #Treemap semi-quantitative (minor)
 treemap_neg_peaks_min
 
 #Treemap both (major and minor)
 treemap_neg_special
 
# ===========================================================================================================================================================================================================
 
# ================================= GENERAL TABLE RESULTS GENERATION ========================================================================================================================================
 
 #5.  Create the general results tables
 
 #5.1.  Run generate_tables in pos
 
 table_list_pos<- generate_tables(
   annotations = annotations_pos,    # Path to the annotations file
   file_negative = neg_features_informed_path,             # Provide negative file if available
   file_positive = pos_features_informed_path,    # Positive processed peaks file
   min_confidence = 0.4,             # Minimum confidence level
   show_example = FALSE,             # Use your data, not example data
   export_csv = TRUE,                # Export as CSV
   export_html = TRUE,               # Export as HTML
   export_dir = export_path,         # Export directory
   export_name = "cascade_table_results_pos"        # Export file name
 )
 
 # Access the 'pretty_table' element
 pretty_table_pos <- table_list_pos$pretty_table #interactive table
 pretty_table_pos <- table_list_pos$csv_table #csv version
 
 #5.2.  Run generate_tables in neg
 
 table_list_neg <- generate_tables(
   annotations = annotations_neg,    # Path to the annotations file
   file_negative = NULL,             # Provide negative file if available
   file_positive =  neg_features_informed_path,    # Positive processed peaks file
   min_confidence = 0.4,             # Minimum confidence level
   show_example = FALSE,             # Use your data, not example data
   export_csv = TRUE,                # Export as CSV
   export_html = TRUE,               # Export as HTML
   export_dir = export_path,         # Export directory
   export_name = "cascade_table_results_neg"        # Export file name
 )
 
 # Access the 'pretty_table' element
 pretty_table_neg <- table_list_neg$pretty_table #interactive table
 pretty_table_neg <- table_list_neg$csv_table #csv version
 
 # treat pos and neg together 
 
 table_list_d_pos <- generate_tables(
   annotations = annotations_pos,    # Path to the annotations file
   file_negative = neg_features_informed_path,             # Provide negative file if available
   file_positive = pos_features_informed_path,    # Positive processed peaks file
   min_confidence = 0.4,             # Minimum confidence level
   show_example = FALSE,             # Use your data, not example data
   export_csv = TRUE,                # Export as CSV
   export_html = TRUE,               # Export as HTML
   export_dir = export_path,         # Export directory
   export_name = "cascade_table_full_results_neg"        # Export file name
 )
 
 # Access the 'pretty_table' element
 pretty_table_neg <- table_list_neg$pretty_table #interactive table
 pretty_table_neg <- table_list_neg$csv_table #csv version
 
 
# if lotus doesn't work use: 
 # Manually download the LOTUS file and place it in the directory
 #download.file(
 #  url = "https://doi.org/10.5281/zenodo.5794106",
 #  destfile = "data/source/libraries/230106_frozen_metadata.csv.gz"
 # )
 
# ===========================================================================================================================================================================================================
# ================================= THAT'S ALL FOLKS ========================================================================================================================================================
# ===========================================================================================================================================================================================================

 # ================================ BONUS STAGE =============================================================================================================================================================
 #Bonus.  Quick overview of the literature for the species of interest
 
 bonus_plots_list <- generate_ids(
   taxa = c("Passiflora", "Calendula offici", "Kopsia", "Ginkgo"),
   comparison = c("Swertia", "Kopsia"),
   no_stereo = TRUE,
   filter_ms_conditions = TRUE,
   start = "1950",
   end = "2025"
 )
 #Compounds found in Genus
 bonus_plots_list$plots$Calendula

 #Compounds found in Genus (per species)
 bonus_plots_list$plots$Calendula_grouped
# ===========================================================================================================================================================================================================

 
