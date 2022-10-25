source("./GCMClimTool_function_library_CMIP6.R")
######### Modifiable data for particular application #####################################################

mFolder = "H:/Shared drives/GCM_Climate_Tool"
rFolder <- "C:/Users/dmf/projects/gcm-project/test-results"
shp_file_background <- "/required_files/OBSERVATIONS/LLdM_AOI2/SHP/Basin_LldM.shp" # Basin Shapefile
BootstrapRuns <- 1 # number of iterations to perform of the bootstrapping routine
CLIMATE_VARS = c("pr")
# setwd(mFolder)

### ------------------------------
### SITE SPECIFIC DATA SOURCES: climate observations
### ------------------------------
observation_files_lock = c(daily = "/required_files/OBSERVATIONS/LLdM_AOI2/series_pr_diario.csv",
                           monthly = "/required_files/OBSERVATIONS/LLdM_AOI2/series_pr_mensual.csv",
                           stationsCatalog = "/required_files/OBSERVATIONS/LLdM_AOI2/Catalogo.csv")


# in the "BoundingBoxClimaZones.csv". add rows to include different regions
bboxes <- read.csv(
  file.path(mFolder, "/required_files/OBSERVATIONS/LLdM_AOI2/BoundingBoxClimaZones.csv"),
  stringsAsFactors=F)

# rFolderMain = paste0("/Results/LLdM_test/")  # Results Folder

gcmFolders = "GCMs/" #df3$File[1]
enso_signal_file = "/required_files/Other/ONI_SIGNAL.csv" # ENSO EL NINO-LA NINA CYCLE DATA
svd_analysis_file = "/required_files/Other/SVD_Table_Comp1_Pac_Atl.csv" #Pacific and Atlantic SST Oscillation Data

df_dumm = read.csv(
  file.path(mFolder, 'required_files/Other/Archive_list_batch_LLdM_Nov2021 - CMIP6.csv'),
  header=T)
print(df_dumm)
df = subset(df_dumm, df_dumm$DoIt==1)
df$BaseDay <- NA
df$BaseMonth <- NA
df$BaseYear <- NA
df$LeapFix <- NA

models = unique(df$Model)
ensembles = unique(df$Ensemble)
experiments = unique(df$Experiment)
print(experiments)
fut_all = unique(df$Experiment)
fut_all = fut_all[fut_all != "Historical"]
var_all = c("pr")

if (!dir.exists(file.path(rFolder))) {
  dir.create(file.path(rFolder))
}
save(df, file = "df.Rda")

for(r in 1:1) { #1:nrow(bboxes) 
  region <- bboxes$Descripcion[r] # name of your basin or region for use in output graphics
  print(paste("REGION", region))
  bBox = c(
    latMin = bboxes$ymin[r],
    latMax = bboxes$ymax[r],
    lonMin = bboxes$xmin[r],
    lonMax = bboxes$xmax[r])  # Bounding Box of Study Region
  ## enter -999 for the bounding box values to use the shapefile to define the bounding box.
  
  #-----------------------------------------------------------
  # PART ONE: Loops through the catalog and reads the netcdfs and the local data (pr and tas):
  #-----------------------------------------------------------
  start1 <- Sys.time()
  
  for (model in models) {
    print(paste("MODEL", model))
    model_results_path = file.path(rFolder, model)
    # rFolder = paste0(paste0(rFolderMain, model, ""))
    if (!dir.exists(model_results_path)) {
      dir.create(model_results_path)
    }
    save(df, file = file.path(model_results_path, "df.Rda")) # copied main df to model folder?
    
    df1 = subset(df, Model == model)
    exp_model = unique(df1$Experiment)
    
    for (ensemble in ensembles) {
      print(paste("ENSEMBLE", ensemble))
      df2 = subset(df1, ensemble == Ensemble)
      
      for (varName in CLIMATE_VARS) {
        print(paste("VARIABLE", varName))
        df3 = varName #subset(df2, Variable == varName)
        mv = 1 #df3$Moving.Year[1]  #moving year to calculate transition matrix
        varLabel = varName ## (lazy, fix if wanted). .. "Precipitation" #df3$Alias[1]
        observation_files = observation_files_lock
        observation_files[1] = gsub("VARIABLE", varName, observation_files_lock[1])
        observation_files[2] = gsub("VARIABLE", varName, observation_files_lock[2])
        load(file.path(model_results_path, "df.Rda"), verbose = TRUE)
        
        ## Temp
        # experiment = experiment[1] #"historical"
        
        # res <- try(
        print("READ GCM NCDF")
        read_GCM_NCDF(boundingBox = bBox,
                               movingDateStep = mv,
                               main_folder = mFolder,
                               results_folder = model_results_path,
                               Experiments = experiments,
                               Ensemble = ensemble,
                               variableName = varName,
                               variableLabel = varLabel,
                               shp_file = shp_file_background,
                               obs_files = observation_files,
                               sta_bbox = "region", #historical data from the bbox of the selected GCM pixel ("gcm pixel"), or from the user-defined bbox ("region")
                               minObsYear = 1981,
                               minGCMYear = as.numeric(substr(as.character(df$Begin),1,4)))
        # )
        
        # if(inherits(res, "try-error")) {
        #   try(dev.off())
        # }
        load(file.path(rFolder, "df.Rda"))
        # print(df)
      } # next variable
    } # next ensemble
 
    
    load(file.path(rFolder, "df.Rda"))
    # print(df)
    df1 = subset(df, Model == model)
    exp_model = unique(df1$Experiment)
    exp_model <- exp_model[exp_model != "historical"] # excludes historical experiment from following functions to just have futures
    
    #---------------------------------------------------------------------------
    # comparison reports between GCM and observed climate
    # for (ensemble in ensembles){
    #   print(paste("ENSEMBLE", ensemble))
    #   #Ensemble <- levels(ensemble)[1]
    #   rcp_pr_ind <- which(df$Model == model & df$Ensemble == ensemble) # just to find base date of a future scenario
    #   print(experiment)
    #   rcp_tas_ind <- which(df$Model == model & df$Experiment == experiment)
    #   
    #   for (experiment2 in exp_model) {
    #     print(paste("EXPERIMENT", experiment2))
    #     
    #     print("COMPARE GCM")
    #     # try(
    #     compare_GCM(refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind]),
    #                   main_folder = mFolder,
    #                   results_folder = rFolder,
    #                   enso_signal_file = enso_signal_file,
    #                   svd_analysis_file = svd_analysis_file,
    #                   modelName = model,
    #                   historical = "historical",
    #                   futures = experiment2, #experiment,
    #                   varName = "pr",
    #                   varLabel = "Precipitation (mm)",
    #                   minObsYear = 1981,
    #                   maxObsYear = 2016,
    #                   minGCMYear = 2015,
    #                   maxGCMYear = 2099,
    #                   alignHistYears = TRUE)
    #     # )
    # 
    #   } # next experiment
    # } # next n ensemble member
  } #next model     
  
  #-----------------------------------------------------------
  # PART TWO: Runs the downscaling routines
  #-----------------------------------------------------------
  
  print("PART TWO DOWNSCALING")
  for (model in models) {
    print(paste("MODEL", model))
    model_results_path = file.path(rFolder, model)
    load(file.path(rFolder, "df.Rda"))
    df1 = subset(df, Model == model)
    exp_model = unique(df1$Experiment)
    exp_model <- exp_model[exp_model != "historical"] # excludes historical experiment from following functions to just have futures

    # rFolder = paste0(paste0(rFolderMain, model, ""))

    # Optional: update of extremes Obs based (data_d) of GCM Results.
    # The future reference year represents the upper boundary year for which new extremes are derived

    for (experiment2 in exp_model) {
      print(paste("EXPERIMENT", experiment2))
      # refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind])
      experiment_date = df[
        df$Variable == "pr" & df$Model == model & df$Experiment == experiment2,
        c("BaseDay", "BaseMonth", "BaseYear")]
      refDate = c(
        month = experiment_date$BaseMonth,
        day = experiment_date$BaseDay,
        year = experiment_date$BaseYear)

     # Legacy function. now the extreme value correction is done inside the KNN_BOOTSTRAP function
     #
     # try(extremePR_GPD(modelName = model,
     #                  futures = experiment2,
     #                  main_folder = mFolder,
     #                  results_folder = rFolder,
     #                  minObsYear = 1980,
     #                  minGCMYear = 2015,
     #                  maxGCMYear = 2099))

      for (n in 1:BootstrapRuns) {
        if ("pr" %in% CLIMATE_VARS) {
          print(paste("KNN BOOTSTRAP n =", n))
          # try(
          knn_bootstrap(refDate = refDate,
                        main_folder = mFolder,
                        results_folder = model_results_path,
                        modelName = model,              #"CCSM4",
                        futures = experiment2,
                        varName = "pr",
                        varLabel = "Precipitation [mm]",
                        minObsYear = 1980,
                        maxObsYear = 2015,
                        nearWindow = 30,
                        minGCMYear = 2020,
                        maxGCMYear = 2050,
                        expNumber = n,
                        JPmode = "Window",
                        alignHistYears = TRUE,
                        HistRepro = FALSE)
          # )
        }
        if ("tas" %in% CLIMATE_VARS) {
          print("HEAT SIGNAL BOOTSTRAP n =", n)
          # try(
          heatSignal(refDate = refDate,
                     main_folder = mFolder,
                     results_folder = model_results_path,
                     modelName = model,
                     futures = experiment2,
                     varName = "tas",
                     varLabel = "Temperature [c]",
                     minObsYear = 1980,
                     maxObsYear = 2015,
                     minGCMYear = 2020,
                     maxGCMYear = 2050)
          # )  
        }
        


        # try(sinteticSeries(refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind]),
        #                    modelName = model,
        #                    futures = experiment2,
        #                    main_folder = mFolder,
        #                    results_folder = rFolder,
        #                    expNumber = n,
        #                    fullObsCatalog = FALSE,
        #                    HistRepro = FALSE,
        #                    minObsYear = 1980,
        #                    maxObsYear = 2015))

      }  # next bootstrap
    } #next experiment
  } # next model
} # next region