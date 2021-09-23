MOFAobject <- create_mofa(CRCdata)

## plot sample overall info
plot_data_overview(MOFAobject,show_dimensions = T,
                   colors = c("#2ca25f",'#756bb1',"#2c7fb8","#feb24c","#f03b20")) + 
  labs(title = "Sample Distribution of Omics data") +
  theme(plot.title = element_text(hjust = 0.5, size = 16,face = "bold"))+
  annotate("text", x=30, y=1, label="N = 56") +
  annotate("text", x=30, y=2, label="N = 50") +
  annotate("text", x=30, y=3, label="N = 54") +
  annotate("text", x=30, y=4, label="N = 55") +
  annotate("text", x=30, y=5, label="N = 45")


## define MOFA options
dataOpts <- get_default_data_options(MOFAobject)
dataOpts$scale_views <- FALSE
dataOpts

modelopts <- get_default_model_options(MOFAobject)
modelopts$ard_factors <- TRUE
modelopts$ard_weights <- FALSE
modelopts$num_factors <- 6
modelopts

trainOpts <- get_default_training_options(MOFAobject)
# trainOpts$drop_factor_threshold <- 0.02
trainOpts$seed <- 2
trainOpts$maxiter <- 3000
trainOpts$convergence_mode <- "fast"
trainOpts

MOFAobject1 <- prepare_mofa(MOFAobject,
                            data_options = dataOpts,
                            model_options = modelopts,
                            training_options = trainOpts) 
MOFAobject1

MOFAobject2 <- run_mofa(MOFAobject1,outfile = "5_MOFA2/CRC_seed2.hdf5",
                        save_data = TRUE)

saveRDS(MOFAobject2,"5_MOFA2/CRC_seed2.hdf5")
