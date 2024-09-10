#1 Read data files 
##Control 
ClinicalTestControl <- read.csv("Supplementary\\Control\\clinical_tests_input_control.csv")
CytokineControl <- as.matrix(read.csv("Supplementary\\Control\\cytokine_abundance_input_control.csv", row.names = 1))
Gut16Control <- as.matrix(read.csv("Supplementary\\Control\\gut_16s_abundance_input_control.csv", row.names = 1))
MetabolomeControl <- read.csv("Supplementary\\Control\\metabolome_abundance_input_control.csv")
Nares16Control <-as.matrix(read.csv("Supplementary\\Control\\nares_16s_abundance_input_control.csv", row.names = 1))
proteomeControl <-read.csv("Supplementary\\Control\\proteome_abundance_input_control.csv")
RNAseqControl <-as.matrix(read.csv("Supplementary\\Control\\RNAseq_abundance_input_control.csv", row.names = 1))
##prediabetic
ClinicalPre <- read.csv("Supplementary\\Prediabetic\\clinical_tests_input.csv")
CytokinePre<- as.matrix(read.csv("Supplementary\\Prediabetic\\cytokine_abundance_input.csv", row.names = 1))
Gut16Pre<- as.matrix(read.csv("Supplementary\\Prediabetic\\gut_16s_abundance_input.csv", row.names = 1))
MetabolomePre <-read.csv("Supplementary\\Prediabetic\\metabolome_abundance_input.csv")
Nares16Pre <-as.matrix(read.csv("Supplementary\\Prediabetic\\nares_16s_abundance_input.csv", row.names = 1))
ProteomePre <-read.csv("Supplementary\\Prediabetic\\proteome_abundance_input.csv")
RNAseqpre <-as.matrix(read.csv("Supplementary\\Prediabetic\\RNAseq_abundance_input.csv", row.names = 1))

#duplicated row names in clinical, Metabolome, proteome
duplicated_rows <- duplicated(ClinicalTestControl[, 1])  
sum(duplicated_rows) 
#to know which row is duplicated
duplicated_rows <- duplicated(rownames(ModifedClinicalTest))
print(duplicated_rows)

#duplicated
ModifedClinicalTest <-head(ClinicalTestControl, -1)
write.csv(ModifedClinicalTest, "ModifedData\\Control\\ModifedClinicalTest.csv", row.names = FALSE)
ModifedClinicalPre <-head(ClinicalPre,-1)
write.csv(ModifedClinicalPre, "ModifedData\\Control\\ModifedClinicalPre.csv", row.names = FALSE)

ModifedMetabolomeControl <-head(MetabolomeControl, -1)
write.csv(ModifedMetabolomeControl, "ModifedData\\Control\\ModifedMetabolomeControl.csv", row.names = FALSE)
ModifedMetabolomePre <-head(MetabolomePre, -1)
write.csv(ModifedMetabolomePre, "ModifedData\\Control\\ModifedMetabolomePre.csv", row.names = FALSE)

ModifedProteomeControl <-head(proteomeControl,-1)
write.csv(ModifedProteomeControl, "ModifedData\\Control\\ModifedProteomeControl.csv", row.names = FALSE)
ModifedProteomePre <-head(ProteomePre,-1)
write.csv(ModifedProteomePre, "ModifedData\\Control\\ModifedProteomePre.csv" , row.names = FALSE)


#no of rows in pre and control must be the same
nrow(ModifedClinicalTest)
nrow(ModifedClinicalPre)
min_rows <- min(nrow(ModifedClinicalTest), nrow(ModifedClinicalPre))
ModifedClinicalTest <- ModifedClinicalTest[1:min_rows, ]
ModifedClinicalPre <- ModifedClinicalPre[1:min_rows, ]


ModifedClinicalTest <- as.matrix(read.csv("ModifedData\\Control\\ModifedClinicalTest.csv", row.names = 1))
ModifedClinicalPre<-as.matrix(read.csv("ModifedData\\Control\\ModifedClinicalPre.csv", row.names = 1))
ModifedMetabolomeControl<-as.matrix(read.csv("ModifedData\\Control\\ModifedMetabolomeControl.csv", row.names = 1))
ModifedMetabolomePre<-as.matrix(read.csv("ModifedData\\Control\\ModifedMetabolomePre.csv", row.names = 1))
ModifedProteomeControl<-as.matrix(read.csv("ModifedData\\Control\\ModifedProteomeControl.csv", row.names = 1))
ModifedProteomePre<-as.matrix(read.csv("ModifedData\\Control\\ModifedProteomePre.csv", row.names = 1))

#bind
Clinical <-cbind(ModifedClinicalTest,ModifedClinicalPre)
Cytokine <- cbind(CytokineControl, CytokinePre )
Gut16 <- cbind(Gut16Control, Gut16Pre)
Metabolome <-cbind(ModifedMetabolomeControl, ModifedMetabolomePre)
Nares16 <-cbind(Nares16Control,Nares16Pre)
Proteome <-cbind(ModifedProteomeControl, ModifedProteomePre)
RNAseq <- cbind(RNAseqControl, RNAseqpre)

View(Clinical)
View(Cytokine)
View(Gut16)
View(Metabolome)
View(Nares16)
View(Proteome)
View(RNAseq)

library(MOFA2)

data <- make_example_data(
  n_views = 6, 
  n_samples = 200, 
  n_features = 1000, 
  n_factors = 6
)[[1]]
View(data)
names(data)<-c("Cytokine", "Gut16", "Metabolome", "Nares16", "Proteome", "RNAseq")
data[["Cytokine"]]<-Cytokine
data[["Gut16"]]<-Gut16
data[["Metabolome"]] <-Metabolome
data[["Nares16"]]<-Nares16
data[["Proteome"]]<-Proteome
data[["RNAseq"]]<-RNAseq

lapply(data,dim)

MOFAobject <- create_mofa(data)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

reticulate::py_config()
reticulate::use_python()

# Example of running MOFA with basilisk-enabled environment
MOFAmodel <- run_mofa(MOFAobject, use_basilisk = TRUE)

library(rhdf5)
h5createFile("myhdf5file.h5")
MOFAmodel<-run_mofa(MOFAobject,"myhdf5file.h5" )
