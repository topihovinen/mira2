##    This is an R code for statistical analysis and generating figures of MIRA2 study            ##
# ------------------------------------------------------------------------------------------------ #
#             .-"''-.  _                                                                           #
# .'       `( \                                    __  __ ___ ____      _    __                    #
#         @/            ')   ,--,__,-"            |  \/  |_ _|  _ \    / \   __)                   #
#         /        /      \ /     /   _/          | |\/| || || |_) |  / _ \ (__                    #
#       __|           ,   |/         /            | |  | || ||  _ <  / ___ \                       #
#     .~  `\   / \ ,  |   /                       |_|  |_|___|_| \_\/_/   \_\                      #
#   .~      `\    `  /  _/   _/                                                                    #
# .~          `\  ~~`__/    /                                                                      #
# ~             `--'/                                                                              #
#              /   /    /                                                                          #
#             /  /'    /jgs                                                                        #
#                                                              http://github.com/topihovinen/mira2 #
# ------------------------------------------------------------------------------------------------ #

{
  # -------------------------------------------0. PREFACE------------------------------------------- #
  {
    # Introduce needed libraries for R
    message("-------This file includes data import, curation and analyses for MIRA2 study.-------")
    # Import required libraries
    {
      message("Importing needed libraries....", appendLF = FALSE)
      {
        setwd(dirname(rstudioapi::getSourceEditorContext()$path))  # This requires RStudio as your API to work
        library("ggplot2")
        library("ggbiplot")
        library("Hmisc")
        library("ggbreak")
        library("dplyr")
        library("ggsignif")
        library("png")
        library("openxlsx")
        library("showtext")
        source("mira_cgdiet_function.R")
        source("mira2_perform_permutation_fdr.R")
      }
      message("done.")
    }
    # Import and preprocess data
    {
      # Import the "targeted" metabolomics data and master metadata, and use exclusion criteria
      message("Importing metadata and targeted metabolomics....", appendLF = FALSE)
      {
        mira2data <- read.csv("mira2_huslab_R.txt", sep='\t', dec=',', header=T) # Start with data frame including main blood values
        mira2data <- subset(mira2data, select = -c(Comments)) #Discard the column with text comments
        # Import, merge and curate masterdata, use masterdata for exclusion
        {
          # Import and merge masterdata (not including clinical data which is analysed separately due to sensitive data)
          {
            mira2master <- read.csv("master/participants_master.csv", sep=';', dec=',', header=T)
            mira2master <- subset(mira2master, select = -c(age_blood, sex, diet_reported))
            mira2data <- merge(mira2data, mira2master, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            rm(mira2master)
          }
          # Exclude, criteria: cancelled participation, no food record or no blood sample
          {
            mira2data <- mira2data[((mira2data$is_participant) & (mira2data$diet_group_datadriven != "") & (!is.na(mira2data$age_blood))),] #Discard all subjects who cancelled their participation (no date of birth as a marker for that)
          }
          # Assign new variables for pregnancy, breastfeeding and preterm (manually imported), and discard unused variables (e.g. multiple different diet groups with slightly different logic)
          {
            mira2data$pregnant <- !is.na(mira2data$pregnancy_week)
            mira2data$breastfeed <- !is.na(mira2data$breastfeeding_frequency)
            mira2data$preterm_weeks[mira2data$id %in% c("M3333", "M3416", "M3497")] = c(32, 36, 35)
            # Convert diets during pregnancy to same groups as diet_group_datadriven. Dictionary for this translation was validated manually based on additional open answers in original questionnaire
            {
              mira2data$diet_preg_ticks <- as.character(mira2data$diet_preg_ticks)
              mira2data$diet_preg_ticks[mira2data$diet_preg_ticks %in% c("1000000000", "1000000010", "100000000", "100")] <- "mixed diet"
              mira2data$diet_preg_ticks[mira2data$diet_preg_ticks %in% c("100010000", "1010000", "10000")] <- "vegetarian"
              mira2data$diet_preg_ticks[mira2data$diet_preg_ticks %in% c("10100000", "100010", "100000")] <- "vegan"
              mira2data$diet_preg_ticks[mira2data$id %in% c("M3552", "M3607")] <- "mixed diet"
              mira2data$diet_preg_ticks[mira2data$id %in% c("M3638")] <- "vegan"
            }
            # Convert diet variables to factors
            {
              mira2data$diet_group_datadriven <- factor(mira2data$diet_group_datadriven, levels = c("mixed diet", "vegetarian", "vegan"))
              mira2data$diet_preg_ticks <- factor(mira2data$diet_preg_ticks, levels = c("mixed diet", "vegetarian", "vegan"))
            }
            mira2data <- subset(mira2data, select = c(id, id_family, sex, diet_group_datadriven, diet_preg_ticks, age_blood, blood_time_of_day, blood_month_of_year, height, weight, muac, head_circ, height_cg, weight_cg, education, education_family_max, diet_home_always, physical_activity, smoking, alcohol_use, B.Leuk, B.Eryt, B.Hb, B.HKR, E.MCV, E.RDW, E.MCH, E.MCHC, B.Trom, fP.Kol.HDL, fP.Kol.LDL, fP.Kol, fP.Trigly, P.Alb, fP.LipoA1, fP.LipoB, fP.ApoBperA1, S.Prealb, S.B12.TC2, pregnant, breastfeed, preterm_weeks))
          }
          # Assign cohorts (child, father, mother)
          {
            mira2data$cohort <- (mira2data$age_blood < 10) # First, all under 10yo are "TRUE" and other "FALSE"
            mira2data$cohort[mira2data$cohort == TRUE] = "child" # change "TRUE" to "child"
            mira2data$cohort[mira2data$cohort == FALSE] = "adult, F" # change "FALSE" first to "adult, F"
            mira2data$cohort[((mira2data$cohort == "adult, F") & (mira2data$sex == "M"))] = "adult, M" # check which "adult, F"'s are considered male and change to "adult, M"
            mira2data$cohort <- factor(mira2data$cohort, levels = c("child", "adult, F", "adult, M"))
          }
        }
        # Import rest of the separate data frames and merge to main data
        {
          # Combine different amino acid data and merge
          {
            mira2aminos <- read.csv("mira2_aminoacids_R.txt", sep='\t', dec=',', header=T)
            mira2taurine <- read.csv("mira2_taurine_R.txt", sep='\t', dec=',', header=T)
            mira2aminos_adult <- read.csv("mira2_aminoacids_adult_R.txt", sep = '\t', dec=',', header = T)
            mira2aminos <- merge(mira2taurine, mira2aminos, by.x = "id", by.y = "ID", all = TRUE)
            mira2aminos <- merge(mira2aminos, mira2aminos_adult, all = TRUE)
            mira2data <- merge(mira2data, mira2aminos, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE) # Included adult measurements as of 2.4.2023
            rm(mira2aminos, mira2aminos_adult, mira2taurine)
          }
          # Import and merge the rest
          {
            mira2nadgsh <- read.csv("mira2_nadgsh_R.txt", sep='\t', dec=',', header=T)
            mira2vitd <- read.csv("mira2_vitd_R.txt", sep='\t', dec=',', header=T)
            mira2ironrbp <- read.csv("mira2_ironrbp_R.txt", sep='\t', dec=',', header=T)
            mira2bas <- read.csv("mira2_bileacids_R.txt", sep='\t', dec=',', header=T)
            mira2vitae <- read.csv("mira2_vit_ae_R.txt", sep='\t', dec=',', header=T)
            mira2sterols <- read.csv("mira2_sterols_R.txt", sep='\t', dec=',', header=T)
            mira2intake_suppl <- read.csv("mira2_suppl_intake_R.txt", sep='\t', dec=',', header = T)
            mira2intake_suppl <- subset(mira2intake_suppl, select = -c(X))
            colnames(mira2intake_suppl) <- paste("intake_", colnames(mira2intake_suppl), "_s", sep = '')
            mira2intake_food <- read.csv("mira2_asep_intake_R.txt", sep='\t', dec=',', header=T)
            mira2fattyacids <- read.csv("mira2_fattyacids_R.txt", sep='\t', dec=',', header=T)
            mira2iodine <- read.csv("mira2_iodine_R.txt", sep='\t', dec=',', header=T)
            mira2bone <- read.csv("mira2_bonemarker_R.txt", sep='\t', dec=',', header=T)
            mira2data <- merge(mira2data, mira2nadgsh, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2vitd, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2ironrbp, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2bas, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2vitae, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2sterols, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2intake_suppl, by.x = "id", by.y = "intake_id_s", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2intake_food, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2fattyacids, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2iodine, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            mira2data <- merge(mira2data, mira2bone, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
            rm(mira2nadgsh, mira2vitd, mira2ironrbp, mira2bas, mira2vitae, mira2sterols, mira2intake_suppl, mira2intake_food, mira2fattyacids, mira2iodine, mira2bone)
          }
        }
      }
      message("done. Found ", length(colnames(mira2data)), " variables and ", length(mira2data$id), " participants:\n   ",
              sum(mira2data$age_blood < 10), " children\n   ",
              sum(((mira2data$age_blood > 10) & (mira2data$sex == "F"))), " adult women, of which\n   ...",
              sum(((mira2data$age_blood > 10) & (mira2data$sex == "F") & (mira2data$pregnant))), " pregnant\n   ...",
              sum(((mira2data$age_blood > 10) & (mira2data$sex == "F") & (mira2data$breastfeed))), " breastfeeding\n   ",
              sum(((mira2data$age_blood > 10) & (mira2data$sex == "M"))), " adult men.")
      # Calculate secondary variables
      message("Curating targeted data....", appendLF = FALSE)
      {
        # Amino acid combination variables
        {
          mira2data$aa.essential <- mira2data$HIS + mira2data$ILE + mira2data$LEU + mira2data$LYS + 
            mira2data$MET + mira2data$PHE + mira2data$THR + mira2data$TRP + mira2data$VAL
          for(id in mira2data$id[!is.na(mira2data$GLN)]) {
            mira2data$aa.nonessential[mira2data$id == id] <- mira2data$ALA[mira2data$id == id] + mira2data$ARG[mira2data$id == id] + mira2data$ASN[mira2data$id == id] +
              mira2data$ASP[mira2data$id == id] + mira2data$CYS[mira2data$id == id] + mira2data$GLU[mira2data$id == id] + mira2data$GLN[mira2data$id == id] +
              mira2data$GLY[mira2data$id == id] + mira2data$PRO[mira2data$id == id] + mira2data$SER[mira2data$id == id] + mira2data$TYR[mira2data$id == id]
          }
          for(id in mira2data$id[!is.na(mira2data$GLN.ARG.comb)]) {
            mira2data$aa.nonessential[mira2data$id == id] <- mira2data$ALA[mira2data$id == id] + mira2data$ASN[mira2data$id == id] +
              mira2data$ASP[mira2data$id == id] + mira2data$CYS[mira2data$id == id] + mira2data$GLU[mira2data$id == id] + mira2data$GLN.ARG.comb[mira2data$id == id] +
              mira2data$GLY[mira2data$id == id] + mira2data$PRO[mira2data$id == id] + mira2data$SER[mira2data$id == id] + mira2data$TYR[mira2data$id == id]
          }
          rm(id)
          mira2data$aa.protgen.total <- mira2data$aa.essential + mira2data$aa.nonessential
          mira2data$aa.sulf.total <- mira2data$MET + mira2data$CYS
          mira2data$aa.bcaa <- mira2data$ILE + mira2data$LEU + mira2data$VAL
          mira2data$aa.ess.per.ness <- mira2data$aa.essential / mira2data$aa.nonessential
          mira2data$aa.catechol.pc <- mira2data$TYR + mira2data$PHE
          mira2data$aa.charge.pos <- mira2data$ARG + mira2data$HIS + mira2data$LYS
          mira2data$aa.charge.neg <- mira2data$ASP + mira2data$GLU
          mira2data$aa.aromatic <- mira2data$PHE + mira2data$HIS + mira2data$TRP + mira2data$TYR
          mira2data$aa.aliphatic <- mira2data$ALA + mira2data$ILE + mira2data$LEU + mira2data$VAL + mira2data$MET
          mira2data$aa.glypertau <- mira2data$GLY / mira2data$TAUR
          mira2data$aa.urea <- mira2data$ARG + mira2data$ORN
        }
        # Add combined bile acid information
        {
          mira2data$CApool <- mira2data$CA + mira2data$TCA + mira2data$GCA
          mira2data$CDCApool <- mira2data$CDCA + mira2data$TCDCA + mira2data$GCDCA
          mira2data$DCApool <- mira2data$DCA + mira2data$TDCA + mira2data$GDCA
          mira2data$LCApool <- mira2data$LCA + mira2data$TLCA + mira2data$GLCA
          mira2data$UDCApool <- mira2data$UDCA + mira2data$TUDCA + mira2data$GUDCA
          mira2data$BA_primary_total <- mira2data$CA + mira2data$TCA + mira2data$GCA + mira2data$CDCA + mira2data$TCDCA + mira2data$GCDCA
          mira2data$BA_primary_unconj <- mira2data$CA + mira2data$CDCA
          mira2data$BA_primary_conj <- mira2data$TCA + mira2data$GCA + mira2data$TCDCA + mira2data$GCDCA
          mira2data$BA_primary_tauro <- mira2data$TCA + mira2data$TCDCA
          mira2data$BA_primary_glyco <- mira2data$GCA + mira2data$GCDCA
          mira2data$BA_primary_taupergly <- mira2data$BA_primary_tauro/mira2data$BA_primary_glyco
          mira2data$BA_secondary_total <- mira2data$DCA + mira2data$TDCA + mira2data$GDCA + mira2data$LCA + mira2data$TLCA + mira2data$GLCA + mira2data$UDCA + mira2data$TUDCA + mira2data$GUDCA
          mira2data$BA_secondary_unconj <- mira2data$DCA + mira2data$LCA + mira2data$UDCA
          mira2data$BA_secondary_conj <- mira2data$TDCA + mira2data$GDCA + mira2data$TLCA + mira2data$GLCA + mira2data$TUDCA + mira2data$GUDCA
          mira2data$BA_secondary_tauro <- mira2data$TDCA + mira2data$TLCA + mira2data$TUDCA
          mira2data$BA_secondary_glyco <- mira2data$GDCA + mira2data$GLCA + mira2data$GUDCA
          mira2data$BA_secondary_taupergly <- mira2data$BA_secondary_tauro/mira2data$BA_secondary_glyco
          mira2data$BA_total_total <- mira2data$CA + mira2data$TCA + mira2data$GCA + mira2data$CDCA + mira2data$TCDCA + mira2data$GCDCA + mira2data$DCA + mira2data$TDCA + mira2data$GDCA + mira2data$LCA + mira2data$TLCA + mira2data$GLCA + mira2data$UDCA + mira2data$TUDCA + mira2data$GUDCA
          mira2data$BA_total_unconj <- mira2data$CA + mira2data$CDCA + mira2data$DCA + mira2data$LCA + mira2data$UDCA
          mira2data$BA_total_conj <- mira2data$TCA + mira2data$GCA + mira2data$TCDCA + mira2data$GCDCA + mira2data$TDCA + mira2data$GDCA + mira2data$TLCA + mira2data$GLCA + mira2data$TUDCA + mira2data$GUDCA
          mira2data$BA_total_tauro <- mira2data$TCA + mira2data$TCDCA + mira2data$TLCA + mira2data$TDCA + mira2data$TUDCA
          mira2data$BA_total_glyco <- mira2data$GCA + mira2data$GCDCA + mira2data$GLCA + mira2data$GDCA + mira2data$GUDCA
          mira2data$BA_total_taupergly <- mira2data$BA_total_tauro/mira2data$BA_total_glyco
        }
        # Calculate creatinine corrected urine iodine
        {
          mira2data$U.I.per.Crea <- mira2data$U.I / mira2data$U_CREA_MMOL_L
        }
      }
      # Curation of data
      {
        # For B12 and prealbumin/transthyretin there are one or two values exceeding the measurement range, replace them with the limit value.
        {
          mira2data$S.B12.TC2[mira2data$S.B12.TC2 == ">438"] <- "438"
          mira2data$S.B12.TC2 <- as.numeric(mira2data$S.B12.TC2)
          mira2data$S.Prealb[mira2data$S.Prealb == "<100"] <- "100"
          mira2data$S.Prealb <- as.numeric(mira2data$S.Prealb)
        }
        # Curate CRP data (includes a lot of "<0.2" by default)
        {
          mira2data$S.CRP[mira2data$S.CRP == "<0.2"] <- 0
          mira2data$S.CRP <- as.numeric(as.numeric(gsub(",", ".", mira2data$S.CRP)))
        }
        # RBP, iron and inflammation measurements were mixed in the original data file with two adjacent samples (confirmed from correlations on 11.1.2023) so this will swap the values
        {
          marker_id <- which(colnames(mira2data) %in% c("S.Ferrit", "S.TfR", "Iron.Storage.mg.kgbw", "S.RBP", "S.CRP", "S.AGP"))
          person_id <- which(mira2data$id %in% c("M3608", "M3609"))
          mira2data[person_id, marker_id] <- mira2data[rev(person_id), marker_id]
          rm(person_id, marker_id)
        }
      }
      message("done.")
      # Import untargeted metabolomics and merge to main data
      message("Importing untargeted metabolomics....", appendLF = FALSE)
      {
        mira2_metab <- read.csv("mira2_metabolomics_Rnew.csv", sep=';', dec = ",", header = T, stringsAsFactors=FALSE)
        mira2_metab_idkey <- read.csv("mira2_metabolomics_id_key_R.txt", sep = '\t', dec = ',', header = T)
        mira2_metab <- merge(mira2_metab, mira2_metab_idkey, by.x = "eth_id", by.y = "eth_id")
        # Use mean of duplicate measurements as final data
        {
          mira2_metab_temp <- mira2_metab_temp2 <- data.frame(id = unique(mira2_metab$id))
          mira2_metab_temp2 <- (mira2_metab[seq(1,length(mira2_metab$eth_id),2),5:(length(mira2_metab[1,])-1)] + mira2_metab[seq(2,length(mira2_metab$eth_id),2),5:(length(mira2_metab[1,])-1)]) / 2
          mira2_metab <- cbind(mira2_metab_temp, mira2_metab_temp2)
        }
        mira2data <- merge(mira2data, mira2_metab, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
        metab_n = length(colnames(mira2_metab)) - 1
        rm(mira2_metab, mira2_metab_temp, mira2_metab_temp2, mira2_metab_idkey)
      }
      message("done. Found ", metab_n, " features.")
      rm(metab_n)
      # Import growth data to a separate data frame due to the shape of data
      message("Importing growth/anthropometrics data to a separate data frame....", appendLF = FALSE)
      {
        mira2growth <- read.csv("mira2_growthdata_R_160525.txt", sep ="\t", dec = ",", header = T)
        mira2growth <- mira2growth[mira2growth$id %in% mira2data$id,]
        for(i in 1:length(mira2growth$id)) {
          if(mira2growth$id[i] %in% mira2data$id) {
            mira2growth$diet_group_datadriven[i] <- mira2data$diet_group_datadriven[mira2data$id == mira2growth$id[i]]
            mira2growth$diet_pregnancy[i] <- mira2data$diet_preg_ticks[mira2data$id == mira2growth$id[i]]
            mira2growth$preterm[i] <- !is.na(mira2data$preterm_weeks[mira2data$id == mira2growth$id[i]])
            mira2growth$sex[i] <- mira2data$sex[mira2data$id == mira2growth$id[i]]
            mira2growth$obs_n[i] <- sum(mira2growth$id == mira2growth$id[i])
          }
        }
        mira2growth$diet_group_datadriven <- factor(mira2growth$diet_group_datadriven, labels = levels(mira2data$diet_group_datadriven))
        mira2growth$diet_pregnancy <- factor(mira2growth$diet_pregnancy, labels = levels(mira2data$diet_preg_ticks))
        mira2growth$sex <- factor(mira2growth$sex)
      }
      message("done. Found ", i, " data points.")
      rm(i)
    }
  }
  # --------------------------------------1. ANTHROPOMETRIC DATA----------------------------------- #
  {
    # Linear regressions of anthropometric data
    {
      # Build linear models
      { # Linear weighted (by observation number per participant) regression models for height z-values corrected by parents' heights
        lm_heightz_curated_w <- lm(height_z ~ diet_group_datadriven + height_est_adult, mira2growth[(mira2growth$age_anthro < 10),], weights = obs_n)
        # Linear regression model of birthweight, predicted by mothers' diet during pregnancy
        lm_birthweight <- lm(weight ~ diet_pregnancy, mira2growth[(mira2growth$age_anthro == 0) & !mira2growth$preterm,])
      }
      # Form predicted lines for figures
      {
        lm_heightzcorr_mix_w <- lm(height_z_corr ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "mixed diet") & (mira2growth$age_anthro < 10),], weights = obs_n, na.action = "na.omit")
        lm_heightzcorr_veg_w <- lm(height_z_corr ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegetarian") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_heightzcorr_vgn_w <- lm(height_z_corr ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegan") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_pred_heightzcorr_mix_w <- data.frame(height_z_corr = predict(lm_heightzcorr_mix_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_heightzcorr_veg_w <- data.frame(height_z_corr = predict(lm_heightzcorr_veg_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_heightzcorr_vgn_w <- data.frame(height_z_corr = predict(lm_heightzcorr_vgn_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pipa_mix_w <- lm(pipa ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "mixed diet") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_pipa_veg_w <- lm(pipa ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegetarian") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_pipa_vgn_w <- lm(pipa ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegan") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_pred_pipa_mix_w <- data.frame(pipa = predict(lm_heightzcorr_mix_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_pipa_veg_w <- data.frame(pipa = predict(lm_heightzcorr_veg_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_pipa_vgn_w <- data.frame(pipa = predict(lm_heightzcorr_vgn_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_headcfz_mix_w <- lm(headcf_z ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "mixed diet") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_headcfz_veg_w <- lm(headcf_z ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegetarian") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_headcfz_vgn_w <- lm(headcf_z ~ age_anthro, mira2growth[(mira2growth$diet_group_datadriven == "vegan") & (mira2growth$age_anthro < 10),], weights = obs_n)
        lm_pred_headcfz_mix_w <- data.frame(headcfz = predict(lm_headcfz_mix_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_headcfz_veg_w <- data.frame(headcfz = predict(lm_headcfz_veg_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
        lm_pred_headcfz_vgn_w <- data.frame(headcfz = predict(lm_headcfz_vgn_w, data.frame(age_anthro = seq(0,7,0.01))), age_anthro = seq(0,7,0.01))
      }
    }
  }
  # --------------------------------------2. TARGETED BIOMARKERS----------------------------------- #
  {
    # Form meaningful groups of measurements for analyses
    {
      mira2_stat <- data.frame(i = 1:length(colnames(mira2data)),
                               group = c(rep("", length(colnames(mira2data)))),
                               name = colnames(mira2data),
                               kstest_p = c(rep(NA, length(colnames(mira2data)))),
                               shapirotest_p = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_child_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultF_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultM_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_child_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultF_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultM_vgn.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_child_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultF_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultM_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_child_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultF_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultM_vgn.vgt = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_child_vgt.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultF_vgt.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p_adultM_vgt.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_child_vgt.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultF_vgt.mix = c(rep(NA, length(colnames(mira2data)))),
                               wilcox_p.adj_adultM_vgt.mix = c(rep(NA, length(colnames(mira2data)))))
      mira2_stat$group[c(21:39, 79:87, 105:106, 208:211)] <- "clin"
      mira2_stat$group[c(44:66, 216:228)] <- "aa"
      mira2_stat$group[c(68:79)] <- "nad"
      mira2_stat$group[c(89:105, 229:251)] <- "ba"
      mira2_stat$group[c(110:125)] <- "cholmarkers"
      mira2_stat$group[c(159:207)] <- "fa"
      mira2_stat$group[c(253:length(mira2_stat$group))] <- "metab"
    }
    
    # Testing the normality of the data with Kolmogorov-Smirnov, using children as sufficiently indicative subsample
    {
      message("\nAnalyzing data normality using Kolmogorov-Smirnov test and Shapiro-Wilk test using subsample of child participants....", appendLF = FALSE)
      mira2data_normtest <- mira2data[mira2data$age_blood < 10,]
      i_totest <- mira2_stat$i[mira2_stat$group != ""]
      for(i in i_totest) {
        mira2_stat$kstest_p[i] = ks.test(mira2data_normtest[,i], "pnorm", mean(mira2data_normtest[,i], na.rm = T), sd(mira2data_normtest[,i], na.rm = T))$p.value
        if(sum(!is.na(mira2data_normtest[,i]) & (mira2data_normtest[,i] != 0)) >= 3) {
          mira2_stat$shapirotest_p[i] = shapiro.test(mira2data_normtest[,i])$p.value
        }
        if(i %% 1000 == 0) {
          message("past next 1000 tests....", appendLF = FALSE)
        }
      }
      # Summarize results, print them out and remove temps
      {
        mira2_stat_summary <- data.frame(group = c("All", "Clinical biomarkers", "Amino acids", "NAD and glutathiones", "Bile acids", "Cholesterol biomarkers", "Fatty acids", "Untargeted metabolomics"),
                                         n = c(sum(mira2_stat$group != ""), sum(mira2_stat$group == "clin", na.rm = TRUE), sum(mira2_stat$group == "aa", na.rm = TRUE), sum(mira2_stat$group == "nad", na.rm = TRUE), sum(mira2_stat$group == "ba", na.rm = TRUE), sum(mira2_stat$group == "cholmarkers", na.rm = TRUE), sum(mira2_stat$group == "fa", na.rm = TRUE), sum(mira2_stat$group == "metab", na.rm = TRUE)))
        {
          mira2_stat_summary$normality_ks <- c(sum(mira2_stat$kstest_p[mira2_stat$group != ""] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "All"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "clin"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Clinical biomarkers"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "aa"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Amino acids"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "nad"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "NAD and glutathiones"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "ba"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Bile acids"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "cholmarkers"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Cholesterol biomarkers"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "fa"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Fatty acids"],
                                               sum(mira2_stat$kstest_p[mira2_stat$group == "metab"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Untargeted metabolomics"])
          mira2_stat_summary$normality_sw <- c(sum(mira2_stat$shapirotest_p[mira2_stat$group != ""] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "All"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "clin"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Clinical biomarkers"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "aa"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Amino acids"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "nad"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "NAD and glutathiones"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "ba"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Bile acids"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "cholmarkers"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Cholesterol biomarkers"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "fa"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Fatty acids"],
                                               sum(mira2_stat$shapirotest_p[mira2_stat$group == "metab"] > 0.05, na.rm = T) / mira2_stat_summary$n[mira2_stat_summary$group == "Untargeted metabolomics"])
        }
        message("done.\n\nTotal ", mira2_stat_summary$n[1], " variables tested, of which ",
                round(mira2_stat_summary$normality_ks[1]*100, 0), "% or ", round(mira2_stat_summary$normality_sw[1]*100, 0), "% did not significantly differ from normal distribution in KS test and SW test, respectively. Specifically:")
        message("...", round(mira2_stat_summary$normality_ks[2]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[2]*100, 0), "% (SW) normality in clinical targeted biomarkers (n = ", mira2_stat_summary$n[2], ")")
        message("...", round(mira2_stat_summary$normality_ks[3]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[3]*100, 0), "% (SW) normality in amino acid measurements (n = ", mira2_stat_summary$n[3], ")")
        message("...", round(mira2_stat_summary$normality_ks[4]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[4]*100, 0), "% (SW) normality in NAD measurements (n = ", mira2_stat_summary$n[4], ")")
        message("...", round(mira2_stat_summary$normality_ks[5]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[5]*100, 0), "% (SW) normality in bile acid measurements (n = ", mira2_stat_summary$n[5], ")")
        message("...", round(mira2_stat_summary$normality_ks[6]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[6]*100, 0), "% (SW) normality in cholesterol metabolism biomarkers (n = ", mira2_stat_summary$n[6], ")")
        message("...", round(mira2_stat_summary$normality_ks[7]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[7]*100, 0), "% (SW) normality in fatty acid measurements (n = ", mira2_stat_summary$n[7], ")")
        message("...", round(mira2_stat_summary$normality_ks[8]*100, 0), "% (KS) or ", round(mira2_stat_summary$normality_sw[8]*100, 0), "% (SW) normality in untargeted metabolomics (n = ", mira2_stat_summary$n[8], ")")
        rm(mira2data_normtest, i, i_totest)
      }
    }
    
    # Proceed to statistical tests
    # Use Benjamini-Hochberg FDR method for multiple testing correction
    {
      message("\nProceeding to statistical testing between diet groups, separately testing for each cohort (children, women and men).")
      #for(cohort in c("child", "adult, F", "adult, M")) {
      for(markergroup in unique(mira2_stat$group)) {
        if(markergroup == "") {
          next
        }
        message("Testing for biomarkers in group ", markergroup, ".....", appendLF = FALSE)
        for(cohort in c("child", "adult, F", "adult, M")) {
          if(cohort == "child") {
            for(i in which(mira2_stat$group == markergroup)) {
              if(sum(!is.na(mira2data[(mira2data$cohort == cohort),i])) < 3*length(unique(mira2data$diet_group_datadriven))) {
                next
              }
              mira2_stat$wilcox_p_child_vgn.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_child_vgn.vgt[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_child_vgt.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
            }
          }
          if(cohort == "adult, F") {
            for(i in which(mira2_stat$group == markergroup)) {
              if(sum(!is.na(mira2data[(mira2data$cohort == cohort),i])) < 3*length(unique(mira2data$diet_group_datadriven))) {
                next
              }
              mira2_stat$wilcox_p_adultF_vgn.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_adultF_vgn.vgt[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_adultF_vgt.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
            }
          }
          if(cohort == "adult, M") {
            for(i in which(mira2_stat$group == markergroup)) {
              if(sum(!is.na(mira2data[(mira2data$cohort == cohort),i])) < 3*length(unique(mira2data$diet_group_datadriven))) {
                next
              }
              mira2_stat$wilcox_p_adultM_vgn.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_adultM_vgn.vgt[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegan") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i])$p.value
              mira2_stat$wilcox_p_adultM_vgt.mix[i] = wilcox.test(mira2data[((mira2data$diet_group_datadriven == "vegetarian") & (mira2data$cohort == cohort)),i], mira2data[((mira2data$diet_group_datadriven == "mixed diet") & (mira2data$cohort == cohort)),i])$p.value
            }
          }
        }
        # Adjust with B-H
        message("adjusting with Benjamini-Hochberg....", appendLF = FALSE)
        mira2_stat[mira2_stat$group == markergroup,c(9:11,15:17,21:23)] = p.adjust(unlist(mira2_stat[mira2_stat$group == markergroup,c(6:8,12:14,18:20)]), method = "BH")
        message("done.")
      }
      write.xlsx(mira2_stat, "Article misc for appendix/Statistical test results.xlsx")
      write.xlsx(mira2_metab_hmdb$hmdb_id[mira2_metab_hmdb$name %in% mira2_stat$name[mira2_stat$group == "metab" & mira2_stat$wilcox_p_child_vgn.mix < 0.05]], "Article misc for appendix/Pathway analysis, compoundlist for MetaboAnalyst, child, vgn-omn.xlsx")
    }
  }
  # --------------------------------------3. UNTARGETED METABOLOMICS----------------------------------- #
  {
    # Metabolomics analyses
    {
      #Import HMDB IDs
      {
        mira2data_mbx_hmdbdict <- read.xlsx("Metabolomics Zamboni/2022_MIRA2_Metabolomics_ETHZ_AO/MIRA2_Metabolomics_hmdb_key_R.xlsx", rowNames = TRUE)
      }
      mira2data_mbx_norm <- mira2data[,c(1:4,6:10,17:20,40:42,157,252:2376)]
      #Optional: Show HMDB IDs instead of variable names
      for(i in 1:length(colnames(mira2data_mbx_norm))) {
        if(colnames(mira2data_mbx_norm)[i] %in% rownames(mira2data_mbx_hmdbdict)) {
          colnames(mira2data_mbx_norm)[i] = mira2data_mbx_hmdbdict$hmdb_id[which(rownames(mira2data_mbx_hmdbdict) == colnames(mira2data_mbx_norm)[i])]
          next
        }
        if(colnames(mira2data_mbx_norm)[i] %in% paste("X", rownames(mira2data_mbx_hmdbdict), sep = '')) {
          colnames(mira2data_mbx_norm)[i] = mira2data_mbx_hmdbdict$hmdb_id[which(paste("X",rownames(mira2data_mbx_hmdbdict), sep = '') == colnames(mira2data_mbx_norm)[i])]
        }
      }
      mira2data_mbx_norm <- mira2data_mbx_norm[-c(which(is.na(mira2data_mbx_norm[,18]))),]
      mira2data_mbx_norm[,18:2141] <- scale(mira2data_mbx_norm[,18:2141])
      mira2data_mbx_pca <- prcomp(mira2data_mbx_norm[,18:2141])
      g <- ggbiplot(mira2data_mbx_pca,
                    choices = 2:3,
                    obs.scale = 1,
                    var.scale = 1,
                    groups = mira2data_mbx_norm$diet_group_datadriven,
                    ellipse = TRUE,
                    circle = TRUE,
                    varname.size = 0,
                    var.axes = FALSE,
                    ellipse.prob = 0.68) +
        scale_color_discrete(name = '') +
        theme(legend.direction = 'horizontal',
              legend.position = 'top')
      write.csv(mira2data_mbx_norm[mira2data_mbx_norm$age_blood < 10,c(1,17,18:2142)], "mira2_mbx_metaboanalyst_diet_child_hmdb.csv", row.names=FALSE)
      write.csv(mira2data_mbx_norm[mira2data_mbx_norm$age_blood > 10 & mira2data_mbx_norm$sex == "M",c(1,17,18:2142)], "mira2_mbx_metaboanalyst_diet_adultM_hmdb.csv", row.names=FALSE)
      write.csv(mira2data_mbx_norm[,c(1,17,18:2142)], "mira2_mbx_metaboanalyst_diet_all_hmdb.csv", row.names=FALSE)
      write.csv(mira2data_mbx_norm[mira2data_mbx_norm$age_blood > 10 & mira2data_mbx_norm$sex == "F",c(1,14,18:2142)], "mira2_mbx_metaboanalyst_diet_adultF_preggo_hmdb.csv", row.names=FALSE)
      
      # Wilcoxon tests for significant feature selection
      {
        for(c in c("child", "adult, F", "adult, M")) {
          message("Running for cohort: ", c, ".....", appendLF = F)
          mira2data_mbx_testhmdb <- mira2data_mbx_norm[mira2data_mbx_norm$cohort == c,]
          tests <- apply(mira2data_mbx_testhmdb[,18:2142], 2, pairwise.wilcox.test, g = mira2data_mbx_testhmdb$diet_group_datadriven, p.adjust.methods = "none", use = "complete.obs")
          means <- aggregate(mira2data_mbx_testhmdb[,18:2142], by = list(mira2data_mbx_testhmdb$diet_group_datadriven), FUN = mean, na.rm = T)
          tests_signs <- means[2:length(means[1,])]
          rownames(tests_signs) <- c(paste(means[1,1], "v", means[2,1]),
                                     paste(means[1,1], "v", means[3,1]),
                                     paste(means[3,1], "v", means[2,1]))
          means <- means[,2:length(means[1,])]
          tests_signs[1,] <- sign(means[2,] - means[1,])
          tests_signs[2,] <- sign(means[3,] - means[1,])
          tests_signs[3,] <- sign(means[2,] - means[3,])
          tests_pvalues <- sapply(tests, function(x) {
            p <- x$p.value
            n <- outer(rownames(p), colnames(p), paste, sep='v')
            p <- as.vector(p)
            names(p) <- n
            p
          })
          tests_signs[1,] <- tests_signs[1,] * tests_pvalues[1,]
          tests_signs[2,] <- tests_signs[2,] * tests_pvalues[2,]
          tests_signs[3,] <- tests_signs[3,] * tests_pvalues[4,] # Note: tests_pvalues will include row 3 for "vegan v vegan" comparison, yielding obvious row of NAs. Hence, the seeming discrepancy in indexes
          tests_signs <- data.frame(t(tests_signs))
          write.xlsx(tests_signs, paste("Article misc for appendix/Metabolomics p-values with signs, ", c, ".xlsx", sep = ""), rowNames = TRUE)
          message("done.")
        }
        rm(tests, means, tests_signs, mira2data_mbx_testhmdb,c)
      }
      # Use Robust Rank Aggregation for PLSDA VIP and Random Forest significant features -rank lists to select list of potentially interesting single metabolites
      {
        library("RobustRankAggreg")
        mira2_mbx_ranklists <- read.csv("MetaboAnalyst files/results/Results 10.4.2025/mira2_mbx_plsda_rf_ranks.csv", header = T, sep = ';')
        # Create list of vectors
        mira2_mbx_ranklists_rra <- list(mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$plsda_vip_child, decreasing = FALSE)],
                                        mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$plsda_vip_adultF, decreasing = FALSE)],
                                        mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$plsda_vip_adultM, decreasing = FALSE)],
                                        mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$rf_sigf_child, decreasing = FALSE)],
                                        mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$rf_sigf_adultF, decreasing = FALSE)],
                                        mira2_mbx_ranklists$metabolite[order(mira2_mbx_ranklists$rf_sigf_adultM, decreasing = FALSE)])
        mira2_mbx_rraresults <- aggregateRanks(mira2_mbx_ranklists_rra, method = "RRA", exact = T)
        write.xlsx(mira2_mbx_rraresults, "Article misc for appendix/Top 30 metabolites from RRA.xlsx")
      }
    }
  }
  # --------------------------------------MISCELLANEOUS----------------------------------- #
  {
    # --------------------------------------A. SELECTED EXAMPLES OF MANUAL COLLECTION FOR TABLE 1----------------------------------- #
    {
      # Pairs of siblings in mixed diet children
      sum(summary(factor(mira2data$id_family[mira2data$age_blood < 10 & mira2data$diet_group_datadriven == "mixed diet"])) > 1)
      # Mother's diet for children with mixed diet
      summary(factor(mira2data$diet_group_datadriven[mira2data$id_family %in% mira2data$id_family[mira2data$age_blood < 10 & mira2data$diet_group_datadriven == "mixed diet"] & mira2data$age_blood > 10 & mira2data$sex == "F"]))
      # Diet during pregnancy for children with mixed diet
      summary(factor(mira2data$diet_preg_ticks[mira2data$age_blood < 10 & mira2data$diet_group_datadriven == "mixed diet"]))
      # Height adjusted z-score at time of study for children with mixed diet
      summary(mira2growth$height_z_corr[mira2growth$is.neuvola == 0 & mira2growth$diet_group_datadriven == "mixed diet"])
      # Mother's education for children with mixed diet
      summary(factor(mira2data$education[mira2data$id_family %in% mira2data$id_family[mira2data$age_blood < 10 & mira2data$diet_group_datadriven == "mixed diet"] & mira2data$age_blood > 10 & mira2data$sex == "F"]))
      # Number of overweight women with mixed diet
      sum((mira2data$weight_cg[mira2data$age_blood > 10 & mira2data$sex == "F" & mira2data$diet_group_datadriven == "mixed diet"] / ((mira2data$height_cg[mira2data$age_blood > 10 & mira2data$sex == "F" & mira2data$diet_group_datadriven == "mixed diet"]/100)^2)) > 25)
    }
    
    # -------------------------------------B. MEDICAL HISTORY------------------------------------- #
    {
      mira2clin = read.csv("mira2_clinical_dg.txt", sep = '\t', dec = ',', header = T)
      mira2med = read.csv("mira2_clinical_med.txt", sep = '\t', dec = ',', header = T)
      mira2master <- read.csv("master/participants_master.csv", sep=';', dec=',', header=T)
      mira2clin = merge(mira2clin, mira2med, by = "id", all.x = TRUE)
      mira2clin = merge(mira2data[,c("id", "sex", "diet_group_datadriven", "age_blood")], mira2clin, by = "id", all.x = TRUE, all.y = FALSE)
      mira2clin <- merge(mira2clin, mira2master[,c("id","alcohol_use", "alcohol_use_frequency", "alcohol_doses_at_time", "alcohol_binge_drinking", "smoking", "physical_activity")], by = "id", all.x = TRUE, all.y = FALSE)
      mira2clin$alcohol_doses_at_time[mira2clin$alcohol_doses_at_time == 1] = 2
      mira2clin$alcohol_doses_at_time = mira2clin$alcohol_doses_at_time - 2
      mira2clin$alcohol_binge_drinking = mira2clin$alcohol_binge_drinking - 1
      mira2clin$alcohol_auditc = mira2clin$alcohol_use_frequency + mira2clin$alcohol_doses_at_time + mira2clin$alcohol_binge_drinking
      mira2clin$alcohol_auditc_threshold = FALSE
      mira2clin$alcohol_auditc_threshold[mira2clin$sex == "M" & mira2clin$age_blood > 10 & mira2clin$alcohol_auditc >= 6] = TRUE
      mira2clin$alcohol_auditc_threshold[mira2clin$sex == "F" & mira2clin$age_blood > 10 & mira2clin$alcohol_auditc >= 5] = TRUE
      rm(mira2med, mira2master)
    }
    # --------------------------------------------C. PLOTS------------------------------------------- #
    {
      #Figure 1A
      plotheightz_corr_w <- ggplot(data = mira2growth, aes(x=age_anthro, y=height_z_corr)) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#343434FF", "#2673ABFF", "#68A323FF"),
                            labels = c("Mixed diet", "Vegetarian", "Vegan"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(2,2,2),
                          labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        geom_line(aes(colour = diet_group_datadriven, group = id), size = 0.5, alpha = 0.3) +
        labs(y = "z-value", x = "Age (years)", title = "Corrected age-height") +
        geom_line(color="#343434FF",data = lm_pred_heightzcorr_mix_w, aes(x=age_anthro, y=height_z_corr), size = 1.4) +
        geom_line(color="#2673ABFF",data = lm_pred_heightzcorr_veg_w, aes(x=age_anthro, y=height_z_corr), size = 1.4) +
        geom_line(color="#68A323FF",data = lm_pred_heightzcorr_vgn_w, aes(x=age_anthro, y=height_z_corr), size = 1.4) +
        coord_cartesian(xlim = c(0,7), ylim = c(-4,4), expand = FALSE) +
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,20,0)),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              strip.text = element_text(colour = 'black'),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "top",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,-10,0),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()) +
        guides(linetype=F, size=F, color = guide_legend(override.aes = list(alpha = c(1,1,1))))
      filename = "plots/Article/Fig1A_height_z.svg"
      ggsave(file=filename, plot=plotheightz_corr_w, width=9, height=5, dpi = 300)
      # Output from Export --> SVG of size 1800 x 1000
      
      #Figure 1B
      plotpipa_w <- ggplot(data = mira2growth, aes(x=age_anthro, y=pipa)) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#343434FF", "#2673ABFF", "#68A323FF"),
                            labels = c("Mixed diet", "Vegetarian", "Vegan"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(2,2,2),
                          labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        geom_line(aes(colour = diet_group_datadriven, group = id), size = 0.5, alpha = 0.3) +
        labs(y = "%", x = "Age (years)", title = "Height-corrected weight") +
        geom_line(color="#343434FF",data = lm_pred_pipa_mix_w, aes(x=age_anthro, y=pipa), size = 1.4) +
        geom_line(color="#2673ABFF",data = lm_pred_pipa_veg_w, aes(x=age_anthro, y=pipa), size = 1.4) +
        geom_line(color="#68A323FF",data = lm_pred_pipa_vgn_w, aes(x=age_anthro, y=pipa), size = 1.4) +
        coord_cartesian(xlim = c(0,7), ylim = c(-25,40), expand = FALSE) +
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,20,0)),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              strip.text = element_text(colour = 'black'),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "top",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,-10,0),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()) +
        guides(linetype=F, size=F, color = guide_legend(override.aes = list(alpha = c(1,1,1))))
      filename = "plots/Article/Fig1B_Pipa.svg"
      ggsave(file=filename, plot=plotpipa_w, width=9, height=5, dpi = 300)
      # Output from Export --> SVG of size 1800 x 1000
      
      #Figure 1C
      plotheadcfz_w <- ggplot(data = mira2growth, aes(x=age_anthro, y=headcf_z)) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#343434FF", "#2673ABFF", "#68A323FF"),
                            labels = c("Mixed diet", "Vegetarian", "Vegan"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(2,2,2),
                          labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        geom_line(aes(colour = diet_group_datadriven, group = id), size = 0.5, alpha = 0.3) +
        labs(y = "z-value", x = "Age (years)", title = "Age-corrected head circumference") +
        geom_line(color="#343434FF",data = lm_pred_headcfz_mix_w, aes(x=age_anthro, y=headcfz), size = 1.4) +
        geom_line(color="#2673ABFF",data = lm_pred_headcfz_veg_w, aes(x=age_anthro, y=headcfz), size = 1.4) +
        geom_line(color="#68A323FF",data = lm_pred_headcfz_vgn_w, aes(x=age_anthro, y=headcfz), size = 1.4) +
        coord_cartesian(xlim = c(0,7), ylim = c(-3,3), expand = FALSE) +
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,20,0)),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              strip.text = element_text(colour = 'black'),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "top",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,-10,10),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()) +
        guides(linetype=F, size=F, color = guide_legend(override.aes = list(alpha = c(1,1,1))))
      filename = "plots/Article/Fig1C_headcf_z.svg"
      ggsave(file=filename, plot=plotheadcfz_w, width=6, height=4, dpi = 300)
      
      #Figure 1D
      plotbirthweight <- ggplot(data = mira2growth[mira2growth$age_anthro == 0 & !mira2growth$preterm,], aes(x=diet_pregnancy, y=weight/1000)) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#343434FF", "#2673ABFF", "#68A323FF"),
                            labels = c("Mixed", "Vegetarian", "Vegan"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Mixed", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Mixed", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Mixed", "Vegetarian", "Vegan")) +
        geom_boxplot(width = 2/3, staplewidth = 0.5, outliers = FALSE, coef = Inf, aes(colour = diet_pregnancy)) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = diet_pregnancy, shape = diet_pregnancy, size = diet_pregnancy)) +
        labs(y = "kg", x = "", title = "Birthweight\nby mother's diet during pregnancy") +
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,30,0)),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title.y = element_text(size = 12, face = "bold"),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              #strip.text = element_text(colour = 'black', size = 14, face = "bold"),
              strip.text = element_blank(),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.spacing = unit(1.5, "lines"),
              legend.position = "none",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,0,-20),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              #legend.background = element_rect(colour = "lightgray"),
              #legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12)) +
        facet_wrap(~sex, strip.position = "bottom") +
        guides(linetype=F, size=F)
      filename = "plots/Birthweight.svg"
      ggsave(file=filename, plot=plotbirthweight, width=5, height=4, dpi = 300)
      
      #Figure 2: All subfigures were generated with this piece of code, picking the variable in question as y and manually choosing ticked lines, if needed. P-values were manually assigned to the final subfigure, based on an extensive table of statistical test results (available in Supplementary Appendix)
      font_add_google("Source Sans Pro", "Source Sans Pro")
      plot_mira2_faceted <- ggplot(data = mira2data, aes(x=as.factor(diet_group_datadriven), y=FA_18.3n.3)) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#343434FF", "#2673ABFF", "#68A323FF"),    #4D6B73FF #009FC3FF #6CB43FFF  
                            labels = c("Mixed", "Vegetarian", "Vegan"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Mixed", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(2,2,2),
                          labels = c("Mixed", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Mixed", "Vegetarian", "Vegan")) +
        #scale_y_continuous(breaks = c(10, 20, 30, 40, 50, 100, 200)) + #for log-transformed Ferritin
        #coord_trans(y='log10') + #for log-transformed Ferritin
        geom_boxplot(width = 2/3, staplewidth = 0.5, outliers = FALSE, coef = Inf, aes(colour = diet_group_datadriven)) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = diet_group_datadriven, shape = diet_group_datadriven, size = diet_group_datadriven)) +
        #geom_hline(data=data.frame(cohort = factor(c("child", "adult, F", "adult, M")), cutoffs = c(3, 3, 3)), aes(yintercept = cutoffs), color = "darkgray", linetype = "dashed") + #for clinical cutoffs
        #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15, fun.args = list(conf.int = .95), size = 1.4) +
        #stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3, size = 1.4) +
        labs(y = "mol % of all fatty acids", x = "", title = "Serum ALA, 18:3n3") +
        #coord_cartesian(xlim = c(0.5,3.5), ylim = c(0,0.5), expand = FALSE) + #for tau/gly -conjugated bile acids
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,30,0)),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title.y = element_text(size = 12, face = "bold"),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              strip.text = element_text(colour = 'black', size = 12, face = "bold"),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.spacing = unit(1.5, "lines"),
              legend.position = "none",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,0,-20),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              #legend.background = element_rect(colour = "lightgray"),
              #legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12)) +
        facet_wrap(~cohort, strip.position = "bottom") +
        guides(linetype="none")
      #print to file
      filename = "plots/Article/Fig2D_ALA.svg"
      ggsave(file=filename, plot=plot_mira2_faceted, width=5, height=4, dpi = 300)
      
      #Figure 3A
      mira2_pathwaydata <- read.csv("mira2_pathwayanalysis_tophits.csv", sep = ';', dec = ',', header = T)
      mira2_pathwaydata$cohort <- factor(mira2_pathwaydata$cohort, levels = c("children", "adult, F", "adult, M"))
      plotpathwayanalysis <- ggplot(data = mira2_pathwaydata, aes(x=impact, y=p.logneg)) +
        theme_light() +
        #scale_colour_manual(name = "Cohort",
        #                    values = c("grey10", "grey30", "grey50"),
        #                    labels = c("Adult, F", "Adult, M", "Children"),
        #                    guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Cohort", values=c(17,15,19),
                           labels = c("Children", "Adult, F", "Adult, M")) +
        #scale_size_manual(name = "Dietary group", values=c(2,2,2),
        #                  labels = c("Mixed diet", "Vegetarian", "Vegan")) +
        geom_point(aes(size = 4, shape = cohort), colour = "grey10") +
        labs(y = "p-value (-log10 scale)", x = "Pathway impact", title = "") +
        coord_cartesian(xlim = c(-0.01,1.1), ylim = c(0,2.5), expand = FALSE) +
        geom_hline(yintercept = -log10(0.05), color = "darkgray", linetype = "dashed") +
        geom_vline(xintercept = 0.1, color = "darkgray", linetype = "dashed") +
        theme(text = element_text(family = "Source Sans Pro"),
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(0,0,20,0)),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = 'white'),
              strip.text = element_text(colour = 'black'),
              strip.placement = "outside",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "top",
              legend.justification = "left",
              legend.direction = "horizontal",
              legend.box.margin = margin(-20,0,-10,0),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(0.4, "cm"),
              legend.title = element_blank()) +
        guides(linetype=F, size = F, color = F)
      filename = "plots/Article/Fig3A_Pathwayanalysis.svg"
      ggsave(file=filename, plot=plotpathwayanalysis, width=4, height=4, dpi = 300)
    }
  }

#------END OF FILE------#