
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(MatchIt)
library(viridis)
library(mice)
library(cetcolor)

#region Load data and merge

# MRI modality availability
df_mri_avail = as.data.frame(fread("../../UKB/QC/mri_file_availability.csv"))
df_mri_avail = df_mri_avail %>%
    mutate(ID = as.numeric(substr(ID,5,11))) %>%
    select(ID, t1w_ses2, flair_ses2, dMRI_ses2) %>%
    glimpse()

# Demo
df_demo = as.data.frame(fread("../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
df_demo = df_demo %>%
    filter(InstanceID == 2 & is.na(Age_when_attended_assessment_centre_21003_0) == FALSE) %>%
    rename(ID = "SubjectID", Sex = "Sex_31_0") %>%
    select(ID, Sex) %>%
    glimpse()

# Menopause etc.
df_menopause = as.data.frame(fread("../../UKB/tabular/df_menopause_denise/UKBB_menopause_denise_wider.tsv"))
df_menopause = df_menopause %>%
    rename(
        ID = "SubjectID",
        Waist_circumference = "Waist_circumference_48_0",
        Hip_circumference = "Hip_circumference_49_0",
        Income = "Average_total_household_income_before_tax_738_0",
        Days_walked = "Number_of_days/week_walked_10+_minutes_864_0",
        Duration_walks = "Duration_of_walks_874_0",
        Days_moderate_activity = "Number_of_days/week_of_moderate_physical_activity_10+_minutes_884_0",
        Duration_moderate_activity = "Duration_of_moderate_activity_894_0",
        Days_vigorous_activity = "Number_of_days/week_of_vigorous_physical_activity_10+_minutes_904_0",
        Duration_vigorous_activity = "Duration_of_vigorous_activity_914_0",
        Past_smoking = "Past_tobacco_smoking_1249_0",
        Alcohol_frequency = "Alcohol_intake_frequency._1558_0",
        Had_menopause = "Had_menopause_2724_0",
        MHT_used = "Ever_used_hormone-replacement_therapy_(HRT)_2814_0",
        Oophorectomy = "Bilateral_oophorectomy_(both_ovaries_removed)_2834_0",
        Age_menopause = "Age_at_menopause_(last_menstrual_period)_3581_0",
        Diastolic_blood_pressure_1 = "Diastolic_blood_pressure,_automated_reading_4079_0",
        Diastolic_blood_pressure_2 = "Diastolic_blood_pressure,_automated_reading_4079_1",
        Systolic_blood_pressure_1 = "Systolic_blood_pressure,_automated_reading_4080_0",
        Systolic_blood_pressure_2 = "Systolic_blood_pressure,_automated_reading_4080_1",
        BP_med_1 = "Medication_for_cholesterol,_blood_pressure,_diabetes,_or_take_exogenous_hormones_6153_0",
        BP_med_2 = "Medication_for_cholesterol,_blood_pressure,_diabetes,_or_take_exogenous_hormones_6153_1",
        BP_med_3 = "Medication_for_cholesterol,_blood_pressure,_diabetes,_or_take_exogenous_hormones_6153_2",
        BP_med_4 = "Medication_for_cholesterol,_blood_pressure,_diabetes,_or_take_exogenous_hormones_6153_3",
        Smoking_status = "Smoking_status_20116_0",
        Alcohol_status = "Alcohol_drinker_status_20117_0",
        Ever_smoked = "Ever_smoked_20160_0",
        Age = "Age_when_attended_assessment_centre_21003_0",
        BMI = "Body_mass_index_(BMI)_23104_0",
        Pack_years_smoking = "Pack_years_of_smoking_20161_0",
        MHT_age_ended = "Age_last_used_hormone-replacement_therapy_(HRT)_3546_0",
        MHT_age_started = "Age_started_hormone-replacement_therapy_(HRT)_3536_0",
        Age_high_BP = "Age_high_blood_pressure_diagnosed_2966_0",
        Age_oophorectomy = "Age_at_bilateral_oophorectomy_(both_ovaries_removed)_3882_0"
    ) %>%
    select(-c("InstanceID", "Diabetes_diagnosed_by_doctor_2443_0", "Systolic_blood_pressure,_manual_reading_93_0", "Systolic_blood_pressure,_manual_reading_93_1", "Diastolic_blood_pressure,_manual_reading_94_0", "Diastolic_blood_pressure,_manual_reading_94_1")) %>%
    glimpse()

# MRI site
df_site = as.data.frame(fread("../../UKB/tabular/df_mri_site/UKBB_mri_site_wide.tsv"))
df_site = df_site %>%
    filter(InstanceID == 2) %>%
    rename(ID = "SubjectID", site="UK Biobank assessment centre_54") %>%
    select(ID, site) %>%
    glimpse()

# MRI QC
df_mri_qc = as.data.frame(fread("../../UKB/QC/all_qc_merged_no_dwi_final.csv"))
df_dwi_qc = as.data.frame(fread("../../UKB/tabular/dMRI_ASL_derived/ukb_tab_dwi.csv"))
df_dwi_qc = subset(df_dwi_qc, select=c("eid", "24450-2.0", "24453-2.0"))
colnames(df_dwi_qc) = c("ID", "abs_motion", "rel_motion")
df_dwi_qc = df_dwi_qc[complete.cases(df_dwi_qc),]
df_qc_surg = as.data.frame(fread("../data/surgical_qc_copy.csv"))

df_qc_surg = df_qc_surg %>%
    rename(ID = "V1", t1w_surg = "V2", civet_surg = "V3") %>%
    mutate(ID = as.numeric(substr(ID, 6, 12))) %>%
    mutate(
        t1w_surg = pmin(2, max(t1w_surg, na.rm = TRUE) - t1w_surg),
        civet_surg = pmin(2, max(civet_surg, na.rm = TRUE) - civet_surg)
    ) %>%
    glimpse()

df_mri_qc = df_mri_qc %>%
    mutate(ID = as.numeric(substr(ID,5,11))) %>%
    filter(ses == 2) %>%
    distinct(ID, .keep_all = TRUE) %>%
    full_join(df_qc_surg, by="ID") %>%
    left_join(df_dwi_qc, by="ID") %>%
    mutate(
        # Structural pass: 2 on anybody's rating
        t1_fail = factor(if_else(
            rowSums(across(c(t1w_motion_olivier, t1w_motion_manuela, t1w_motion_grace, t1w_motion_daniela, t1w_surg)) == 2, na.rm = TRUE) > 0,
            "Pass", "Fail")),
        flair_fail = factor(if_else(
            rowSums(across(c(flair_motion_olivier, flair_motion_daniela)) == 2, na.rm = TRUE) > 0,
            "Pass", "Fail")),
        # # Structural fail: 0 on anybody's rating OR NA on everybody
        # t1_fail = factor(if_else(
        #     rowSums(across(c(t1w_motion_olivier, t1w_motion_manuela, t1w_motion_grace, t1w_motion_daniela, t1w_surg)) == 0, na.rm = TRUE) > 0,
        #     "Fail", "Pass")),
        # flair_fail = factor(if_else(
        #     rowSums(across(c(flair_motion_olivier, flair_motion_daniela)) == 0, na.rm = TRUE) > 0,
        #     "Fail", "Pass")),
        # DWI fail: abs_motion > 3
        dwi_fail = factor(if_else(abs_motion > 3, "Fail", "Pass"))
    ) %>%
    select(ID, t1_fail, flair_fail, dwi_fail) %>%
    glimpse()

# Bison volumes
df_bison = as.data.frame(fread("../../UKB/BISON/BISON_volumes.csv"))
df_bison = df_bison %>%
    filter(ses == 2) %>%
    # TBV variable as sum of all brain tissues
    mutate(TBV = Cerebellum_GM + Cerebellum_WM + Brainstem + Subcortical_GM + Cortical_GM + Cerebral_WM + WMH) %>%
    select(ID, WMH, TBV) %>%
    glimpse()

# UKB-derived brain variables (WMH, TBV, WM volume)
df_brain_ukb = as.data.frame(fread("../../UKB/Analyses/clean_brain_tab/results/brain_tab_clean_wide.tsv"))
df_brain_ukb = df_brain_ukb %>% select(
    c("ID", 
    "Total volume of white matter hyperintensities (from T1 and T2_FLAIR images)",
    "Volume of brain, grey+white matter",
    "Volume of white matter"
    )) %>%
    rename(
        UKB_WMH = "Total volume of white matter hyperintensities (from T1 and T2_FLAIR images)",
        UKB_TBV = "Volume of brain, grey+white matter",
        UKB_WM = "Volume of white matter"
    ) %>%
    glimpse()

# Diagnoses
df_diagnoses = as.data.frame(fread("../../UKB/Analyses/clean_firstocc/results/firstocc_categ.tsv"))
df_diagnoses = df_diagnoses %>%
    # Only keep stroke, MS, diabetes
    filter(icd_code %in% c("I60", "I61", "I63", "I64", "G35", "E10", "E11", "E13", "E14")) %>%
    glimpse()

# Genetic sex
df_genetic_sex = as.data.frame(fread("../../UKB/tabular/df_genetic_sex/UKBB_genetic_sex_wide.tsv"))
df_genetic_sex = df_genetic_sex %>%
    filter(InstanceID == 2) %>%
    rename(ID = "SubjectID", Genetic_sex = "Genetic sex_22001") %>%
    select(ID, Genetic_sex) %>%
    glimpse()

#### Merge ####
df_merge = df_mri_avail %>%
    left_join(df_demo, by="ID") %>%
    left_join(df_genetic_sex, by="ID") %>%
    left_join(df_site, by="ID") %>%
    left_join(df_mri_qc, by="ID") %>%
    left_join(df_menopause, by="ID") %>%
    left_join(df_bison, by="ID") %>%
    left_join(df_brain_ukb, by="ID") %>%
    glimpse()

#endregion

#region Clean and transform data

ids_diabetes = df_diagnoses %>% filter(icd_code %in% c("E10", "E11", "E13", "E14") & days_mri_dx < 0) %>% pull(ID)

df_clean = df_merge %>%
    # Transformations of lifestyle variables
    mutate(
        Alcohol_status = if_else(is.na(Alcohol_status), "NA", if_else(Alcohol_status %in% c("Never", "Previous"), "No", "Yes")),
        Alcohol_frequency = as.numeric(factor(Alcohol_frequency, levels=c("Never", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week", "Daily or almost daily"), labels=c(1,2,3,4,5,6))),
        BP_med = if_else(rowSums(across(c(BP_med_1, BP_med_2, BP_med_3, BP_med_4)) == "Blood pressure medication", na.rm = TRUE) > 0, "Yes", "No"),
        Avg_diastolic_BP = (Diastolic_blood_pressure_1 + Diastolic_blood_pressure_2) / 2,
        Avg_systolic_BP = (Systolic_blood_pressure_1 + Systolic_blood_pressure_2) / 2,
        Diabetes = if_else(ID %in% ids_diabetes, "Yes", "No"),
        WHR = Waist_circumference / Hip_circumference,
        Duration_walks = log1p(pmin(Duration_walks, mean(Duration_walks, na.rm = TRUE) + (3.5 * sd(Duration_walks, na.rm = TRUE)))),
        Duration_moderate_activity = log1p(pmin(Duration_moderate_activity, mean(Duration_moderate_activity, na.rm = TRUE) + (3.5 * sd(Duration_moderate_activity, na.rm = TRUE)))),
        Duration_vigorous_activity = log1p(pmin(Duration_vigorous_activity, mean(Duration_vigorous_activity, na.rm = TRUE) + (3.5 * sd(Duration_vigorous_activity, na.rm = TRUE)))),
        Smoking_status = if_else(is.na(Smoking_status), "NA", if_else(Smoking_status %in% c("Never", "Previous"), "No", "Yes")),
        Past_smoking = as.numeric(factor(Past_smoking, levels=c("I have never smoked", "Just tried once or twice", "Smoked occasionally", "Smoked on most or all days"), labels=c(1,2,3,4))),
        Pack_years_smoking = log1p(if_else(Ever_smoked == "No", 0, Pack_years_smoking)),
        High_BP = if_else(is.na(Age_high_BP), "No", "Yes")
    ) %>%
    select(-c(BP_med_1, BP_med_2, BP_med_3, BP_med_4, Diastolic_blood_pressure_1, Diastolic_blood_pressure_2, Systolic_blood_pressure_1, Systolic_blood_pressure_2)) %>%
    # Transformations of WMH variables
    mutate(
        WMH_divTBV_log = log(WMH/TBV),
        UKB_WMH_divTBV_log = log(UKB_WMH/UKB_TBV),
        UKB_WMH_divWM_log = log(UKB_WMH/UKB_WM)
    ) %>%
    # Transformations of menopause variables
    mutate(
        Hysterectomy = if_else(is.na(Had_menopause), "NA", if_else(Had_menopause == "Not sure - had a hysterectomy", "Yes", "No")),
        Had_menopause = if_else(Had_menopause == "" | Had_menopause == "Not sure - other reason" | Had_menopause == "Not sure - had a hysterectomy", "NA", as.character(Had_menopause)),
        Oophorectomy = if_else(Oophorectomy == "Not sure", "NA", Oophorectomy),
        hyst_no_ooph = if_else(Hysterectomy == "Yes" & Oophorectomy == "No", "Yes", "No"),
        MHT_age_ended = as.numeric(if_else(MHT_age_ended == "Still taking MHT", as.character(Age), as.character(MHT_age_ended))),
        MHT_years = log1p(if_else(MHT_used == "No", 0, MHT_age_ended - MHT_age_started))
    ) %>%
    # Fix NA strings
    mutate(across(everything(), ~na_if(.x, "NA")), across(everything(), ~na_if(.x, ""))) %>%
    # Make menopause groups
    mutate(
        # Always exclude Hysterectomy without Oophorectomy
        Menopause_group = factor(case_when(
            # PRE = No menopause
            Had_menopause == "No" & Oophorectomy == "No" ~ "PRE",
            # POST = Menopause and no oophorectomy
            Had_menopause == "Yes" & Oophorectomy == "No" ~ "POST",
            # POST = Menopause, oophorectomy after menopause
            Had_menopause == "Yes" & Oophorectomy == "Yes" & Age_menopause < Age_oophorectomy ~ "POST",
            # POST = Menopause, no record of oophorectomy
            Had_menopause == "Yes" & is.na(Oophorectomy) ~ "POST",
            # SURG = Menopause, oophorectomy before menopause
            Had_menopause == "Yes" & Oophorectomy == "Yes" & Age_menopause >= Age_oophorectomy ~ "SURG",
            # SURG = No menopause, oophorectomy
            Had_menopause == "No" & Oophorectomy == "Yes" ~ "SURG",
            # SURG = Menopause, oophorectomy, no data on age at menopause
            Had_menopause == "Yes" & Oophorectomy == "Yes" & is.na(Age_menopause) ~ "SURG",
            # SURG = No data on menopause, oophorectomy
            is.na(Had_menopause) & Oophorectomy == "Yes" ~ "SURG"
        ))
    ) %>%
    # Lohner combine POST and SURG
    mutate(
        Age_menopause = if_else(Oophorectomy == "Yes", Age_oophorectomy, Age_menopause),
        Time_since_menopause = Age - Age_menopause,
        lohner_Menopause_group = factor(if_else(Menopause_group == "SURG", "POST", as.character(Menopause_group)), levels=c("PRE", "POST"))
    ) %>%
    # Change variable types
    mutate(
        Sex = factor(Sex),
        Genetic_sex = factor(Genetic_sex),
        site = factor(site),
        Income = as.numeric(factor(Income, levels=c("Less than 18,000", "18,000 to 30,999", "31,000 to 51,999", "52,000 to 100,000", "Greater than 100,000"), labels=c(1,2,3,4,5))),
        Days_walked = as.numeric(Days_walked),
        Had_menopause = factor(Had_menopause),
        MHT_used = factor(MHT_used),
        Oophorectomy = factor(Oophorectomy),
        Smoking_status = factor(Smoking_status),
        Alcohol_status = factor(Alcohol_status),
        Ever_smoked = factor(Ever_smoked),
        MHT_age_ended = as.numeric(MHT_age_ended),
        BP_med = factor(BP_med),
        Diabetes = factor(Diabetes),
        High_BP = factor(High_BP),
        Hysterectomy = factor(Hysterectomy),
        hyst_no_ooph = factor(hyst_no_ooph),
        Menopause_group = factor(Menopause_group, levels=c("PRE", "POST", "SURG"))
    ) %>%
    glimpse()

#endregion

#region Exclusions

ids_diabetes = df_diagnoses %>% filter(icd_code %in% c("E10", "E11", "E13", "E14") & days_mri_dx < 0) %>% pull(ID)
ids_stroke = df_diagnoses %>% filter(icd_code %in% c("I60", "I61", "I63", "I64") & days_mri_dx < 0) %>% pull(ID)
ids_ms = df_diagnoses %>% filter(icd_code %in% c("G35")) %>% pull(ID)

df_final = df_clean %>%
    # Remove all men
    filter(Sex == "Female") %T>%
    {print(paste("Remove all men: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Remove missing T1w, FLAIR, or DWI
    filter(t1w_ses2 == 1 & flair_ses2 == 1 & dMRI_ses2 == 1) %T>%
    {print(paste("Remove missing T1w, FLAIR, or DWI: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Remove motion in T1w
    filter(t1_fail == "Pass" & is.na(WMH_divTBV_log) == FALSE & is.na(UKB_WMH_divTBV_log) == FALSE) %T>%
    {print(paste("Remove motion in T1w: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Diagnosis of stroke or MS
    filter(!ID %in% ids_stroke) %T>%
    filter(!ID %in% ids_ms) %T>%
    {print(paste("Diagnosis of stroke or MS: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Hysterectomy without oophorectomy
    filter(hyst_no_ooph == "No") %T>%
    {print(paste("Hysterectomy without oophorectomys: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Unmacthed menopause group
    filter(is.na(Menopause_group) == FALSE) %T>%
    {print(paste("Unmacthed menopause group: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Missing data on income, MHT use
    filter(is.na(Income) == FALSE & is.na(MHT_used) == FALSE) %T>%
    {print(paste("Missing data on income, age at menopause, MHT use: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Missing data on age at menopause if POST or SURF
    filter(!(Menopause_group %in% c("SURG", "POST") & is.na(Age_menopause))) %T>%
    {print(paste("Missing data on age at menopause if POST or SURF: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # # Pre-menopausal reporting MHT use
    filter(!(Menopause_group %in% c("PRE") & MHT_used == "Yes")) %T>%
    {print(paste("Pre-menopausal reporting MHT use: n = ", nrow(.))); print(summary(.$Menopause_group))} %>%
    # Women with male genetic sex
    filter(Genetic_sex != "Male" | is.na(Genetic_sex)) %T>%
    {print(paste("Women with male genetic sex: n = ", nrow(.))); print(summary(.$Menopause_group))}

#endregion

#region Age match groups

colors = viridis(3, option = "D")
color_scale_groups = c(PRE = colors[1], POST = colors[2], SURG = colors[3])

# Function to generate summary tables (Table 1)
table_df = function(df_func, output) {
    df_table = df_func %>%
        group_by(Menopause_group) %>%
        summarise(
            Sample_size = n(),
            Age_mean = mean(Age, na.rm = TRUE),
            Age_sd = sd(Age, na.rm = TRUE),
            Age_min = min(Age, na.rm = TRUE),
            Age_max = max(Age, na.rm = TRUE),
            Age_menopause_mean = mean(Age_menopause, na.rm = TRUE),
            Age_menopause_sd = sd(Age_menopause, na.rm = TRUE),
            Time_since_menopause_mean = mean(Time_since_menopause, na.rm = TRUE),
            Time_since_menopause_sd = sd(Time_since_menopause, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        left_join(
            df_func %>%
            count(Menopause_group, Income) %>%
            pivot_wider(names_from = Income, values_from = n, values_fill = 0) %>%
            mutate(Income_counts = pmap_chr(across(-Menopause_group), ~toString(na.omit(c(...))))) %>%
            select(Menopause_group, Income_counts),
            by = "Menopause_group"
        ) %>%
        left_join(
            df_func %>%
            count(Menopause_group, MHT_used) %>%
            pivot_wider(names_from = MHT_used, values_from = n, values_fill = 0) %>%
            mutate(MHT_counts = pmap_chr(across(-Menopause_group), ~toString(na.omit(c(...))))) %>%
            select(Menopause_group, MHT_counts),
            by = "Menopause_group"
        )

    fwrite(df_table, output, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    print(output)
}

# Function to generate age and WMH density plots (Figure 2BC)
age_wmh_density = function(df_func, output, bounds_age) {
    plt_age = ggplot(df_func, aes(x=Age, color=Menopause_group)) +
        geom_density(size=2, adjust=1.2) +
        theme_classic() +
        scale_color_manual(values = color_scale_groups) +
        scale_y_continuous(name="Density") +
        scale_x_continuous(name="Age") +
        guides(color = "none") + 
        theme(text=element_text(size=20))
    ggsave(paste0(output, "_age.png"), height=5, width=6)
    print(paste0(output, "_age.png"))

    plt_wmh_bison = ggplot(df_func, aes(x=WMH_divTBV_log, color=Menopause_group)) +
        geom_density(size=2, adjust=1.2) +
        theme_classic() +
        scale_color_manual(values = color_scale_groups) +
        scale_y_continuous(name="Density") +
        scale_x_continuous(name="log(WMH/TBV)") +
        guides(color = "none") + 
        theme(text=element_text(size=20))
    ggsave(paste0(output, "_WMH_bison.png"), height=5, width=6)

    plt_wmh_UKB = ggplot(df_func, aes(x=UKB_WMH_divTBV_log, color=Menopause_group)) +
        geom_density(size=2, adjust=1.2) +
        theme_classic() +
        scale_color_manual(values = color_scale_groups) +
        scale_y_continuous(name="Density") +
        scale_x_continuous(name="log(WMH/TBV)") +
        guides(color = "none") + 
        theme(text=element_text(size=20))
    ggsave(paste0(output, "_WMH_UKB.png"), height=5, width=6)

    plt_wmh_UKB = ggplot(df_func, aes(x=UKB_WMH_divWM_log, color=Menopause_group)) +
        geom_density(size=2, adjust=1.2) +
        theme_classic() +
        scale_color_manual(values = color_scale_groups) +
        scale_y_continuous(name="Density") +
        scale_x_continuous(name="log(WMH/WMV)") +
        guides(color = "none") + 
        theme(text=element_text(size=20))
    ggsave(paste0(output, "_WMH_UKB_divWM.png"), height=5, width=6)
}

############ A: Full unmatched sample
df_A = df_final 

table_df(df_A, "./results/df_A_summary.tsv")
age_wmh_density(df_A, "./visualization/density_A", c(min(df_A$Age), max(df_A$Age)))
summary(df_A$Menopause_group)

############ B: Nearest-neighbour age-matched sample
set.seed(123)

# Match PRE and SURG (remove >59 years old)
df_tmp = df_final %>% 
    filter(Age < 60 & Menopause_group %in% c("SURG", "PRE")) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("SURG", "PRE")))

match_B_pre_surg = matchit(Menopause_group ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=2, replace = FALSE)
df_match_B_pre_surg = match.data(match_B_pre_surg) %>%
    select(-c(distance, weights, subclass)) %>%
    bind_rows(df_final %>% filter(Menopause_group == "POST")) %>%
    mutate(PRESURG_POST = factor(if_else(Menopause_group %in% c("SURG", "PRE"), "PRESURG", "POST")))

# Match POST and PRESURG (remove >59 years old)
df_tmp = df_match_B_pre_surg %>%
    filter(Age < 60) %>%
    mutate(PRESURG_POST = factor(PRESURG_POST, levels=c("PRESURG", "POST")))

match_B_post = matchit(PRESURG_POST ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=0, replace = FALSE)
df_B = match.data(match_B_post) %>%
    select(-c(distance, weights, subclass)) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("PRE", "POST", "SURG")))

summary(df_B$Menopause_group)

table_df(df_B, "./results/df_B_summary.tsv")
age_wmh_density(df_B, "./visualization/density_B")

############ C: Exact age-matched sample
set.seed(123)

# Match PRE and SURG (remove >59 years old)
df_tmp = df_final %>% 
    filter(Age < 60 & Menopause_group %in% c("PRE", "SURG")) %>% 
    mutate(Menopause_group = factor(Menopause_group, levels=c("PRE", "SURG")))

match_C_pre_surg = matchit(Menopause_group ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=0, replace = FALSE, exact = "Age")
df_match_C_pre_surg = match.data(match_C_pre_surg) %>%
    select(-c(distance, weights, subclass)) %>%
    bind_rows(df_final %>% filter(Menopause_group == "POST")) %>%
    mutate(PRESURG_POST = factor(if_else(Menopause_group %in% c("SURG", "PRE"), "PRESURG", "POST")))

# Match POST and PRESURG (remove >59 years old)
df_tmp = df_match_C_pre_surg %>%
    filter(Age < 60) %>%
    mutate(PRESURG_POST = factor(PRESURG_POST, levels=c("PRESURG", "POST")))

match_C_post = matchit(PRESURG_POST ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=0, replace = FALSE, exact = "Age")
df_C = match.data(match_C_post) %>%
    select(-c(distance, weights, subclass)) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("PRE", "POST", "SURG")))

summary(df_C$Menopause_group)

table_df(df_C, "./results/df_C_summary.tsv")
age_wmh_density(df_C, "./visualization/density_C")

############ D: POST-SURG age-matched sample
set.seed(123)

df_tmp = df_final %>% 
    filter(Menopause_group %in% c("POST", "SURG")) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("POST", "SURG")))

match_D_post_surg = matchit(Menopause_group ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=0, replace = FALSE, exact = "Age")
df_D = match.data(match_D_post_surg) %>%
    select(-c(distance, weights, subclass)) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("POST", "SURG")))

summary(df_D$Menopause_group)

table_df(df_D, "./results/df_D_summary.tsv")
age_wmh_density(df_D, "./visualization/density_D")

############ E: Lohner et al. (2022) replicated sample
set.seed(123)

df_E = df_final %>%
    filter(Age < 60) %>%
    mutate(Menopause_group = factor(if_else(Menopause_group == "SURG", "POST", as.character(Menopause_group)), levels=c("PRE", "POST")))

summary(df_E$Menopause_group)

table_df(df_E, "./results/df_E_summary.tsv")
age_wmh_density(df_E, "./visualization/density_E")

############ F: >5 years after menopause nearest neighbour age-matched sample
set.seed(123)

# Match PRE and SURG
df_tmp = df_final %>% 
    filter(Age < 60 & Menopause_group %in% c("SURG", "PRE") & (Time_since_menopause > 5 | is.na(Time_since_menopause))) %>% 
    mutate(Menopause_group = factor(Menopause_group, levels=c("SURG", "PRE")))

match_F_pre_surg = matchit(Menopause_group ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=2, replace = FALSE)
df_match_F_pre_surg = match.data(match_F_pre_surg) %>%
    select(-c(distance, weights, subclass)) %>%
    bind_rows(df_final %>% filter(Menopause_group == "POST")) %>%
    mutate(PRESURG_POST = factor(if_else(Menopause_group %in% c("SURG", "PRE"), "PRESURG", "POST")))

# Match POST and PRESURG
df_tmp = df_match_F_pre_surg %>% 
    filter(Age < 60 & Time_since_menopause > 5 | is.na(Time_since_menopause)) %>% 
    mutate(PRESURG_POST = factor(PRESURG_POST, levels=c("PRESURG", "POST")))

match_F_post = matchit(PRESURG_POST ~ Age, data=df_tmp, method = "nearest", distance = "glm", m.order = "random", caliper=0, replace = FALSE)
df_F = match.data(match_F_post) %>%
    select(-c(distance, weights, subclass)) %>%
    mutate(Menopause_group = factor(Menopause_group, levels=c("PRE", "POST", "SURG")))

summary(df_F$Menopause_group)

table_df(df_F, "./results/df_F_summary.tsv")
age_wmh_density(df_F, "./visualization/density_F")

############ G: Extra: full unmatched, same IDs as Denise's

df_denise = as.data.frame(fread("../data/data_full.csv"))
df_G = df_clean %>%
    full_join(
        df_denise %>%
            select(ID, group) %>%
            rename(Menopause_group_denise = "group"),
        by="ID"
    ) %>%
    mutate(
        Menopause_group_denise = factor(Menopause_group_denise, levels=c("premen", "postmen", "surgical"), labels=c("PRE", "POST", "SURG"))
    ) %>%
    filter(!is.na(Menopause_group) & !is.na(Menopause_group_denise)) %>%
    glimpse()

summary(df_G$Menopause_group)
summary(df_G$Menopause_group_denise)

table_df(df_G, "./results/df_G_summary.tsv")
age_wmh_density(df_G, "./visualization/density_G")

#endregion

#region Impute cardiometabolic factors

factor_list = c(
    "Alcohol_frequency",
    "Alcohol_status",
    "Avg_diastolic_BP",
    "Avg_systolic_BP",
    "High_BP",
    "Diabetes",
    "BMI",
    "WHR",
    "Days_moderate_activity",
    "Duration_moderate_activity",
    "Days_vigorous_activity",
    "Duration_vigorous_activity",
    "Days_walked",
    "Duration_walks",
    "Pack_years_smoking",
    "Ever_smoked",
    "Smoking_status"
)

# Impute full unmatched sample
meth = make.method(df_A)
meth[!(names(meth) %in% factor_list)] = ""

pred = make.predictorMatrix(df_A)
pred[,] = 0
pred[factor_list, c(factor_list, "Menopause_group")] = 0 # Don't use menopause group for imputation

meth["Menopause_group"] = ""
pred["Menopause_group", ] = 0 # Don't impute menopause group

df_A_imputed = mice(df_A, m = 5, method = meth, predictorMatrix = pred, seed = 123)

#endregion

#region Differences in WMHV across the menopausal transition

# Make list of cardiometabolic factors to add in lm formulas

factor_list = c(
    "Alcohol_frequency",
    "Alcohol_status",
    "Avg_diastolic_BP",
    "Avg_systolic_BP",
    "High_BP",
    "Diabetes",
    "BMI",
    "WHR",
    "Days_moderate_activity",
    "Duration_moderate_activity",
    "Days_vigorous_activity",
    "Duration_vigorous_activity",
    "Days_walked",
    "Duration_walks",
    "Pack_years_smoking",
    "Ever_smoked",
    "Smoking_status"
)

cardio_factors <- paste(factor_list, collapse = " + ")

# Init results list
results = list()

############ Sample A: unmatched

# Poly age
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_A))
results[["sampleA_polyAge"]] = data.frame(
        Name = "sampleA_polyAge",
        Sample = "A",
        coef_POST = lm_result$coefficients["Menopause_groupPOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_groupPOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleA_polyAge"]]))

# Linear age
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + Age + Income + MHT_used + site, data=df_A))
results[["sampleA_linAge"]] = data.frame(
        Name = "sampleA_linAge",
        Sample = "A",
        coef_POST = lm_result$coefficients["Menopause_groupPOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_groupPOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleA_linAge"]]))

# Poly age + cov cardiometabolic factors
lm_result = summary(pool(with(df_A_imputed, lm(as.formula(paste0("scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site +", cardio_factors)), subset = ID %in% df_A$ID))))
results[["sampleA_polyAge_covCardio"]] = data.frame(
        Name = "sampleA_polyAge_covCardio",
        Sample = "A",
        coef_POST = lm_result$estimate[2],
        coef_SURG = lm_result$estimate[3],
        pval_POST = lm_result$p.value[2],
        pval_SURG = lm_result$p.value[3]
    )

print(t(results[["sampleA_polyAge_covCardio"]]))

############ Sample B: Matched (NN)

# Main
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_B))
results[["sampleB"]] = data.frame(
        Name = "sampleB",
        Sample = "B",
        coef_POST = lm_result$coefficients["Menopause_groupPOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_groupPOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleB"]]))

# Cov cardiometabolic factors
lm_result = summary(pool(with(df_A_imputed, lm(as.formula(paste0("scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site +", cardio_factors)), subset = ID %in% df_B$ID))))
results[["sampleB_covCardio"]] = data.frame(
        Name = "sampleB_covCardio",
        Sample = "B",
        coef_POST = lm_result$estimate[2],
        coef_SURG = lm_result$estimate[3],
        pval_POST = lm_result$p.value[2],
        pval_SURG = lm_result$p.value[3]
    )

print(t(results[["sampleB_covCardio"]]))

############ Sample C: Matched (exact)

# Main
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_C))
results[["sampleC"]] = data.frame(
        Name = "sampleC",
        Sample = "C",
        coef_POST = lm_result$coefficients["Menopause_groupPOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_groupPOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleC"]]))

# Cov cardiometabolic factors
lm_result = summary(pool(with(df_A_imputed, lm(as.formula(paste0("scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site +", cardio_factors)), subset = ID %in% df_C$ID))))
results[["sampleC_covCardio"]] = data.frame(
        Name = "sampleC_covCardio",
        Sample = "C",
        coef_POST = lm_result$estimate[2],
        coef_SURG = lm_result$estimate[3],
        pval_POST = lm_result$p.value[2],
        pval_SURG = lm_result$p.value[3]
    )

print(t(results[["sampleC_covCardio"]]))

############ Sample D: POST-SURG

# Main
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_D))
results[["sampleD"]] = data.frame(
        Name = "sampleD",
        Sample = "D",
        coef_POST = NA,
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = NA,
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleD"]]))

# Cov cardiometabolic factors
lm_result = summary(pool(with(df_A_imputed, lm(as.formula(paste0("scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site +", cardio_factors)), subset = ID %in% df_D$ID))))
results[["sampleD_covCardio"]] = data.frame(
        Name = "sampleD_covCardio",
        Sample = "D",
        coef_POST = NA,
        coef_SURG = lm_result$estimate[2],
        pval_POST = NA,
        pval_SURG = lm_result$p.value[2]
    )

print(t(results[["sampleD_covCardio"]]))

############ Sample E: Lohner (UKB / WM WMH volumes, linear age cov, and add cardiovascular covariates)

# Main
lm_result = summary(pool(with(df_A_imputed, lm(scale(UKB_WMH_divWM_log) ~ lohner_Menopause_group + scale(Age) + Income + site + High_BP + Smoking_status + Ever_smoked + Pack_years_smoking + BMI, subset = ID %in% df_E$ID))))
results[["sampleE_covCardio_UKBWMH"]] = data.frame(
        Name = "sampleE_covCardio_UKBWMH",
        Sample = "E",
        coef_POST = lm_result$estimate[2],
        coef_SURG = NA,
        pval_POST = lm_result$p.value[2],
        pval_SURG = NA
    )

print(t(results[["sampleE_covCardio_UKBWMH"]]))

############ Sample F: Time since menopause > 5 years

# Main
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_F))
results[["sampleF"]] = data.frame(
        Name = "sampleF",
        Sample = "F",
        coef_POST = lm_result$coefficients["Menopause_groupPOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_groupPOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleF"]]))

# Cov cardiometabolic factors
lm_result = summary(pool(with(df_A_imputed, lm(as.formula(paste0("scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site +", cardio_factors)), subset = ID %in% df_F$ID))))
results[["sampleF_covCardio"]] = data.frame(
        Name = "sampleF_covCardio",
        Sample = "F",
        coef_POST = lm_result$estimate[2],
        coef_SURG = lm_result$estimate[3],
        pval_POST = lm_result$p.value[2],
        pval_SURG = lm_result$p.value[3]
    )

print(t(results[["sampleF_covCardio"]]))

############ Sample G: Same as Denise's unmatched sample from paper

# Main
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group_denise + scale(poly(Age,2)) + Income + MHT_used + site, data=df_G))
results[["sampleG"]] = data.frame(
        Name = "sampleG",
        Sample = "G",
        coef_POST = lm_result$coefficients["Menopause_group_denisePOST", "Estimate"],
        coef_SURG = lm_result$coefficients["Menopause_group_deniseSURG", "Estimate"],
        pval_POST = lm_result$coefficients["Menopause_group_denisePOST", "Pr(>|t|)"],
        pval_SURG = lm_result$coefficients["Menopause_group_deniseSURG", "Pr(>|t|)"]
    )

print(t(results[["sampleG"]]))

############ Merge results
results = do.call(rbind, results)
fwrite(results, "./results/group_diffs.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#endregion

#region Effect of age at menopause and MHT use on WMHV

results = list()

# Age (POST-SURG)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age) + Income + MHT_used + site, data=df_D))
results[["Age"]] = data.frame(
        Name = "Age",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "age",
        coef = lm_result$coefficients["scale(Age)", "Estimate"],
        pval = lm_result$coefficients["scale(Age)", "Pr(>|t|)"]
    )

print(t(results[["Age"]]))

# Age (POST)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age) + Income + MHT_used + site, data=df_D %>% filter(Menopause_group == "POST")))
results[["Age_POST"]] = data.frame(
        Name = "Age_POST",
        Sample = "D",
        POST_SURG = "POST",
        coef_meaning = "age",
        coef = lm_result$coefficients["scale(Age)", "Estimate"],
        pval = lm_result$coefficients["scale(Age)", "Pr(>|t|)"]
    )

print(t(results[["Age_POST"]]))

# Age (SURG)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age) + Income + MHT_used + site, data=df_D %>% filter(Menopause_group == "SURG")))
results[["Age_SURG"]] = data.frame(
        Name = "Age_SURG",
        Sample = "D",
        POST_SURG = "SURG",
        coef_meaning = "age",
        coef = lm_result$coefficients["scale(Age)", "Estimate"],
        pval = lm_result$coefficients["scale(Age)", "Pr(>|t|)"]
    )

print(t(results[["Age_SURG"]]))

# Age at menopause (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age_menopause) + scale(poly(Age, 2)) + Income + MHT_used + site, data=df_D))
results[["Age_meno"]] = data.frame(
        Name = "Age_meno",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "age_menopause",
        coef = lm_result$coefficients["scale(Age_menopause)", "Estimate"],
        pval = lm_result$coefficients["scale(Age_menopause)", "Pr(>|t|)"]
    )

print(t(results[["Age_meno"]]))

# Age at menopause in POST (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age_menopause) + scale(poly(Age, 2)) + Income + MHT_used + site, data=df_D %>% filter(Menopause_group == "POST")))
results[["Age_meno_POST"]] = data.frame(
        Name = "Age_meno_POST",
        Sample = "D",
        POST_SURG = "POST",
        coef_meaning = "age_menopause",
        coef = lm_result$coefficients["scale(Age_menopause)", "Estimate"],
        pval = lm_result$coefficients["scale(Age_menopause)", "Pr(>|t|)"]
    )

print(t(results[["Age_meno_POST"]]))

# Age at menopause in SURG (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age_menopause) + scale(poly(Age, 2)) + Income + MHT_used + site, data=df_D %>% filter(Menopause_group == "SURG")))
results[["Age_meno_SURG"]] = data.frame(
        Name = "Age_meno_SURG",
        Sample = "D",
        POST_SURG = "SURG",
        coef_meaning = "age_menopause",
        coef = lm_result$coefficients["scale(Age_menopause)", "Estimate"],
        pval = lm_result$coefficients["scale(Age_menopause)", "Pr(>|t|)"]
    )

print(t(results[["Age_meno_SURG"]]))

# Age at menopause by group interaction (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(Age_menopause) * Menopause_group + scale(poly(Age, 2)) + Income + MHT_used + site, data=df_D))
results[["Age_meno_groupInter"]] = data.frame(
        Name = "Age_meno_groupInter",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "age_menopause_by_group",
        coef = lm_result$coefficients["scale(Age_menopause):Menopause_groupSURG", "Estimate"],
        pval = lm_result$coefficients["scale(Age_menopause):Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["Age_meno_groupInter"]]))

# MHT (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used + Menopause_group + scale(poly(Age,2)) + Income + Age_menopause + site, data=df_D))
results[["MHT_covGroup"]] = data.frame(
        Name = "MHT_covGroup",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "MHT",
        coef = lm_result$coefficients["MHT_usedYes", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes", "Pr(>|t|)"]
    )

print(t(results[["MHT_covGroup"]]))

# MHT in POST (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used + scale(poly(Age,2)) + Income + Age_menopause + site, data=df_D %>% filter(Menopause_group == "POST")))
results[["MHT_inPOST"]] = data.frame(
        Name = "MHT_inPOST",
        Sample = "D",
        POST_SURG = "POST",
        coef_meaning = "MHT",
        coef = lm_result$coefficients["MHT_usedYes", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes", "Pr(>|t|)"]
    )

print(t(results[["MHT_inPOST"]]))

# MHT in SURG (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used + scale(poly(Age,2)) + Income + Age_menopause + site, data=df_D %>% filter(Menopause_group == "SURG")))
results[["MHT_inSURG"]] = data.frame(
        Name = "MHT_inSURG",
        Sample = "D",
        POST_SURG = "SURG",
        coef_meaning = "MHT",
        coef = lm_result$coefficients["MHT_usedYes", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes", "Pr(>|t|)"]
    )

print(t(results[["MHT_inSURG"]]))

# MHT group interaction (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used * Menopause_group + scale(poly(Age,2)) + Income + Age_menopause + site, data=df_D))
results[["MHT_groupInter"]] = data.frame(
        Name = "MHT_groupInter",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "MHT_by_group",
        coef = lm_result$coefficients["MHT_usedYes:Menopause_groupSURG", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes:Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["MHT_groupInter"]]))

# MHT by age at menopause interaction in POST (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used * Age_menopause + scale(poly(Age,2)) + Income + site, data=df_D  %>% filter(Menopause_group == "POST")))
results[["MHT_ageMenoInter_POST"]] = data.frame(
        Name = "MHT_ageMenoInter_POST",
        Sample = "D",
        POST_SURG = "POST",
        coef_meaning = "MHT_by_age_menopause",
        coef = lm_result$coefficients["MHT_usedYes:Age_menopause", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes:Age_menopause", "Pr(>|t|)"]
    )

print(t(results[["MHT_ageMenoInter_POST"]]))

# MHT by age at menopause interaction in SURG (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ MHT_used * Age_menopause + scale(poly(Age,2)) + Income + site, data=df_D  %>% filter(Menopause_group == "SURG")))
results[["MHT_ageMenoInter_SURG"]] = data.frame(
        Name = "MHT_ageMenoInter_SURG",
        Sample = "D",
        POST_SURG = "SURG",
        coef_meaning = "MHT_by_age_menopause",
        coef = lm_result$coefficients["MHT_usedYes:Age_menopause", "Estimate"],
        pval = lm_result$coefficients["MHT_usedYes:Age_menopause", "Pr(>|t|)"]
    )

print(t(results[["MHT_ageMenoInter_SURG"]]))

# Menopause status (sample D)
lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data=df_D))
results[["Meno_group"]] = data.frame(
        Name = "Meno_group",
        Sample = "D",
        POST_SURG = "POST_SURG",
        coef_meaning = "Meno_group",
        coef = lm_result$coefficients["Menopause_groupSURG", "Estimate"],
        pval = lm_result$coefficients["Menopause_groupSURG", "Pr(>|t|)"]
    )

print(t(results[["Meno_group"]]))

# Merge results
results = do.call(rbind, results)
fwrite(results, "./results/MHT_ageMeno.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#endregion

#region Effect of menopausal status on the relationships between cardiometabolic factors and WMHV

factor_list = c(
    "Alcohol_frequency",
    "Alcohol_status",
    "Avg_diastolic_BP",
    "Avg_systolic_BP",
    "High_BP",
    "Diabetes",
    "BMI",
    "WHR",
    "Days_moderate_activity",
    "Duration_moderate_activity",
    "Days_vigorous_activity",
    "Duration_vigorous_activity",
    "Days_walked",
    "Duration_walks",
    "Pack_years_smoking",
    "Ever_smoked",
    "Smoking_status"
)

# Simple effect of cardiometabolic factor on WMHV

results = list()

for (f in 1:length(factor_list)) {
    print(factor_list[f])

    df_tmp = df_A %>%
        rename(cardio_factor = factor_list[f]) %>%
        filter(is.na(cardio_factor) == FALSE)
    
    lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(as.numeric(cardio_factor)) + scale(poly(Age,2)) + Income + MHT_used + site, data = df_tmp))

    results[[factor_list[f]]] = data.frame(
            Factor = factor_list[f],
            Sample = "A",
            coef = lm_result$coefficients[2,"Estimate"],
            pval = lm_result$coefficients[2,"Pr(>|t|)"],
            n_PRE = nrow(df_tmp %>% filter(Menopause_group == "PRE")),
            n_POST = nrow(df_tmp %>% filter(Menopause_group == "POST")),
            n_SURG = nrow(df_tmp %>% filter(Menopause_group == "SURG"))
        )
}

results = do.call(rbind, results)
fwrite(results, "./results/cardio_on_WMHV.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Cardiometabolic factors: group interaction

results = list()

for (f in 1:length(factor_list)) {
    print(factor_list[f])

    df_tmp = df_B %>%
        rename(cardio_factor = factor_list[f]) %>%
        filter(is.na(cardio_factor) == FALSE)

    lm_result = summary(lm(scale(WMH_divTBV_log) ~ scale(as.numeric(cardio_factor)) * Menopause_group + scale(poly(Age,2)) + Income + MHT_used + site, data = df_tmp))

    results[[factor_list[f]]] = data.frame(
            Factor = factor_list[f],
            Sample = "B",
            coef_post = lm_result$coefficients[11,"Estimate"],
            coef_surg = lm_result$coefficients[12,"Estimate"],
            pval_post = lm_result$coefficients[11,"Pr(>|t|)"],
            pval_surg = lm_result$coefficients[12,"Pr(>|t|)"],
            n_PRE = nrow(df_tmp %>% filter(Menopause_group == "PRE")),
            n_POST = nrow(df_tmp %>% filter(Menopause_group == "POST")),
            n_SURG = nrow(df_tmp %>% filter(Menopause_group == "SURG"))
        )
}

results = do.call(rbind, results)
fwrite(results, "./results/cardio_WMHV_group_inter.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# BP medication by group interaction
df_tmp = df_B %>%
    filter(!is.na(BP_med) & !is.na(Avg_systolic_BP) & !is.na(Avg_diastolic_BP))

lm_result = summary(lm(scale(WMH_divTBV_log) ~ Menopause_group * scale(as.numeric(BP_med)) + Avg_systolic_BP + Avg_diastolic_BP + scale(poly(Age,2)) + Income + site, data=df_tmp))
results = data.frame(
        Name = "BPmed_byGroup",
        Sample = "B",
        coef_post = lm_result$coefficients["Menopause_groupPOST:scale(as.numeric(BP_med))","Estimate"],
        coef_surg = lm_result$coefficients["Menopause_groupSURG:scale(as.numeric(BP_med))","Estimate"],
        pval_post = lm_result$coefficients["Menopause_groupPOST:scale(as.numeric(BP_med))","Pr(>|t|)"],
        pval_surg = lm_result$coefficients["Menopause_groupSURG:scale(as.numeric(BP_med))","Pr(>|t|)"],
        n_PRE = nrow(df_tmp %>% filter(Menopause_group == "PRE")),
        n_POST = nrow(df_tmp %>% filter(Menopause_group == "POST")),
        n_SURG = nrow(df_tmp %>% filter(Menopause_group == "SURG"))
    )

fwrite(results, "./results/BPmed_group_inter.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#endregion

#region Figures

########### Figure 2

# Figure 2A: WMHV ~ Age - Unmatched sample (sample A)
ggplot(df_A, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group)) +
    geom_point(size=0.5, alpha=0.2) +
    geom_smooth(method="lm", formula=y ~ poly(x, 1)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_A_lin.png", height=5, width=5)

ggplot(df_A, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group)) +
    geom_point(size=0.5, alpha=0.2) +
    geom_smooth(method="lm", formula=y ~ poly(x, 2)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_A_poly.png", height=5, width=5)

# Figure 2B: WMHV ~ Age - Matched sample (sample B)
ggplot(df_B, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group)) +
    geom_point(size=1, alpha=0.5) +
    geom_smooth(method="lm", formula=y ~ poly(x, 2)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_B_poly.png", height=5, width=5)

# Figure 2C: WMHV ~ Age - Lohner sample (sample E)
ggplot(df_E, aes(x=Age, y=UKB_WMH_divWM_log, color=lohner_Menopause_group)) +
    geom_point(size=0.5, alpha=0.2) +
    geom_smooth(method="lm", formula=y ~ poly(x, 2)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_E_poly.png", height=5, width=5)

# Figure 2D: WMHV ~ Age - >5 years after menopause (sample F)
ggplot(df_F, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group)) +
    geom_point(size=1, alpha=0.5) +
    geom_smooth(method="lm", formula=y ~ poly(x, 2)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_F_poly.png", height=5, width=5)

# Figure 2E: WMHV - Matched sample (sample B)
ggplot(df_B, aes(x=Menopause_group, y=WMH_divTBV_log, fill=Menopause_group)) +
    geom_violin(trim=TRUE, scale="width") + 
    geom_boxplot(width=0.2, color="black", fill="white") + 
    scale_fill_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_discrete(name="Menopausal group") +
    guides(fill = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/WMHV_group.png", height=5, width=5)

# Figure 2F: MHT and age at menopause
results = as.data.frame(fread("./results/MHT_ageMeno.tsv"))
empty_row_1 = data.frame(Name = NA, Sample = NA, POST_SURG = "POST", coef_meaning = "Meno_group", coef = NA, pval = NA)
empty_row_2 = data.frame(Name = NA, Sample = NA, POST_SURG = "SURG", coef_meaning = "Meno_group", coef = NA, pval = NA)
empty_row_3 = data.frame(Name = NA, Sample = NA, POST_SURG = "POST_SURG", coef_meaning = "MHT_by_age_menopause", coef = NA, pval = NA)
empty_row_4 = data.frame(Name = NA, Sample = NA, POST_SURG = "POST", coef_meaning = "age_menopause_by_group", coef = NA, pval = NA)
empty_row_5 = data.frame(Name = NA, Sample = NA, POST_SURG = "SURG", coef_meaning = "age_menopause_by_group", coef = NA, pval = NA)

results = results %>%
    # Add empty rows
    bind_rows(list(empty_row_1, empty_row_2, empty_row_3, empty_row_4, empty_row_5)) %>%
    filter(coef_meaning != "MHT_by_group") %>%
    mutate(
        POST_SURG = factor(POST_SURG, levels=c("POST_SURG", "POST", "SURG"), labels=c("Matched POST-SURG", "POST", "SURG")),
        coef_meaning = factor(coef_meaning, levels=c(
                "age", "age_menopause", "age_menopause_by_group", "MHT", "MHT_by_age_menopause", "Meno_group"
            ), labels=c(
                "Age", "Age at menopause", "Group * Age at menopause", "MHT", "MHT * Age at menopause", "Menopausal status"
            ))
    ) %>%
    mutate(
        pval_fdr = p.adjust(pval, method = "fdr"),
        pval_fdr_sig = ifelse(pval_fdr < 0.05, pval_fdr, NA),
        coef_sig = ifelse(pval_fdr < 0.05, coef, NA)
    ) %>%
    mutate(
        Label = paste0("b = ", round(coef,2))
    ) %>%
    glimpse()

ggplot(results, aes(x=POST_SURG, y=coef_meaning, fill=coef_sig, label=Label)) +
    geom_tile(color = "white", lwd = 2, linetype = 1) +
    geom_text(size = 5) +
    scale_fill_gradientn(colours=cet_pal(256, name = "d1", alpha = 1), limits=c(-0.6,0.6), breaks=c(-0.6,0,0.6), name="std(B)") +
    scale_x_discrete(labels=c("WMHV", "WMHV", "WMHV")) +
    facet_grid(cols = vars(POST_SURG), scales = "free", space = "free") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=0.1))
ggsave("./visualization/MHT_age_meno.png", height = 5, width = 9)


# TEST: sample G (matched IDs of Denise's unmatched sample from paper)
ggplot(df_G, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group_denise)) +
    geom_point(size=0.5, alpha=0.2) +
    geom_smooth(method="lm", formula=y ~ poly(x, 1)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_G_lin.png", height=5, width=5)

ggplot(df_G, aes(x=Age, y=WMH_divTBV_log, color=Menopause_group_denise)) +
    geom_point(size=0.5, alpha=0.2) +
    geom_smooth(method="lm", formula=y ~ poly(x, 2)) +
    scale_color_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_continuous(name="Age") +
    guides(color = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/age_G_poly.png", height=5, width=5)

########### Figure 3

# Figure 3A: WMHV ~ cardiometabolic factors
results = as.data.frame(fread("./results/cardio_on_WMHV.tsv"))
results = results %>%
    mutate(
        Factor = factor(Factor, 
            levels=c(
                "Alcohol_frequency", "Alcohol_status",
                "Avg_diastolic_BP", "Avg_systolic_BP", "High_BP",
                "Diabetes",
                "BMI", "WHR",
                "Days_walked", "Duration_walks", "Days_moderate_activity", "Duration_moderate_activity", "Days_vigorous_activity", "Duration_vigorous_activity",
                "Smoking_status", "Ever_smoked", "Pack_years_smoking"
            ),
            labels = c(
                "Alcohol frequency", "Alcohol status",
                "Diastolic BP", "Systolic BP", "Hypertension",
                "Diabetes",
                "BMI", "WHR",
                "Days walked", "Duration walks", "Days moderate activity", "Duration moderate activity", "Days vigorous activity", "Duration vigorous activity",
                "Smoking status", "Ever smoked", "Pack years smoking"
            ))
    ) %>%
    mutate(
        categ = factor(case_when(
            Factor %in% c("Alcohol frequency", "Alcohol status") ~ "Alcohol",
            Factor %in% c("Diastolic BP", "Systolic BP", "Hypertension") ~ "Blood pressure (BP)",
            Factor %in% c("Diabetes") ~ "Diabetes",
            Factor %in% c("BMI", "WHR") ~ "Body size",
            Factor %in% c("Days walked", "Duration walks", "Days moderate activity", "Duration moderate activity", "Days vigorous activity", "Duration vigorous activity") ~ "Physical activity",
            Factor %in% c("Smoking status", "Ever smoked", "Pack years smoking") ~ "Smoking",
        ))
    ) %>%
    mutate(
        dummy_var = factor("WMHV")
    ) %>%
    # FDR correction
    mutate(
        pval_fdr = p.adjust(pval, method = "fdr"),
        pval_fdr_sig = ifelse(pval_fdr < 0.05, pval_fdr, NA),
        coef_sig = ifelse(pval_fdr < 0.05, coef, NA)
    ) %>%
    mutate(
        Label = paste0("b = ", round(coef,2))
    ) %>%
    glimpse()

ggplot(results, aes(x=Factor, y=dummy_var, fill=coef_sig, label=Label)) +
    geom_tile(color = "white", lwd = 0.1, linetype = 1) +
    geom_text(angle = 90, size = 5) +
    # scale_fill_gradientn(colors = viridis_pal(direction = -1, option = "A", begin=0.3, end=1)(9), name = "std(B)", limits = c(0, 0.15)) +
    scale_fill_gradientn(colours=cet_pal(256, name = "d1", alpha = 1), limits=c(-0.15,0.15), breaks=c(-0.15,0,0.15), name="std(B)") +
    facet_grid(cols = vars(categ), scales = "free", space = "free") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 52, vjust = 1, hjust=1, size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size=12))
ggsave("./visualization/cardio_WMHV.png", height = 4, width = 14)

# Figure 3B: WNHV ~ cardio * group
results = as.data.frame(fread("./results/cardio_WMHV_group_inter.tsv"))
results = results %>%
    mutate(
        Factor = factor(Factor, 
            levels=c(
                "Alcohol_frequency", "Alcohol_status",
                "Avg_diastolic_BP", "Avg_systolic_BP", "High_BP",
                "Diabetes",
                "BMI", "WHR",
                "Days_walked", "Duration_walks", "Days_moderate_activity", "Duration_moderate_activity", "Days_vigorous_activity", "Duration_vigorous_activity",
                "Smoking_status", "Ever_smoked", "Pack_years_smoking"
            ),
            labels = c(
                "Alcohol frequency", "Alcohol status",
                "Diastolic BP", "Systolic BP", "Hypertension",
                "Diabetes",
                "BMI", "WHR",
                "Days walked", "Duration walks", "Days moderate activity", "Duration moderate activity", "Days vigorous activity", "Duration vigorous activity",
                "Smoking status", "Ever smoked", "Pack years smoking"
            ))
    ) %>%
    mutate(
        categ = factor(case_when(
            Factor %in% c("Alcohol frequency", "Alcohol status") ~ "Alcohol",
            Factor %in% c("Diastolic BP", "Systolic BP", "Hypertension") ~ "Blood pressure (BP)",
            Factor %in% c("Diabetes") ~ "Diabetes",
            Factor %in% c("BMI", "WHR") ~ "Body size",
            Factor %in% c("Days walked", "Duration walks", "Days moderate activity", "Duration moderate activity", "Days vigorous activity", "Duration vigorous activity") ~ "Physical activity",
            Factor %in% c("Smoking status", "Ever smoked", "Pack years smoking") ~ "Smoking",
        ))
    ) %>%
    pivot_longer(cols=c(3,4,5,6), names_sep = "_", names_to = c(".value", "group")) %>%
    mutate(group = factor(group, levels=c("post", "surg"), labels=c("POST", "SURG"))) %>%
    # FDR correction
    mutate(
        pval_fdr = p.adjust(pval, method = "fdr"),
        pval_fdr_sig = ifelse(pval_fdr < 0.05, pval_fdr, NA),
        coef_sig = ifelse(pval_fdr < 0.05, coef, NA)
    ) %>%
    mutate(
        Label = paste0("b = ", round(coef,2))
    ) %>%
    glimpse()

# Handle case when all coefficients non-significant (NAs)
if(is.logical(results$coef_sig)) {
    plt = ggplot(results, aes(x=Factor, y=group, label=Label)) +
        geom_tile(fill = "#7f7f7f", color = "white", lwd = 0.1, linetype = 1)
} else {
    plt = ggplot(results, aes(x=Factor, y=group, fill=coef_sig, label=Label)) +
        geom_tile(color = "white", lwd = 0.1, linetype = 1)
}

plt +
    geom_text(angle = 90, size = 5) +
    scale_fill_gradientn(colours=cet_pal(256, name = "d1", alpha = 1), limits=c(-0.3,0.3), breaks=c(-0.3,0,0.3), name="std(B)") +
    facet_grid(rows = vars(group), cols = vars(categ), scales = "free", space = "free") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 52, vjust = 1, hjust=1, size = 12),
        legend.text = element_text(size=12), legend.title = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=0.1))
ggsave("./visualization/cardio_WMHV_byGroup.png", height = 6, width = 14)

# Figure 3C: Interaction between BP medication and menopause group on WMHV
df_tmp = df_B %>%
    filter(!is.na(BP_med) & !is.na(Avg_systolic_BP) & !is.na(Avg_diastolic_BP))

ggplot(df_tmp, aes(x=BP_med, y=WMH_divTBV_log, fill=Menopause_group)) +
    geom_violin(trim=TRUE, scale="width") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="No", Menopause_group=="PRE"), width=0.05, color="black", position=position_nudge(x=-0.3), fill="white") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="No", Menopause_group=="POST"), width=0.05, color="black", position=position_nudge(x=0), fill="white") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="No", Menopause_group=="SURG"), width=0.05, color="black", position=position_nudge(x=0.3), fill="white") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="Yes", Menopause_group=="PRE"), width=0.05, color="black", position=position_nudge(x=-0.3), fill="white") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="Yes", Menopause_group=="POST"), width=0.05, color="black", position=position_nudge(x=0), fill="white") + 
    geom_boxplot(data = df_tmp %>% filter(BP_med=="Yes", Menopause_group=="SURG"), width=0.05, color="black", position=position_nudge(x=0.3), fill="white") + 
    scale_fill_manual(values = color_scale_groups) +
    scale_y_continuous(name="White Matter Hyperintensity Volume") +
    scale_x_discrete(name="Blood Pressure Medication") +
    guides(fill = "none") + 
    theme_classic() +
    theme(text=element_text(size=15))
ggsave("./visualization/BPmed_by_group.png", width=5, height=5)

#endregion
