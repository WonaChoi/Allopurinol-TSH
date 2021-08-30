

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------- Analysis of the relationship between allopurinol intake and elevated thyroid-stimulating hormone level ---------------------- # 
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# rm(list = ls(all = TRUE))

## Installation of required packages ####
wants <- c("RPostgreSQL", "dplyr", "readr", "SqlRender", "lubridate", "RDocumentation", "survival", "reshape2", "stringr", "purrr", "broom", "tidyr")
has <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

# Packages loading
invisible(lapply(wants, library, character.only = TRUE))

############### ★ Please modify the below section manually ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ ###############

HOSPITAL = "CMCS"

#### ... 1. Database connection setting ####
mydbtype = "postgresql"     # db type: "postgresql", "oracle", "pdw", "redshift", "impala", "netezza", "bigquery", "redshift", "sqlite", "mssql", "sql server"
mydbname = "dbname"         # database name
myschemaname = "schema"     # schema name
myhost = ""      # db host
myport = 5432               # db port
myuserid = ""       # ID    
mypassword = ""  # password
myvocabschemaname='vocab'   # name of the vocabulary schema

#### ... 2. EDI code adjustment ####
procedure_edi_code ='procedure_source_value' # column name of lab EDI code in 'procedure occurrence' table


#### ... 3. Measurement code adjustment ####  
source_column = 'measurement_source_value'
# Thyroid hormones
T3code <- c(3010340, 3030897, 3008304) # LOINC
T4code <- c(3016451, 3016991, 3032367, 3024908, 3014620, 3003045) # LOINC
freeT4code <- c(1175900, 3008486, 3008598, 3032600, 3016363, 40761932) # LOINC

#### 0. Preparation #########################################################################################################################################
#### ... 1. DB connection ####

drv = dbDriver("PostgreSQL") # "MySQL", "PostgreSQL", "RSQLite"
con = dbConnect(drv, dbname= mydbname, host= myhost, port= myport, user= myuserid, password= mypassword)


#### ... 2. Loading reference files ####
allocode     <- read.csv("data/allopurinol_list_2.csv")
ccicode      <- read.csv("data/CCI_list.csv")
covdiscode   <- read.csv("data/cov_disease_snomed_list_2.csv")
covdrcode    <- read.csv("data/cov_drug_list_2.csv") # RxNorm
covmcode     <- read.csv("data/cov_meas_list.csv")
expcode_con  <- read.csv("data/exception_disease_snomed_list.csv")
expcode_drug <- read.csv("data/exception_drug_list.csv")
TSH_code     <- read.csv("data/Thyrotropin_final_list_2.csv")
TSHlist      <- TSH_code[,1] %>% as.vector
thyroid_operation <- read.csv("data/thyroidectomy.csv")


##### 1. Data extracting ###################################################################################################################################

sql <- "SELECT COUNT(DISTINCT person_id)
        FROM @A.person"
diagram_1_entire_patient <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                   A= myschemaname)) %>% .[[1]]


sql <- "SELECT drug_exposure_id, person_id,drug_concept_id, drug_exposure_start_date, drug_exposure_end_date, days_supply
        FROM @A.drug_exposure
        WHERE drug_concept_id IN (@B);"
allo_drug_exposure_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                  A= myschemaname,
                                                  B= allocode[,1])) 


## ... 1. Date of the first TSH examination : cohort-in date ####
sql <- "SELECT *
        FROM @A.measurement
        WHERE measurement_concept_id IN (@B);"
TSH <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
                              A= myschemaname, 
                              B= TSHlist))

TSH <- TSH %>% as_tibble() %>% mutate(value_source_value= word(value_source_value, 1) %>% as.numeric) 
TSH <- TSH %>% filter(!is.na(value_source_value))

TSH_uniquelist <- TSH %>% group_by(person_id) %>% select(person_id, measurement_date) %>% rename('enroll_date' = 'measurement_date') %>% arrange(enroll_date) %>% slice(1)

diagram_2_TSH <- TSH_uniquelist %>% nrow

enroll0 <- TSH

## ... 2. Person table: List of persons with at least one TSH examination record  #####
sql <- "SELECT person_id, year_of_birth, month_of_birth, day_of_birth, gender_concept_id
        FROM @A.person
        WHERE person_id IN (@B)"
person <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                 A= myschemaname,
                                 B= TSH_uniquelist %>% pull(person_id))) %>% as_tibble(); remove(sql)


TSH_all <- left_join(TSH_uniquelist, person, by='person_id') %>% as_tibble() %>% 
  mutate(birth = make_date(year_of_birth,month_of_birth,day_of_birth),
         age = time_length(enroll_date - birth, "year")) %>% filter(age >= 19)


diagram_3_adult <- TSH_all %>% pull(person_id) %>% length
#### 2. Exclusion criteria  #####

#### ... 1. ID of person whose TSH test results are less than 0.5 mμIU/ml or more than 5 mIU/ml for one year after enroll date ####
exclusion1 <- 
  TSH %>% select(person_id, measurement_concept_id, measurement_date, value_source_value) %>% 
  left_join(TSH_all %>% select(person_id,enroll_date,age), by='person_id') %>% as_tibble %>% 
  filter(measurement_date-enroll_date <= 365-1 & (value_source_value > 4.5 | value_source_value < .5)) %>% 
  select(person_id)

enroll2 <- TSH_all %>% anti_join(exclusion1,by="person_id")

diagram_4_exclusion1 <- enroll2 %>% pull(person_id) %>% length

##### ... 2. Excluding people diagnosed with a condition that could affect thyroid function at least once before 'cohort-in date + 365 days' ####
#  graves disease, ptuitary tumor, hypopituitarism, sheehan's syndrome 
sql <- "SELECT person_id, condition_concept_id, condition_start_date, 1 AS expcon
        FROM @A.condition_occurrence
        WHERE condition_concept_id IN (@B);"
ex2_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
                                  A= myschemaname, 
                                  B= expcode_con[,1]))

if (count(enroll2 %>% left_join(ex2_tmp, by='person_id') %>% 
          filter(!is.na(condition_start_date)) %>% 
          filter(condition_start_date<=(enroll_date+365-1)))==0) {
  enroll3 <- enroll2
} else {
  enroll3 <- 
    anti_join(enroll2, 
              enroll2 %>% left_join(ex2_tmp, by='person_id') %>% 
                filter(!is.na(condition_start_date)) %>% filter(condition_start_date<=(enroll_date+365-1)),
              by='person_id')
}

diagram_4_exclusion2 <- enroll3 %>% pull(person_id) %>% length


#### ... 3. Anyone who had a thyroidectomy prior to 'cohort-in date + 365 days' #####
# A person who had a thyroid surgery for a year after enroll date.

# EDI
sql <- "SELECT person_id, @B, procedure_date, 1 AS exppro
        FROM @A.procedure_occurrence
        WHERE @B like '%P4551%' OR @B like '%P4552%' OR @B like '%P4553%' OR @B like '%P4554%';"
ex3_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
                                  A= myschemaname,
                                  B= procedure_edi_code))

# SNOMED
# sql <- "SELECT person_id, procedure_concept_id, procedure_date, 1 AS exppro
#         FROM @A.procedure_occurrence
#         WHERE procedure_concept_id IN (@B);"
# 
# ex3_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
#                                   A= myschemaname,
#                                   B= thyroid_operation[,2]))


if(count(enroll3 %>% inner_join(ex3_tmp,by='person_id') %>% filter(procedure_date <= (enroll_date+365-1)))==0) {
  enroll4 <- enroll3
} else {
  enroll4 <-
    enroll3 %>% anti_join(enroll3 %>% inner_join(ex3_tmp,by='person_id') %>% filter(procedure_date <= (enroll_date+365-1)), 
                          by='person_id')
}

diagram_4_exclusion3 <- enroll4 %>% pull(person_id) %>% length

#### ... 4. Exclude persons who have taken the excluded drugs between 'cohort-in date ~ cohort-in date+365'  #####
# antithyroid drug, Levothyroxine 
sql <- "select person_id, drug_concept_id, drug_exposure_start_date, drug_exposure_end_date, days_supply
        FROM @A.drug_exposure
        WHERE drug_concept_id in (@B);"
ex4_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                  A= myschemaname, 
                                  B= expcode_drug %>% pull(concept_id)))


if(count(ex4_tmp)==0) {
  ex4_tmp2 <- ex4_tmp 
} else {
  ex4_tmp2  <- 
    ex4_tmp %>% left_join(TSH_uniquelist, by = "person_id") %>% 
    mutate(start_diffday = time_length(interval(enroll_date,drug_exposure_start_date ),"day"),
           end_diffday = time_length(interval(enroll_date,drug_exposure_end_date),"day")) %>% 
    filter(between(start_diffday,0,365-1) | between(end_diffday,0,365-1))
}


if (count(ex4_tmp2)==0) {
  enroll4_2 <- enroll4
} else {
  enroll4_2 <- anti_join(enroll4,ex4_tmp2, by='person_id')
}


diagram_4_exclusion4 <- enroll4_2 %>% pull(person_id) %>% length

#####  ... 5.  A person with a history of radiation iodine therapy. ######

iodine131 <- c(4036252, 4333348, 4331674, 4161520, 4040441, 4161255)

sql <- "SELECT person_id, procedure_concept_id, procedure_date
        FROM @A.procedure_occurrence
        WHERE procedure_concept_id in (@B)"
ex4_2_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
                                    A= myschemaname,
                                    B= iodine131)) 


if (count(ex4_2_tmp) == 0) {
  enroll5 <- enroll4_2
} else {
  enroll5 <- enroll4_2 %>% 
    anti_join(enroll4_2 %>% inner_join(ex4_2_tmp, by='person_id') %>% 
                filter(procedure_date <= (enroll_date+365-1)), by='person_id')
}

diagram_4_exclusion5 <- enroll5 %>% pull(person_id) %>% length


#### 3. Extracting case group ###########################################################################################################################

allo_case <- 
  enroll5 %>% inner_join(enroll0 %>% select(person_id, measurement_date, value_source_value) %>% filter(value_source_value > 4.5), 
                         by='person_id') %>% 
  mutate(age = time_length(enroll_date - birth, "year") %>% trunc) %>% 
  arrange(person_id, measurement_date) %>% group_by(person_id) %>% slice(1) %>% 
  select(person_id, enroll_date, measurement_date, birth, age, value_source_value,gender_concept_id) %>% 
  group_by(person_id) %>% mutate(age= time_length(interval(min(birth), min(measurement_date)),'year') %>% trunc, 
                                 TSH = 1)  %>% 
  dplyr::rename('indexdt'='measurement_date')

allo_case$case_number <- 1:nrow(allo_case)

#### 4. Extracting control group ##########################################################################################################################

#### ... 1. The potential subjects for the control group ####

enroll5_m <- 
  enroll5 %>% 
  left_join(enroll0 %>% select(person_id, measurement_date, value_source_value), by="person_id") %>% 
  left_join(allo_case %>% select(person_id, indexdt), by="person_id") %>% 
  select(person_id, enroll_date, birth, age, value_source_value, measurement_date, gender_concept_id, indexdt) %>% 
  mutate(age = age %>% trunc) %>% 
  distinct()

#### ... 2. Making empty dataframe  ####
allo_case_tmp <- allo_case %>% 
  select(person_id, enroll_date, indexdt, birth,age,value_source_value,gender_concept_id,TSH, case_number) 

allo_control <- allo_case[1,] %>% 
  select(person_id,enroll_date,birth,age,value_source_value,gender_concept_id,case_number,indexdt) %>% 
  .[0,]


enroll5_m2 <- enroll5_m %>% group_by(person_id) %>% 
  filter(measurement_date == min(measurement_date)) %>% ungroup()

allo_case_patient_number <- length(allo_case$person_id)

#### ... 3. Control group selection ####
for (i in 1:allo_case_patient_number) {
  # patient
  target_patient <- allo_case_tmp[i,]
  # control
  allo_control_tmp1 <- enroll5_m2 %>%
    filter((gender_concept_id == target_patient[7] %>% as.numeric) &  # gender
             (age >= target_patient[5] %>% as.numeric - 5 & age <= target_patient[5] %>% as.numeric + 5) & # age +- 5
             (enroll_date %within% interval(target_patient[2][[1]] - days(30), target_patient[2][[1]] + days(30))) & # enroll date near 1 month 
             (is.na(indexdt) | (!is.na(indexdt) & target_patient[3][[1]] < indexdt)) &
             (!(person_id %in% (allo_control %>% pull(person_id)))) &
             (!(person_id %in% (allo_case_tmp %>% pull(person_id)))) &
             (value_source_value >= 0.5 & value_source_value <= 4.5)) 
  
  sample_n_number <- if_else(nrow(allo_control_tmp1) >= 4, 4, nrow(allo_control_tmp1) %>% as.double())
  allo_control_tmp2 <- allo_control_tmp1 %>% sample_n(sample_n_number) %>% mutate(case_number = i)
  
  allo_control <- bind_rows(allo_control, allo_control_tmp2)
  # print
  print(paste0("|",strrep('■',round(i/allo_case_patient_number,2)*100) %>% noquote(), strrep('-',round(1- i/allo_case_patient_number,2)*100) %>% noquote(),"| ", format(round(i/allo_case_patient_number,4)*100,nsmall=3),"% [",i,"/",allo_case_patient_number,"]"))
  print(target_patient %>% as.data.frame)
  print(allo_control_tmp2 %>% as.data.frame())
}


allo_control %>% group_by(case_number) %>% summarise(count=n()) %>% group_by(count) %>% summarise(n=n()) # 대조군 쌍의 수 확인

#### ... 4. Selected case & control group ####
allo_control2 <- allo_control %>% 
  left_join(allo_case %>% ungroup %>% select(indexdt,case_number) %>% rename('indexdt2'='indexdt'), 
            by='case_number') %>% 
  distinct()

allo_case2 <- allo_case %>% 
  mutate(indexdt2 = indexdt) %>% 
  distinct()

allo_casecontrol <- 
  bind_rows(allo_case2 %>% select(person_id, enroll_date, age, gender_concept_id, value_source_value, TSH, case_number, indexdt2), 
            allo_control2 %>% mutate(TSH = 0, value_source_value = NA) %>% mutate_at(vars('value_source_value'), as.numeric) %>% 
              select(person_id, enroll_date, age, gender_concept_id, value_source_value, TSH, case_number, indexdt2)) %>% 
  arrange(case_number,TSH, person_id)



#### 5. Target drug definition ############################################################################################################################

sql <- "SELECT drug_exposure_id, person_id,drug_concept_id, drug_exposure_start_date, drug_exposure_end_date, days_supply
        FROM @A.drug_exposure
        WHERE drug_concept_id IN (@B);"
allo_drug_exposure_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                  A= myschemaname,
                                                  B= allocode[,1])) 


# allopurinol patient total
allo_allallo <- allo_drug_exposure_list %>% as_tibble() %>% mutate(drug_end_date = drug_exposure_start_date + days_supply-1) %>% mutate(allo = 1) %>% 
  arrange(person_id, drug_exposure_start_date, drug_end_date) 

allo_ccallo <- allo_casecontrol %>% select(person_id, indexdt2, case_number) %>% 
  inner_join(allo_allallo %>% select(person_id, drug_exposure_start_date, drug_end_date, allo, days_supply), by='person_id') %>%  
  filter((drug_exposure_start_date %within% interval(indexdt2-180,indexdt2-1)) | (drug_end_date %within% interval(indexdt2-180,indexdt2-1))) %>% 
  distinct() %>% 
  mutate(drug_exposure_start_date2 = case_when(drug_exposure_start_date < indexdt2-180 ~ indexdt2-180, 
                                               TRUE ~ drug_exposure_start_date),
         drug_end_date2 = case_when(drug_end_date>=indexdt2 ~ indexdt2-1, 
                                    TRUE ~ drug_end_date)) %>% 
  arrange(person_id, drug_exposure_start_date2, drug_end_date2)


if (nrow(allo_ccallo)==0) { 
  allo_allow <- tibble(person_id=double(), case_number=integer(), allo=double(), w=double()) 
  allo_alldose <- tibble(person_id=double(),case_number=integer(), drugday=double())  
} else {
  allo_allow <- allo_ccallo %>% distinct() %>% 
    mutate(w = case_when((drug_exposure_start_date2 %within% interval(indexdt2-30, indexdt2-1) | drug_end_date2 %within% interval(indexdt2-30, indexdt2-1) ~ 1),
                         (drug_exposure_start_date2 %within% interval(indexdt2-60, indexdt2-30-1) | drug_end_date2 %within% interval(indexdt2-60, indexdt2-30-1) ~ 2),
                         (drug_exposure_start_date2 %within% interval(indexdt2-90, indexdt2-60-1) | drug_end_date2 %within% interval(indexdt2-90, indexdt2-60-1) ~ 3),
                         (drug_exposure_start_date2 %within% interval(indexdt2-180, indexdt2-90-1) | drug_end_date2 %within% interval(indexdt2-180, indexdt2-90-1) ~ 4)) %>% 
             as.numeric()) %>% group_by(person_id, case_number) %>% summarise(allo = max(allo), w=min(w))
  
  allo_allodose <- 
    allo_ccallo %>% 
    select(person_id, case_number, drug_end_date, drug_exposure_start_date2, days_supply, drug_end_date2) %>% 
    group_by(person_id,case_number) %>% 
    mutate(lag = lag(drug_exposure_start_date2),
           gap = drug_end_date-lag(drug_exposure_start_date2),
           cwc= drug_end_date2-lag(drug_exposure_start_date2)) %>% 
    mutate(n= if_else(cwc<=30|is.na(cwc),0,1)) %>% 
    mutate(cs = cumsum(n)) %>% 
    filter(cs==0) %>% 
    select(person_id, case_number, drug_end_date2, drug_exposure_start_date2) %>% 
    summarise(drugday = max(drug_end_date2)-min(drug_exposure_start_date2)+1)
}


### condition covariate - cci

sql <- "SELECT person_id, condition_concept_id, condition_start_date
        FROM @A.condition_occurrence 
        WHERE condition_concept_id IN (@B);"
cci_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                  A= myschemaname,
                                  B= ccicode[,1]))


allo_cci <- 
  allo_casecontrol %>% select(person_id, case_number, indexdt2) %>% 
  left_join(
    cci_tmp %>% left_join(ccicode, by=c("condition_concept_id"="concept_id")) %>% select(person_id, condition_start_date, CCI1:CCI17), 
    by='person_id') %>% 
  filter(condition_start_date %within% interval(indexdt2-180, indexdt2-1)) %>% 
  mutate_at(vars(paste0('CCI',1:17)), ~replace(.,is.na(.), 0)) %>% 
  group_by(person_id, case_number) %>% 
  summarise(CCI1=max(CCI1), CCI2=max(CCI2), CCI3=max(CCI3), CCI4=max(CCI4), 
            CCI5=max(CCI5), CCI6=max(CCI6), CCI7=max(CCI7), CCI8=max(CCI8), 
            CCI9=max(CCI9), CCI10=max(CCI10), CCI11=max(CCI11), CCI12=max(CCI12), 
            CCI13=max(CCI13), CCI14=max(CCI14), CCI15=max(CCI15), CCI16=max(CCI16), CCI17=max(CCI17)) %>% 
  distinct(person_id, case_number, .keep_all = TRUE)



### condition covariate - allo

sql <- "SELECT person_id, condition_concept_id, condition_start_date
        FROM @A.condition_occurrence 
        WHERE condition_concept_id IN (@B);"
covdis_tmp <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                     A= myschemaname,
                                     B= covdiscode[,1]))


allo_covdis <- 
  allo_casecontrol %>% select(person_id, case_number, indexdt2) %>% inner_join(covdis_tmp %>% inner_join(covdiscode, by =c("condition_concept_id"="concept_id")), by='person_id') %>% 
  filter(condition_start_date %within% interval(indexdt2-180, indexdt2-1)) %>% group_by(person_id, case_number) %>% 
  mutate_at(vars(dis1), ~replace(.,is.na(.), 0)) %>% mutate(cdis1= dis1) %>% 
  summarise(cdis1=max(cdis1))


### drug covariate
sql <- "SELECT drug_exposure_id, person_id, drug_concept_id, drug_exposure_start_date, drug_exposure_end_date, days_supply
        FROM @A.drug_exposure
        WHERE drug_concept_id IN (@B);"
cov_drug_exposure_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                 A= myschemaname,
                                                 B= covdrcode[,1])) 


allo_covdr <- 
  allo_casecontrol %>% inner_join(cov_drug_exposure_list %>% inner_join(covdrcode, by= c("drug_concept_id"="concept_id")) %>% mutate(drug_end_date = drug_exposure_start_date+days_supply-1), by='person_id') %>% 
  filter((drug_exposure_start_date %within% interval(indexdt2-180, indexdt2-1)) | (drug_end_date %within% interval(indexdt2-180,indexdt2-1))) %>% 
  group_by(person_id, case_number, drname) %>% summarise(n=n()) %>% mutate(n= if_else(n>=1,1,0)) %>% filter(!n == 0) %>%
  mutate(cdr1= if_else(drname %in% "TKI",1,0),
         cdr2= if_else(drname %in% "cancerimmunotherapy",1,0),
         cdr3= if_else(drname %in% "antituberculosis",1,0),
         cdr4= if_else(drname %in% "dobutamine",1,0),
         cdr5= if_else(drname %in% "octreotide",1,0),
         cdr6= if_else(drname %in% "interferonalfa",1,0),
         cdr7= if_else(drname %in% "lithium",1,0),
         cdr8= if_else(drname %in% "amiodarone",1,0),
         cdr9= if_else(drname %in% "azathioprine",1,0),
         cdr10= if_else(drname %in% "mercaptopurine",1,0),
         cdr11= if_else(drname %in% "warfarin",1,0),
         cdr12= if_else(drname %in% "dopamine",1,0),
         cdr13= if_else(drname %in% "metformin",1,0),
         cdr14= if_else(drname %in% "nsaid",1,0),
         cdr15= if_else(drname %in% "acetaminophen",1,0),
         cdr16= if_else(drname %in% "oxycodone",1,0),
         cdr17= if_else(drname %in% "colchicine",1,0),
         cdr18= if_else(drname %in% "corticosteroid",1,0)) %>%
  select(-n) %>% 
  group_by(person_id, case_number) %>% 
  summarise(cdr1=as.factor(max(cdr1)),
            cdr2=as.factor(max(cdr2)),
            cdr3=as.factor(max(cdr3)),
            cdr4=as.factor(max(cdr4)),
            cdr5=as.factor(max(cdr5)),
            cdr6=as.factor(max(cdr6)),
            cdr7=as.factor(max(cdr7)),
            cdr8=as.factor(max(cdr8)),
            cdr9=as.factor(max(cdr9)),
            cdr10=as.factor(max(cdr10)),
            cdr11=as.factor(max(cdr11)),
            cdr12=as.factor(max(cdr12)),
            cdr13=as.factor(max(cdr13)),
            cdr14=as.factor(max(cdr14)),
            cdr15=as.factor(max(cdr15)),
            cdr16=as.factor(max(cdr16)),
            cdr17=as.factor(max(cdr17)),
            cdr18=as.factor(max(cdr18)))



allo_t01 <- 
  allo_casecontrol %>% 
  left_join(allo_allow, by= c('person_id'='person_id','case_number'='case_number')) %>% 
  left_join(allo_allodose, by= c('person_id'='person_id','case_number'='case_number')) %>% 
  left_join(allo_cci, by= c('person_id'='person_id','case_number'='case_number')) %>% 
  left_join(allo_covdis, by= c('person_id'='person_id','case_number'='case_number')) %>% 
  left_join(allo_covdr, by= c('person_id'='person_id','case_number'='case_number')) 


allo_t1 <- 
  allo_t01 %>% mutate_at(vars(starts_with('CCI'), starts_with('cdis'), starts_with('cdr'), w, allo), ~replace(.,is.na(.),0)) %>% 
  mutate(agegr = case_when(age %>% between(19,29) ~ 1,
                           age %>% between(30,39) ~ 2,
                           age %>% between(40,49) ~ 3,
                           age %>% between(50,59) ~ 4,
                           age %>% between(60,69) ~ 5,
                           age %>% between(70,79) ~ 6,
                           TRUE ~ 7),
         enrollyear = year(enroll_date)) # 오래 걸립니다.

# create variables

t1 <- 
  allo_t1 %>% ungroup %>% 
  mutate(CCI = 1*rowSums(select(.,c(CCI6,CCI7,CCI11,CCI13))) + 2*rowSums(select(.,c(CCI2,CCI5,CCI9,CCI12,CCI14))) + 4*rowSums(select(.,c(CCI15,CCI7))) + 6*rowSums(select(.,CCI16)),
         CCIGR = case_when(CCI < 1 ~ 1,
                           CCI < 2 ~ 2,
                           CCI < 3 ~ 3,
                           TRUE ~ 4),
         gender_concept_id = if_else(gender_concept_id==8507,1,0)) %>% 
  mutate_at(vars(starts_with('CCI'), starts_with('cdis'),starts_with('cdr'), gender_concept_id, agegr, enrollyear, w, CCIGR, allo), as.factor) # 꽤 오래 걸림


# transformation to factor  

# write_rds(t1, "rds/t1.rds") 
# t1 <- read_rds("rds/t1.rds") 

#### 6. OUTPUT ####################################################################################################################################################
##### ... 1.  Frequency table  ##############################

myfreq <- 
  t1 %>% select(gender_concept_id, TSH, allo, CCI1:CCIGR) %>%
  map(~table(.x,t1$TSH) %>% as.data.frame %>% pivot_wider(names_from='Var2',values_from='Freq')) %>% 
  bind_rows(., .id="variable") %>% 
  rename("control"="0", "patient"="1", "value"=".x") %>% 
  mutate(total = control + patient,
         control_ratio = control/(control+patient),
         patient_ratio = patient/(control+patient)) %>%
  # ratio calculation
  relocate(variable, value, total, control, control_ratio, patient, patient_ratio) %>% filter(!variable=='TSH')


write_rds(myfreq, paste0("output/",HOSPITAL,"-TABLE1.rds"))
write_csv(myfreq, paste0("output/",HOSPITAL,"-TABLE1.csv"))




##### ... 2. Logistic analysis  #####################################################

#### ....... 1. Variable list ######################
# Logistic regression variable list

logit_list <- 
  myfreq %>% filter(!variable %in% c("gender_concept_id","agegr","CCI","CCIGR")) %>% 
  filter(value==1) %>% 
  filter(!(control==0 | patient==0)) %>% 
  pull(variable) %>% unique()

#### ....... 2. Crude OR ######################
# conditional logistic regression: crude OR
# (Table2)

crude_logit <- 
  t1 %>% select(TSH, one_of(logit_list), case_number) %>% 
  # needed variables
  map(~clogit(as.formula("t1$TSH ~ .x + strata(t1$case_number)")) %>% tidy) %>% 
  # logit 
  bind_rows(.,.id="variable") %>% 
  filter(!variable %in% c("TSH","case_number")) %>% 
  mutate(OR = exp(estimate),
         lower_ci = exp(estimate - 1.96*std.error),
         upper_ci = exp(estimate + 1.96*std.error)) %>% 
  # OR, CI
  relocate(variable,term, OR,lower_ci,upper_ci)


write_rds(crude_logit, paste0("output/",HOSPITAL,"-TABLE2.rds"))
write_csv(crude_logit, paste0("output/",HOSPITAL,"-TABLE2.csv"))

#### ....... 3. Adjusted OR ######################
# multivariate conditional logistic regression: adjusted OR
# (Table3)

adj_logit_formula <- as.formula(paste0("TSH~", paste(logit_list, collapse = "+")," +strata(case_number)"))

adjusted_logit <- 
  clogit(adj_logit_formula, t1) %>% 
  tidy %>% 
  # coefficients
  mutate(OR = exp(estimate),
         lower_ci = exp(estimate - 1.96*std.error),
         upper_ci = exp(estimate + 1.96*std.error))
# OR, CI

write_rds(adjusted_logit, paste0("output/",HOSPITAL,"-TABLE3.rds"))
write_csv(adjusted_logit, paste0("output/",HOSPITAL,"-TABLE3.csv"))

##### ... 3. diagram ####
diagram_table <- tribble(~no, ~name , ~number,
                         1, "total patient", diagram_1_entire_patient,
                         2, "TSH patient", diagram_2_TSH,
                         3, "adult", diagram_3_adult,
                         4, "exclusion1: TSH abnormal", diagram_4_exclusion1,
                         5, "exclusion1: condition", diagram_4_exclusion2,
                         6, "exclusion1: Thyroidectomy", diagram_4_exclusion3,
                         7, "exclusion1: drug", diagram_4_exclusion4,
                         8, "exclusion1: I-131 therapy", diagram_4_exclusion5,
                         9, "patient", allo_case2 %>% nrow(),
                         10,"control", allo_control2 %>% nrow())

write.csv(diagram_table, paste0("output/",HOSPITAL,"-diagram.csv"))

#### ... 4. Secondary Analysis ###############################################################################


# ...... 1. Analysis of patients taking alopurinol
 ###################

# List of patients who took Allopurinol
secnd_anal_patient_list <- 
  t1 %>% filter((!is.na(value_source_value)) & (allo == 1)) %>% 
  pull(person_id)

# TSH test records for allopurinol users in the patient group

secnd_anal_patient_TSH <- 
  TSH %>% filter(person_id %in% secnd_anal_patient_list) %>% 
  left_join(allo_case2 %>% select(person_id, indexdt2, case_number), 
            by="person_id") %>% 
  select(person_id, measurement_concept_id, measurement_date, value_source_value, indexdt2, case_number)

# Allopurinol dosing records
secnd_anal_patient_allopurinol <- 
  allo_drug_exposure_list %>% filter(person_id %in% secnd_anal_patient_list) %>% 
  left_join(allo_case2 %>% select(person_id, indexdt2, case_number), 
            by="person_id")


secnd_anal_patient_allopurinol2 <- 
  secnd_anal_patient_allopurinol %>% arrange(person_id, drug_exposure_start_date, drug_exposure_end_date) %>% group_by(person_id) %>% 
  mutate(before_end_date = lag(drug_exposure_end_date),
         drugera = if_else(drug_exposure_start_date>before_end_date+7,1,0),
         drugera2 = if_else(is.na(drugera),1,drugera)) %>% ungroup %>% 
  mutate(drugera3 = cumsum(drugera2)) %>% group_by(person_id, drugera3) %>% 
  summarise(start_date = min(drug_exposure_start_date),
            last_intake = max(drug_exposure_end_date),
            drug_days = last_intake-start_date)

secnd_anal_patient_allopurinol3 <- 
  secnd_anal_patient_allopurinol2 %>% 
  left_join(allo_case2 %>% select(person_id, enroll_date, indexdt2), by='person_id') %>% rename("index_date"="indexdt2") 

# Remove medication history 6 months before or after index dates
secnd_anal_patient_allopurinol4 <- 
  secnd_anal_patient_allopurinol3 %>% 
  filter(last_intake > index_date-180 | (start_date <= index_date & last_intake >= index_date))

# cohort out date 
cohort_out_date <- 
  secnd_anal_patient_allopurinol4 %>% group_by(person_id) %>% 
  filter(last_intake ==max(last_intake)) %>%
  mutate(cohort_out_date = last_intake + days(365))

secnd_anal_patient_allopurinol5 <- 
  secnd_anal_patient_allopurinol4 %>% left_join(cohort_out_date %>% select(person_id, cohort_out_date), by='person_id')

entire_allo_period <- 
  secnd_anal_patient_allopurinol5 %>% group_by(person_id) %>% 
  summarise(entire_start = min(start_date), entire_end = max(last_intake))

secnd_anal_patient_TSH2 <- 
  secnd_anal_patient_TSH %>% select(person_id, measurement_date, value_source_value, indexdt2) %>% 
  left_join(entire_allo_period, by='person_id') %>% 
  filter(measurement_date %within% interval(entire_start, entire_end))

TSH_cohortin <- 
  secnd_anal_patient_TSH %>% group_by(person_id) %>% filter(measurement_date == min(measurement_date)) %>%
  ungroup %>% summarise(value_mean = mean(value_source_value), n=n()) %>% mutate(name = "TSH_cohortin")
TSH_indxdt <- 
  secnd_anal_patient_TSH2 %>% group_by(person_id) %>% filter(measurement_date == min(measurement_date)) %>% ungroup %>% 
  summarise(value_mean = mean(value_source_value), n=n()) %>% mutate(name = "TSH_indexdt")
TSH_after_indxdt_max <- 
  secnd_anal_patient_TSH2 %>% group_by(person_id) %>% filter(value_source_value == max(value_source_value)) %>% ungroup %>% 
  summarise(value_mean = mean(value_source_value), n=n()) %>% mutate(name = "TSH_after_indexdt_max")
TSH_after_indxdt_mean <- 
  secnd_anal_patient_TSH2 %>% summarise(value_mean = mean(value_source_value), n=n())  %>% mutate(name = "TSH_after_indexdt_mean") 


TSH_table <- bind_rows(TSH_cohortin, TSH_indxdt, TSH_after_indxdt_max, TSH_after_indxdt_mean)

write.csv(TSH_table, paste0("output/",HOSPITAL,"-TSH_table.csv"),row.names = FALSE)

### ........... 1. T3 ####
sql <- "SELECT *
        FROM @A.measurement
        WHERE (@B) IN (@C);"
T3_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                  A= myschemaname,
                                  B= source_column,
                                  C= T3code))  


# First Measurement value
T3_1 <- 
  T3_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>%  
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5, by='person_id') %>% 
  filter(measurement_date < start_date) %>% filter(measurement_date==min(measurement_date)) %>% 
  summarise(min_value = min(value_source_value)) %>% ungroup %>% 
  summarise(value = mean(min_value %>% as.numeric()),
            n=n()) %>% mutate(name= "T3_1")



# After index date 
T3_2 <- 
  T3_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>% 
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5 %>% select(person_id, index_date) %>% unique, by='person_id') %>% 
  filter(measurement_date >= index_date) %>% arrange(measurement_date %>% desc) %>% slice(1) %>% 
  mutate(source_value2 = as.numeric(gsub("[^0-9.-]", "", value_source_value))) %>% ungroup() %>%  
  summarise(value = mean(source_value2), n=n()) %>% mutate(name= "T3_2")


### ........... 2. T4 ####
sql <- "SELECT *
        FROM @A.measurement
        WHERE (@B) IN (@C);"
T4_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                  A= myschemaname,
                                  B= source_column,
                                  C= T4code))  

# First Measurement value
T4_1 <-
  T4_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>%  
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5, by='person_id') %>% filter(measurement_date < start_date) %>% 
  filter(measurement_date==min(measurement_date)) %>% 
  summarise(min_value = min(value_source_value)) %>% ungroup %>% 
  summarise(value = mean(min_value %>% as.numeric()), n=n()) %>% mutate(name= "T4_1")

# After index date 
T4_2 <-
  T4_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>% 
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5 %>% select(person_id, index_date) %>% unique, by='person_id') %>% 
  filter(measurement_date >= index_date) %>% arrange(measurement_date %>% desc) %>% slice(1) %>% 
  mutate(source_value2 = as.numeric(gsub("[^0-9.-]", "", value_source_value))) %>% ungroup() %>%  
  summarise(value = mean(source_value2), n=n()) %>% mutate(name= "T4_2")


### ........... 3. freeT4 ####
sql <- "SELECT *
        FROM @A.measurement
        WHERE (@B) IN (@C);"
freeT4_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                      A= myschemaname,
                                      B= source_column,
                                      C= freeT4code)) 

# First Measurement value
freeT4_1 <- 
  freeT4_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>%  
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5, by='person_id') %>% filter(measurement_date < start_date) %>% filter(measurement_date==min(measurement_date)) %>% 
  summarise(min_value = min(value_source_value)) %>% ungroup %>% 
  summarise(value = mean(min_value %>% as.numeric()), n=n())  %>% mutate(name= "freeT4_1")

# After index date 
freeT4_2 <-
  freeT4_list %>% filter(person_id %in% (secnd_anal_patient_allopurinol5 %>% pull(person_id))) %>% as_tibble %>% group_by(person_id) %>% 
  select(person_id, measurement_date, value_source_value) %>% 
  left_join(secnd_anal_patient_allopurinol5 %>% select(person_id, index_date) %>% unique, by='person_id') %>% 
  filter(measurement_date >= index_date) %>% arrange(measurement_date %>% desc) %>% slice(1) %>% 
  mutate(source_value2 = as.numeric(gsub("[^0-9.-]", "", value_source_value))) %>% ungroup() %>%  
  summarise(value = mean(source_value2),n=n()) %>% mutate(name="freeT4_2")

t_series <- list()
t_series[[1]] <- T3_1
t_series[[2]] <- T3_2
t_series[[3]] <- T4_1
t_series[[4]] <- T4_2
t_series[[5]] <- freeT4_1
t_series[[6]] <- freeT4_2

t_series2 <- t_series %>% bind_rows()


write.csv(t_series2, paste0("output/",HOSPITAL,"-t_series.csv"),row.names = FALSE)

# ........... 4. Diseases diagnosed after increased thyroid-stimulating hormone (TSH) levels ####

# allopurinol patients
secnd_anal_patient_disease <- secnd_anal_patient_allopurinol5 %>% pull(person_id) %>% unique

sql <- "SELECT person_id, condition_concept_id, condition_start_date 
        FROM @A.condition_occurrence
        WHERE person_id IN (@B);"
secnd_anal_patient_disease_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype),
                                                          A= myschemaname, 
                                                          B= secnd_anal_patient_disease))

# index_date - before/ same/ after
secnd_anal_patient_disease_list2 <- 
  secnd_anal_patient_disease_list %>% 
  left_join(secnd_anal_patient_allopurinol5 %>% select(person_id, index_date) %>% unique, 
            by="person_id") %>% as_tibble() %>% 
  mutate(condition = if_else(condition_start_date < index_date,"before",
                             if_else(condition_start_date > index_date, "after","same")))

# SNOMED 
secnd_anal_patient_disease_list3 <- 
  secnd_anal_patient_disease_list2 %>% filter(condition=="after") %>% select(person_id, condition_concept_id) %>% 
  anti_join(secnd_anal_patient_disease_list2 %>% filter(condition=="before") %>% select(person_id, condition_concept_id), key=c("person_id","condition_concept_id")) %>% unique

# concept relationship information.
sql <- "select concept_id_1, concept_id_2, relationship_id
        FROM @A.concept_relationship 
        WHERE concept_id_1 in (@B);"
secnd_disease_concept_rela_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                          A= myvocabschemaname, 
                                                          B= secnd_anal_patient_disease_list3 %>% pull(condition_concept_id) %>% unique)) %>% as_tibble()

# ICD10
sql <- "select concept_id, concept_name, domain_id, vocabulary_id, concept_class_id, concept_code
        FROM @A.concept
        WHERE concept_id in (@B);"
secnd_disease_concept_icd_list <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                                         A= myvocabschemaname, 
                                                         B= secnd_disease_concept_rela_list %>% pull(concept_id_2) %>% unique)) %>% as_tibble()


# ICD10CM
sql <- "select concept_id, concept_name, domain_id, vocabulary_id, concept_class_id, concept_code
        FROM @A.concept
        WHERE vocabulary_id = 'ICD10CM';"
icd_code_name <- dbGetQuery(con, render(translate(sql, targetDialect = mydbtype), 
                                        A= myvocabschemaname)) %>% as_tibble() %>% filter(str_length(concept_code)==3)

# newly diagnosed diseases according to ICD10
newly_diagnosed_disease_1 <- 
  secnd_anal_patient_disease_list3 %>% 
  left_join(secnd_disease_concept_rela_list, by=c("condition_concept_id"='concept_id_1')) %>% 
  left_join(secnd_disease_concept_icd_list, by=c("concept_id_2"="concept_id")) %>% 
  filter(vocabulary_id =="ICD10CM" & relationship_id =="Mapped from") %>% 
  group_by(person_id, condition_concept_id) %>% 
  slice(1) %>% mutate(ICD = substr(concept_code,1,3)) %>% 
  group_by(ICD) %>% summarise(n=n()) %>% arrange(n %>% desc) %>% 
  left_join(icd_code_name %>% select(concept_code,concept_name), by=c("ICD"="concept_code"))


write.csv(newly_diagnosed_disease_1, paste0("output/",HOSPITAL,"-newdisease_topN.csv"),row.names = FALSE)

ICDchapters <- 
  tribble( ~Chapter,  ~Block,  ~Title, ~ K_title,
           1,  'A00~B99',	'Certain infectious and parasitic diseases', '특정 감염성 및 기생충성 질환',
           2,  'C00~D48',	'Neoplasms', '신생물',
           3,  'D50~D89',	'Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism', '혈액 및 조혈기관의 질환과 면역 메커니즘을 침범한 특정 장애',
           4,  'E00~E90',	'Endocrine, nutritional and metabolic diseases', '내분비, 영양 및 대사 질환',
           5,  'F00~F99',	'Mental and behavioural disorders', '정신 및 행동 장애',
           6,  'G00~G99',	'Diseases of the nervous system', '신경계통의 질환',
           7,  'H00~H59',	'Diseases of the eye and adnexa', '눈 및 눈 부속기의 질환',
           8,  'H60~H95',	'Diseases of the ear and mastoid process', '귀 및 유돌의 질환',
           9,  'I00~I99', 'Diseases of the circulatory system', '순환계통의 질환',
           10, 'J00~J99',	'Diseases of the respiratory system', '호흡계통의 질환',
           11, 'K00~K93',	'Diseases of the digestive system', '소화계통의 질환',
           12, 'L00~L99', 'Diseases of the skin and subcutaneous tissue', '피부 및 피하조직의 질환',
           13, 'M00~M99',	'Diseases of the musculoskeletal system and connective tissue', '근골격계통 및 결합조직의 질환',
           14, 'N00~N99',	'Diseases of the genitourinary system', '비뇨생식계통의 질환',
           15, 'O00~O99',	'Pregnancy, childbirth and the puerperium', '임신, 출산 및 산후기',
           16, 'P00~P96',	'Certain conditions originating in the perinatal period', '출생전후기에 기원한 특정 병태',
           17, 'Q00~Q99',	'Congenital malformations, deformations and chromosomal abnormalities', '선천기형, 변형 및 염색체 이상',
           18, 'R00~R99',	'Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified', '달리 분류되지 않은 증상, 징후와 임상 및 검사의 이상 소견',
           19, 'S00~T98',	'Injury, poisoning and certain other consequences of external causes', '손상, 중독 및 외인에 의한 특정 기타 결과',
           20, 'V01~Y98',	'External causes of morbidity and mortality', '질병이환 및 사망의 외인',
           21, 'Z00~Z99',	'Factors influencing health status and contact with health services', '건강상태 및 보건서비스 접촉에 영향을 주는 요인',
           22, 'U00~U99',	'Codes for special purposes', '특수목적코드') 


newly_diagnosed_disease_2 <- 
  newly_diagnosed_disease_1 %>% mutate(Chapter = case_when(substr(ICD,1,1) %in% c("A", "B") ~ 1,
                                                           substr(ICD,1,1) == "C" | ICD %in% paste0("D",formatC(seq(1,48),width=2,flag="0")) ~ 2,
                                                           ICD %in% paste0("D",formatC(seq(50,89),width=2,flag="0")) ~ 3,
                                                           substr(ICD,1,1) == "E" ~ 4,
                                                           substr(ICD,1,1) == "F" ~ 5,
                                                           substr(ICD,1,1) == "G" ~ 6,
                                                           ICD %in% paste0("H",formatC(seq(00,59),width=2,flag="0")) ~ 7,
                                                           ICD %in% paste0("H",formatC(seq(60,95),width=2,flag="0")) ~ 8,
                                                           substr(ICD,1,1) == "I" ~ 9,
                                                           substr(ICD,1,1) == "J" ~ 10,
                                                           substr(ICD,1,1) == "K" ~ 11,
                                                           substr(ICD,1,1) == "L" ~ 12,
                                                           substr(ICD,1,1) == "M" ~ 13,
                                                           substr(ICD,1,1) == "N" ~ 14,
                                                           substr(ICD,1,1) == "O" ~ 15,
                                                           substr(ICD,1,1) == "P" ~ 16,
                                                           substr(ICD,1,1) == "Q" ~ 17,
                                                           substr(ICD,1,1) == "R" ~ 18,
                                                           substr(ICD,1,1) %in% c("S", "T") ~ 19,
                                                           substr(ICD,1,1) %in% c("V", "Y") ~ 20,
                                                           substr(ICD,1,1) == "Z" ~ 21,
                                                           substr(ICD,1,1) == "Y" ~ 22)) %>% group_by(Chapter) %>% summarise(n =sum(n)) %>% ungroup %>% 
  left_join(ICDchapters,by='Chapter'); newly_diagnosed_disease_2

write.csv(newly_diagnosed_disease_2, paste0("output/",HOSPITAL,"-newdisease_ICDchapter.csv"), row.names = FALSE)

