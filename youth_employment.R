# load libraries
library(oaxaca)
library(readr)
library(dplyr)
library(XML)
library(tidyverse)
library(stargazer)
library(margins)
library(forcats)
library(pROC)
library(VGAM)
library(nnet)


full_codebook <- xmlParse("10852_da10_de_v1_0-ddi.xml")

# Declare the DDI namespace (default in AUSSDA DDI v2 files)
ns <- c(d1 = "http://www.icpsr.umich.edu/DDI")

# Grab all <var> nodes and their variable-level <labl>
vars_nodes <- getNodeSet(full_codebook, "//d1:var", namespaces = ns)

var_names  <- sapply(vars_nodes, function(x) xmlGetAttr(x, "name"))
var_labels <- sapply(vars_nodes, function(x) {
  lab <- getNodeSet(x, "./d1:labl[@level='variable']|./d1:labl", namespaces = ns)
  if (length(lab)) xmlValue(lab[[1]]) else ""
})

full_codebook_df <- data.frame(variable = var_names, label = var_labels, stringsAsFactors = FALSE)
head(full_codebook_df)
nrow(full_codebook_df)

write.csv(full_codebook_df, "full_transcript")
write.csv(full_df, "full_data")

# read the codebook with answers
answers <- read_tsv("10852_vi_de_v2_0.tab")
write.csv(answers, "people_answers")

# read the dataset
df_2024 <- read_tsv("2024_10852_da10_de_v1_0.tab")
df_2023 <- read_tsv("2023_10821_da10_de_v1_0.tab")
df_2022 <- read_tsv("2022_10777_da10_de_v1_0.tab")

combined_df <- bind_rows(
  list(
    "2022" = df_2022,
    "2023" = df_2023,
    "2024" = df_2024
  ),
  .id = "year"   # optional
)

young_df <- combined_df %>%
  filter(balt >= 15 & balt <= 34)
#write_csv(young_df, "young_df.csv")

                     
########################################################################################
#                                     SETUP

# For models. Introduce variables, or rather change and rename them
young_df$bfst<-as.character(young_df$bfst)
young_df$kab11 <-as.character(young_df$kab11)
# Start from your 15–34 sample
young_df_model <- young_df %>%
  mutate(
    # Outcome: employed (ILO)
    employed = as.integer(xerwstat == 1),
    
    # Migration background: 1st or 2nd generation vs natives
    # xmigr_gen: 0 = natives, 1 = 1st gen, 2 = 2nd gen
    migr_bg = as.integer(xmigr_gen %in% c(1, 2)),
    
    # Gender: 1 = female
    female = as.integer(bsex == 2),
    
    # Age and age squared
    age  = balt,
    age2 = balt^2,
    
    # Household structure
    live_parents = as.integer(xeltern1 == 1),  # lives with parents
    hh_size      = bhhgr,                      # household size
    nr_child = xanzkind,
    
    # Region & municipality size
    region_nuts2 = factor(xnuts2),             # NUTS 2 (Bundesland)
    # municipality size class
    
    # Year FE
    year_fe = factor(year),
    
    #Marital status
    marital = (bfst = fct_collapse(bfst, "1" = c("1", "3", "4"))), #4 options (single, married, widowed and divorced)
    
    educ = (kab11=fct_collapse(kab11, "2"=c("2","3","4"))), #8 options, highest completed education
            educ = (kab11=fct_collapse(kab11, "3"=c("5","6")))
            )
# collapsing kab11 into a 3 category variable: compulsury schooling, apprenticeship or matura and higher education

#regions and urbanisation sizes:
young_df_model <- young_df_model %>% mutate(bundesland = case_when(xnuts2==11 ~ "Burgenland",
                                                         xnuts2==12 ~ "Lower Austria",
                                                         xnuts2==13 ~ "Vienna",
                                                         xnuts2==21 ~ "Carinthia",
                                                         xnuts2==22 ~ "Styria",
                                                         xnuts2==31 ~ "Upper Austria",
                                                         xnuts2==32 ~ "Salzburg",
                                                         xnuts2==33 ~ "Tyrol",
                                                         xnuts2==34 ~ "Vorarlberg"),
                                  urb = case_when(xurb==1 ~ "high",
                                                           xurb==2 ~ "medium",
                                                          xurb==3 ~ "low"))

########################################################################################
########################################################################################
#                                     BASE MODEL


#should add marital status, education and
model_emp_all <- glm(
  employed ~ migr_bg + female + age + I(age^2) + educ +
    live_parents + nr_child + marital +
    bundesland +              # regional FE
    urb + 
    year_fe,                   # municipality size FE
  family = binomial(link = "logit"),
  data = young_df_model
)
summary(model_emp_all)



#success - being employed
#the intercept is the base probability of success or the probability of success with all predictors being 0. 
#The slope coefficients are: 1 unit increase in x1, 
#holding all others constant means the log odds of success are increased by beta1
#In order to get the odds, you would need to exponentiate it. 

exp(coef(model_emp_all))


########################################################################################
#                                         TESTS

# Fit the null model (only intercept)
null_model <- glm(employed ~ 1, family = binomial(link = "logit"), data = young_df_model)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_emp_all) / logLik(null_model))) 
#Mcfaddens pseudo R-squared calculates thed difference between just using the intercept to predict (so that every individual is predicted the same probability of 'success')
#vs my model. Sample size for migrants is 18000 vs 31000. 
#migr_bg: For the lighter model its 0.18, for the fuller one its almost 0.2


#Control for overfitting
set.seed(777)
n <- nrow(young_df_model)
train <- sample(seq_len(n), size = floor(0.7 * n))

train <- young_df_model[train, , drop = FALSE]
test  <- young_df_model[setdiff(seq_len(n), train), , drop = FALSE]


model_train <- glm(
  employed ~ migr_bg + female + age + I(age^2) + educ +
    live_parents + nr_child + marital +
    bundesland +              # regional FE
    urb + 
    year_fe,                   # municipality size FE
  family = binomial(link = "logit"),
  data = train
)

# Predicted probabilities
train$pred <- predict(model_train, type = "response") #obtains predictions from a glm model.
test$pred  <- predict(model_train, newdata = test, type = "response") 
#newdata gives the data frame in which it looks for variables to predict. If not given, the fitted predictors are used.
#this means that the predict function chooses whether something is employed or not

print(auc(train$employed, train$pred)) #two vectors, one response and one predictor 
print(auc(test$employed, test$pred)) #the response are the binary variables like 0 and 1. 
#The predictor is the actual model prediction which is like 0.42 or 0.3 and would mean it 
#predicts no employment. If test$employed says it is 0, the model ranked it correctly. 
#If auc is 0.8, 80% of the time, it ranks a positive case above a negative one (which we would want)
#Predict predicts one set of data based on the fitted parameters of the same data. 
#Predict predicts a second set based on the data that the model was not fitted on. 

#actual values from above: 0.7934 vs 0.7938


########################################################################################
#                                         stargazer

stargazer(model_emp_all, )

f<-function(x) exp(x)
stargazer(model_emp_all, apply.coef = f,
          title = "Model1",
          out = "1.3.html",
          align = TRUE)

# Marginal effects or Average partial
margins(model_emp_all, variables = "migr_bg")
print(margins)
########################################################################################
#                                   Model 2:Desired hours
########################################################################################
# dmws=desire to work more(1)/fewer(2)/same(3) hours
emp_df <- young_df_model %>%
  filter(employed == 1, !is.na(dmws))
names(emp_df)

emp_df <- emp_df %>%
  mutate(
    industry_fe = factor(xewzsekt08) # industry fixed effects
  )

model_hours <- multinom(dmws ~ migr_bg + female + age + I(age^2) +
                    educ + live_parents + nr_child + marital
                   + urb + year_fe + industry_fe,
                  data = emp_df)

stargazer(model,
          title = "Model1", apply.coef = f,
          out = "1.3.html",
          align = TRUE)

exp(coef(model))


########################################################################################
#                                   Model 3: Want to work?
########################################################################################
n_emp_df <- young_df_model %>%
  filter(employed == 0, !is.na(dmws))
#sample sizes
#nrow(inactive_df) #23716 vs 32136 if you exclude hawun=NA, vs #77218 that are employed

inactive_df <- young_df_model %>%
  filter(employed == 0, hawun %in% c(1, 2))   # 1 = wants to work, 2 = no. There is also -3, which is a missing value

#changes it to a binary variable

inactive_df <- inactive_df %>%
  mutate(want_work = as.integer(hawun == 1))  # 1 = wants to work
#this just makes it be 0/1. Length still 23761

model_want <- glm(
  want_work ~ migr_bg + female + age + I(age^2) +
    live_parents + nr_child + marital +
    bundesland + urb + year_fe,
  family = binomial("logit"),
  data = inactive_df
)

summary(model_want)
exp(coef(model_want))

########################################################################################
#                                   TESTS
#overfitting
n <- nrow(inactive_df)
train <- sample(seq_len(n), size = floor(0.7 * n))

train <- inactive_df[train, , drop = FALSE]
test  <- inactive_df[setdiff(seq_len(n), train), , drop = FALSE]

model_train <- glm(
  want_work ~ migr_bg + female + age + I(age^2) + educ +
    live_parents + nr_child + marital +
    bundesland +              # regional FE
    urb + 
    year_fe,                   # municipality size FE
  family = binomial(link = "logit"),
  data = train
)

# Predicted probabilities
train$pred <- predict(model_train, type = "response")
test$pred  <- predict(model_train, newdata = test, type = "response")

print(auc(train$want_work, train$pred))
print(auc(test$want_work, test$pred))

null_want <- glm(want_work ~ 1, family = binomial(link = "logit"), data = inactive_df)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_want) / logLik(null_want))) 
#0.086....not great

margins(model_want, variables = "migr_bg")

stargazer(model_want, apply.coef = f,
          title = "Model2",
          out = "1.3.html",
          align = TRUE)

########################################################################################
#                                   Model 4.1: Vienna base
########################################################################################
#Vienna only
vienna<-young_df_model %>% filter(bundesland=="Vienna")
nrow(vienna) #we have a sample size of 16882, instead of 109264, a reduction of 85%
1-nrow(vienna)/nrow(young_df_model)

model_vienna <-glm(employed ~ migr_bg +female + age + I(age^2) + educ +
                     live_parents + nr_child +marital +
                   year_fe , 
                   family = binomial(link="logit"), 
                   data=vienna)

summary(model_vienna)
exp(coef(model_vienna))

########################################################################################
#                                      TESTS

null_vienna <- glm(employed ~ 1, family = binomial(link = "logit"), data = vienna)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_vienna) / logLik(null_vienna))) 
#good R^2 (0.19)

#Control for overfitting
n <- nrow(vienna)
train <- sample(seq_len(n), size = floor(0.7 * n))

train <- vienna[train, , drop = FALSE]
test  <- vienna[setdiff(seq_len(n), train), , drop = FALSE]


model_train <- glm(employed ~ migr_bg +female + age + I(age^2) + educ +
                     live_parents + nr_child +marital +
                     year_fe , 
                   family = binomial(link="logit"), 
                   data=train)

# Predicted probabilities
train$pred <- predict(model_train, type = "response")
test$pred  <- predict(model_train, newdata = test, type = "response")

print(auc(train$employed, train$pred))-print(auc(test$employed, test$pred))


#slight difference, nothing to worry about

########################################################################################
#                               Model 4.2: Vienna desire for more hours
########################################################################################

emp_df2 <- vienna %>%
  filter(employed == 1, !is.na(dmws))

model_hours_vienna <- multinom(dmws ~ migr_bg + female + age + I(age^2) +
                          educ + live_parents + nr_child + marital
                         + year_fe,
                        data = emp_df2)

exp(coef(model_hours_vienna))

n_obs <- nrow(emp_df2)
n_pred <- length(coef(lm(
  dmws ~ migr_bg + female + age + I(age^2) +
    educ + live_parents + nr_child + marital
    + year_fe,
  data = emp_df2)))-1

n_obs / n_pred
n_obs
#we have 10709 observations and about 1070 per predictor(10), which is probably enough.

########################################################################################
#                              Model 4.3: Vienna want to work? - does not work
########################################################################################

inactive_df <- vienna %>%
  filter(employed == 0, hawun %in% c(1, 2))   # 1 = wants to work, 2 = no. There is also -3, which is a missing value

#unemployed but want to work
inactive_df <- inactive_df %>%
  mutate(want_work = as.integer(hawun == 1))

model_want_vienna <- glm(
  want_work ~ migr_bg + female + age + I(age^2) +
    live_parents + nr_child + marital + year_fe,
  family = binomial("logit"),
  data = inactive_df
)

exp(coef(model_want_vienna))

########################################################################################
#                               Model 5.1: Origin base - comparison
########################################################################################
#                                          Setup
#combining AUT and GER
young_df_model <- young_df_model %>%
  filter(xbstaa16!="16")
young_df_model$xbstaa16<-as.character(young_df_model$xbstaa16)

young_df_origin <- young_df_model %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "0" = c("1", "3")))

young_df_origin <- young_df_origin %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "1" = c("2", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")))

summary(young_df_model$xbstaa16)
table(young_df_origin$xbstaa16)
table(young_df_model$xbzzgla16) #Zuzugsland 16 Gruppen f_bzzgland


                     
########################################################################################
#                                       MODEL

model_young_origin_emp <-glm(employed ~ xbstaa16 +female + age + I(age^2) + educ +
                     live_parents + nr_child +marital +
                     year_fe, 
                   family = binomial(link="logit"), 
                   data=young_df_origin)

summary(model_young_origin_emp) #-0.63 instead of -0.53 if we used migr. background
exp(coef(model_young_origin_emp))

########################################################################################
#                                     TESTS

null_origin <- glm(employed ~ 1, family = binomial(link = "logit"), data = young_df_origin)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_young_origin_emp) / logLik(null_origin))) 
#0.194


#overfitting
n <- nrow(young_df_origin)
train <- sample(seq_len(n), size = floor(0.8 * n))

train <- young_df_origin[train, , drop = FALSE]
test  <- young_df_origin[setdiff(seq_len(n), train), , drop = FALSE]


model_train <- glm(employed ~ xbstaa16 +female + age + I(age^2) + educ +
                        live_parents + nr_child +marital +
                        year_fe, 
                      family = binomial(link="logit"), 
                      data=train)

# Predicted probabilities
train$pred <- predict(model_train, type = "response")
test$pred  <- predict(model_train, newdata = test, type = "response")

auc_train <- auc(train$employed, train$pred)
auc_test  <- auc(test$employed, test$pred)

auc_train
auc_test
#0.7845 vs 0.7854



########################################################################################
#                   Model 5.2: Adding western europe and the EU
########################################################################################
#                                     SETUP
#Western europe -eu and EFTA/UK vs the rest

young_df_origin <- young_df_model %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "0" = c("1", "3", "2", "4", "5", "6", "7", "8")))

young_df_origin <- young_df_origin %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "1" = c("9", "10", "11", "12", "13", "14", "15")))



########################################################################################
#                                     MODEL

model_young_west_emp <-glm(employed ~ xbstaa16 +female + age + I(age^2) + educ +
                               live_parents + nr_child +marital +
                               year_fe, 
                             family = binomial(link="logit"), 
                             data=young_df_origin)
summary(model_young_west_emp) #-0.85, crazy
exp(coef(model_young_west_emp))

########################################################################################
#                                     TESTS


null_origin <- glm(employed ~ 1, family = binomial(link = "logit"), data = young_df_origin)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_young_west_emp) / logLik(null_origin))) 
#0.195
table(young_df_origin$xbstaa16)



#overfitting
n <- nrow(young_df_origin)
train <- sample(seq_len(n), size = floor(0.7 * n))

train <- young_df_origin[train, , drop = FALSE]
test  <- young_df_origin[setdiff(seq_len(n), train), , drop = FALSE]


model_train <- glm(employed ~ migr_bg +female + age + I(age^2) + educ +
                     live_parents + nr_child +marital +
                     year_fe , 
                   family = binomial(link="logit"), 
                   data=train)

# Predicted probabilities
train$pred <- predict(model_train, type = "response")
test$pred  <- predict(model_train, newdata = test, type = "response")

print(auc(train$employed, train$pred))
print(auc(test$employed, test$pred))

#0.7867 vs 0.7866









##########################################################################################
#                                     PROBABLY UNNECESSARY
##########################################################################################
#Adding the rest of european countries and america (11 and 13)
                     
young_df_origin <- young_df_model %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "0" = c("1","2",  "3", "4", "5", "6", "7", "8","11", "13")))

young_df_origin <- young_df_origin %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "1" = c("9", "10", "12", "14", "15")))

model_young_expand_emp <-glm(employed ~ xbstaa16 +female + age + I(age^2) + educ +
                             live_parents + nr_child +marital +
                             year_fe, 
                           family = binomial(link="logit"), 
                           data=young_df_origin)
summary(model_young_expand_emp)

###############################################

origin_df <- combined_df %>%
  filter(balt >= 15 & balt <= 65)

nrow(young_df)-nrow(origin_df) #413042 more obs

origin_df <- origin_df %>%
  mutate(
    # Outcome: employed (ILO)
    employed = as.integer(xerwstat == 1),
    
    # Migration background: 1st or 2nd generation vs natives
    # xmigr_gen: 0 = natives, 1 = 1st gen, 2 = 2nd gen
    migr_bg = as.integer(xmigr_gen %in% c(1, 2)),
    
    # Gender: 1 = female
    female = as.integer(bsex == 2),
    
    # Age and age squared
    age  = balt,
    age2 = balt^2,
    
    # Household structure
    live_parents = as.integer(xeltern1 == 1),  # lives with parents
    hh_size      = bhhgr,                      # household size
    
    # Region & municipality size
    region_nuts2 = factor(xnuts2),             # NUTS 2 (Bundesland)
    mun    = factor(xeinw),               # municipality size class
    
    # Year FE
    year_fe = factor(year),
    
    educ = kab11, #8 options, highest completed education
  )

#changing municipality size and marital status
origin_df <- origin_df %>% mutate(urb = case_when(xurb==1 ~ "high",
                                                                     xurb==2 ~ "medium",
                                                                     xurb==3 ~ "low"))
origin_df$bfst<-as.character(origin_df$bfst)
origin_df <- origin_df%>%
  mutate(marital = fct_collapse(bfst, "1" = c("1", "3", "4")))
#single, divorced and widowed all go into single(1)


origin_df <- origin_df %>%
filter(xbstaa16!="16")
table(origin_df$xbstaa16)


origin_df$xbstaa16<-as.character(origin_df$xbstaa16)

origin_old<- origin_df %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "0" = c("1", "3")))

origin_old<- origin_old %>%
  mutate(xbstaa16 = fct_collapse(xbstaa16, "1" = c("2", "4", "5", "6", "7", "8", "11", "13", "9", "10", "12", "14", "15")))


model_old_expand_emp <-glm(employed ~ migr_bg +female + age + I(age^2) + educ +
                               live_parents + nr_child +marital + urb +
                               year_fe, 
                             family = binomial(link="logit"), 
                             data=origin_old)
summary(model_old_expand_emp) #-1.383 pretty good


null_origin <- glm(employed ~ 1, family = binomial(link = "logit"), data = origin_old)

# Calculate McFadden's Pseudo R-squared
print(1 - (logLik(model_old_expand_emp) / logLik(null_origin))) 
#0.2286
#even here, the pseudo R^2 is negative for clustering germans and austrians but positive for migrational background. 
table(origin_old$xbstaa16)


########################################################################################
#                                   OAXACA DECOMPOSITION
########################################################################################

library(oaxaca)
library(dplyr)
library(tidyr)


make_group01 <- function(df) df %>% mutate(group = as.integer(migr_bg == 1))
make_binary01 <- function(x) as.integer(x == 1)

run_oaxaca_auto_blocks <- function(data, formula, blocks_fun, R = 999, robust = TRUE) {
  
  fit <- oaxaca(
    formula = formula,
    data    = data,
    R       = R,
    robust  = robust
  )
  
  var_tab <- as.data.frame(fit$twofold$variables[[3]])
  
  overall <- c(
    gap         = fit$y$y.diff,
    explained   = sum(var_tab[["coef(explained)"]],   na.rm = TRUE),
    unexplained = sum(var_tab[["coef(unexplained)"]], na.rm = TRUE)
  )
  
  blocks <- blocks_fun(var_tab)
  
  extract_block <- function(rows) {
    if (length(rows) == 0) return(c(explained = NA_real_, unexplained = NA_real_))
    rows <- rows[rows %in% rownames(var_tab)]
    if (length(rows) == 0) return(c(explained = 0, unexplained = 0))
    c(
      explained   = sum(var_tab[rows, "coef(explained)"],   na.rm = TRUE),
      unexplained = sum(var_tab[rows, "coef(unexplained)"], na.rm = TRUE)
    )
  }
  
  block_tab <- do.call(rbind, lapply(blocks, extract_block))
  block_tab <- as.data.frame(block_tab)
  block_tab$block <- rownames(block_tab)
  rownames(block_tab) <- NULL
  block_tab <- block_tab[, c("block","explained","unexplained")]
  
  list(
    fit       = fit,
    overall   = overall,
    variables = var_tab,
    blocks    = block_tab
  )
}

make_blocks_full_updated2 <- function(var_tab, educ_name = "educ") {
  rn <- rownames(var_tab)
  educ_pat <- paste0("^factor\\(", educ_name, "\\)")
  list(
    female   = "female",
    age      = c("age", "age2"),
    educ     = c(educ_name, grep(educ_pat, rn, value = TRUE)),
    parents  = "live_parents",
    children = "nr_child",
    marital  = grep("^factor\\(marital\\)", rn, value = TRUE),
    region   = grep("^factor\\(bundesland\\)", rn, value = TRUE),
    urb      = grep("^factor\\(urb\\)", rn, value = TRUE),
    year     = grep("^factor\\(year\\)", rn, value = TRUE)
  )
}

make_blocks_vienna_updated2 <- function(var_tab, educ_name = "educ") {
  rn <- rownames(var_tab)
  educ_pat <- paste0("^factor\\(", educ_name, "\\)")
  list(
    female   = "female",
    age      = c("age", "age2"),
    educ     = c(educ_name, grep(educ_pat, rn, value = TRUE)),
    parents  = "live_parents",
    children = "nr_child",
    marital  = grep("^factor\\(marital\\)", rn, value = TRUE),
    year     = grep("^factor\\(year\\)", rn, value = TRUE)
  )
}


pick_educ_name <- function(df) if ("educ3" %in% names(df)) "educ3" else "educ"


########################################################################################
# 1) EMPLOYMENT — FULL SAMPLE (young_df_model, group = migr_bg)
########################################################################################
educ_name <- pick_educ_name(young_df_model)

df_emp <- young_df_model %>%
  make_group01() %>%
  mutate(
    employed     = as.integer(employed == 1),
    age2         = age^2,
    female       = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    # IMPORTANT: use your actual FE vars here
    marital    = factor(marital),
    bundesland = factor(bundesland),
    urb        = factor(urb),
    year       = factor(year_fe)   # your FE variable is year_fe
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name]], live_parents, nr_child,
          marital, bundesland, urb, year)

form_emp <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name, ") + ",
    "live_parents + nr_child + factor(marital) + factor(bundesland) + ",
    "factor(urb) + factor(year) | group"
  )
)

oax_emp <- run_oaxaca_auto_blocks(
  data = df_emp,
  formula = form_emp,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name),
  R = 999, robust = TRUE
)

oax_emp$overall
oax_emp$blocks
# deep dive:
oax_emp$variables


########################################################################################
# 2) DESIRED HOURS (EMPLOYED ONLY) — FULL SAMPLE
########################################################################################
df_hours <- young_df_model %>%
  filter(employed == 1, !is.na(dmws)) %>%
  make_group01() %>%
  mutate(
    dmws_y = as.numeric(dmws),
    age2   = age^2,
    female = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital    = factor(marital),
    bundesland = factor(bundesland),
    urb        = factor(urb),
    year       = factor(year_fe)
  ) %>%
  drop_na(dmws_y, group, female, age, age2, .data[[educ_name]], live_parents, nr_child,
          marital, bundesland, urb, year)

form_hours <- as.formula(
  paste0(
    "dmws_y ~ female + age + age2 + factor(", educ_name, ") + ",
    "live_parents + nr_child + factor(marital) + factor(bundesland) + ",
    "factor(urb) + factor(year) | group"
  )
)

oax_hours <- run_oaxaca_auto_blocks(
  data = df_hours,
  formula = form_hours,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name),
  R = 999, robust = TRUE
)

oax_hours$overall
oax_hours$blocks
oax_hours$variables


########################################################################################
# 3) WANT TO WORK (NON-EMPLOYED) — FULL SAMPLE
########################################################################################
df_want <- young_df_model %>%
  filter(employed == 0, hawun %in% c(1,2)) %>%
  make_group01() %>%
  mutate(
    want_work = as.integer(hawun == 1),
    age2      = age^2,
    female    = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital    = factor(marital),
    bundesland = factor(bundesland),
    urb        = factor(urb),
    year       = factor(year_fe)
  ) %>%
  drop_na(want_work, group, female, age, age2, .data[[educ_name]], live_parents, nr_child,
          marital, bundesland, urb, year)

form_want <- as.formula(
  paste0(
    "want_work ~ female + age + age2 + factor(", educ_name, ") + ",
    "live_parents + nr_child + factor(marital) + factor(bundesland) + ",
    "factor(urb) + factor(year) | group"
  )
)

oax_want <- run_oaxaca_auto_blocks(
  data = df_want,
  formula = form_want,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name),
  R = 999, robust = TRUE
)

oax_want$overall
oax_want$blocks
oax_want$variables


########################################################################################
# 4) VIENNA — EMPLOYMENT
########################################################################################
educ_name_vie <- pick_educ_name(vienna)

df_emp_vie <- vienna %>%
  make_group01() %>%
  mutate(
    employed     = as.integer(employed == 1),
    age2         = age^2,
    female       = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name_vie]], live_parents, nr_child,
          marital, year)

form_emp_vie <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name_vie, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_emp_vie <- run_oaxaca_auto_blocks(
  data = df_emp_vie,
  formula = form_emp_vie,
  blocks_fun = function(vt) make_blocks_vienna_updated2(vt, educ_name_vie),
  R = 999, robust = TRUE
)

oax_emp_vie$overall
oax_emp_vie$blocks
oax_emp_vie$variables


########################################################################################
# 5) VIENNA — DESIRED HOURS (EMPLOYED ONLY)
########################################################################################
df_hours_vie <- vienna %>%
  filter(employed == 1, !is.na(dmws)) %>%
  make_group01() %>%
  mutate(
    dmws_y = as.numeric(dmws),
    age2   = age^2,
    female = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(dmws_y, group, female, age, age2, .data[[educ_name_vie]], live_parents, nr_child,
          marital, year)

form_hours_vie <- as.formula(
  paste0(
    "dmws_y ~ female + age + age2 + factor(", educ_name_vie, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_hours_vie <- run_oaxaca_auto_blocks(
  data = df_hours_vie,
  formula = form_hours_vie,
  blocks_fun = function(vt) make_blocks_vienna_updated2(vt, educ_name_vie),
  R = 999, robust = TRUE
)

oax_hours_vie$overall
oax_hours_vie$blocks
oax_hours_vie$variables


########################################################################################
# 6) VIENNA — WANT TO WORK (NON-EMPLOYED)
########################################################################################
df_want_vie <- vienna %>%
  filter(employed == 0, hawun %in% c(1,2)) %>%
  make_group01() %>%
  mutate(
    want_work = as.integer(hawun == 1),
    age2      = age^2,
    female    = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(want_work, group, female, age, age2, .data[[educ_name_vie]], live_parents, nr_child,
          marital, year)

form_want_vie <- as.formula(
  paste0(
    "want_work ~ female + age + age2 + factor(", educ_name_vie, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_want_vie <- run_oaxaca_auto_blocks(
  data = df_want_vie,
  formula = form_want_vie,
  blocks_fun = function(vt) make_blocks_vienna_updated2(vt, educ_name_vie),
  R = 999, robust = TRUE
)

oax_want_vie$overall
oax_want_vie$blocks
oax_want_vie$variables


########################################################################################
# 7) ORIGIN SPLIT (young_df_origin) — EMPLOYMENT
########################################################################################
educ_name_org <- pick_educ_name(young_df_origin)

df_origin_emp <- young_df_origin %>%
  mutate(
    group    = as.integer(xbstaa16 == "1"),
    employed = as.integer(employed == 1),
    age2     = age^2,
    female   = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name_org]], live_parents, nr_child,
          marital, year)

form_origin_emp <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name_org, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_origin_emp <- run_oaxaca_auto_blocks(
  data = df_origin_emp,
  formula = form_origin_emp,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name_org),
  R = 999, robust = TRUE
)

oax_origin_emp$overall
oax_origin_emp$blocks
oax_origin_emp$variables


########################################################################################
# 8) WESTERN SPLIT (your model_young_west_emp uses young_df_origin after re-collapse)
########################################################################################

educ_name_west <- pick_educ_name(young_df_origin)

df_west_emp <- young_df_origin %>%
  mutate(
    group    = as.integer(xbstaa16 == "1"),
    employed = as.integer(employed == 1),
    age2     = age^2,
    female   = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name_west]], live_parents, nr_child,
          marital, year)

form_west_emp <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name_west, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_west_emp <- run_oaxaca_auto_blocks(
  data = df_west_emp,
  formula = form_west_emp,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name_west),
  R = 999, robust = TRUE
)

oax_west_emp$overall
oax_west_emp$blocks
oax_west_emp$variables


########################################################################################
# 9) EXPANDED SPLIT (your model_young_expand_emp uses young_df_origin after re-collapse)
########################################################################################

educ_name_exp <- pick_educ_name(young_df_origin)

df_expand_emp <- young_df_origin %>%
  mutate(
    group    = as.integer(xbstaa16 == "1"),
    employed = as.integer(employed == 1),
    age2     = age^2,
    female   = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    year    = factor(year_fe)
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name_exp]], live_parents, nr_child,
          marital, year)

form_expand_emp <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name_exp, ") + ",
    "live_parents + nr_child + factor(marital) + factor(year) | group"
  )
)

oax_expand_emp <- run_oaxaca_auto_blocks(
  data = df_expand_emp,
  formula = form_expand_emp,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name_exp),
  R = 999, robust = TRUE
)

oax_expand_emp$overall
oax_expand_emp$blocks
oax_expand_emp$variables


########################################################################################
# 10) OLDER SAMPLE (origin_old) — EMPLOYMENT (optional)
########################################################################################
educ_name_old <- pick_educ_name(origin_old)

df_old_emp <- origin_old %>%
  make_group01() %>%   # group = migr_bg
  mutate(
    employed     = as.integer(employed == 1),
    age2         = age^2,
    female       = as.integer(female == 1),
    live_parents = as.integer(live_parents == 1),
    
    marital = factor(marital),
    urb     = factor(urb),
    year    = factor(year_fe)
  ) %>%
  drop_na(employed, group, female, age, age2, .data[[educ_name_old]], live_parents, nr_child,
          marital, urb, year)

form_old_emp <- as.formula(
  paste0(
    "employed ~ female + age + age2 + factor(", educ_name_old, ") + ",
    "live_parents + nr_child + factor(marital) + factor(urb) + factor(year) | group"
  )
)

oax_old_emp <- run_oaxaca_auto_blocks(
  data = df_old_emp,
  formula = form_old_emp,
  blocks_fun = function(vt) make_blocks_full_updated2(vt, educ_name_old),
  R = 999, robust = TRUE
)

oax_old_emp$overall
oax_old_emp$blocks
oax_old_emp$variables


########################################################################################
# 11) Collect everything in one list (easy printing)
########################################################################################
oaxaca_results <- list(
  emp_all      = oax_emp,
  hours_all    = oax_hours,
  want_all     = oax_want,
  emp_vienna   = oax_emp_vie,
  hours_vienna = oax_hours_vie,
  want_vienna  = oax_want_vie,
  origin_emp   = oax_origin_emp,
  west_emp     = oax_west_emp,
  expand_emp   = oax_expand_emp,
  old_emp      = oax_old_emp
)

# Quick view: block tables
lapply(oaxaca_results, \(x) x$blocks)

# Deep dive example:
# oaxaca_results$emp_all$variables
########################################################################################
                     
oaxaca(
  employed ~ female + educ + age + I(age^2) + live_parents + nr_child + marital + bundesland + year | migr_bg,
  data = young_df_model,
  model = "logit",
  R = 1000
)

?oaxaca

oaxaca(
  employed ~ educ + age | migr_bg,
  data = young_df_model,
  model = "logit",
  R = 1000
)

oaxaca(
  dmws ~ migr_bg + female + age + 
    educ + live_parents + nr_child + marital
  + bundesland + muni_size + year,
  data = emp_df,
  group.weight = 0   # 0 = group 0 as reference
)


model_hours <- lm(
  dmws ~ migr_bg + female + age + I(age^2) +
    educ + live_parents + nr_child + marital
  + bundesland + muni_size + year,
  data = emp_df
)


model_emp_all <- glm(
  employed ~ migr_bg + female + age + I(age^2) + educ +
    live_parents + nr_child + marital +
    bundesland +              # regional FE
    muni_size + 
    year,                   # municipality size FE
  family = binomial(link = "logit"),
  data = young_df_model
)


m0 <- lm(dmws ~ female + educ + age + I(age^2) + live_parents + nr_child + marital + bundesland + year,   data = young_df_model, subset = migr_bg == 0)
m1 <- lm(dmws ~ female + educ + age + I(age^2) + live_parents + nr_child + marital + bundesland + year,   data = young_df_model, subset = migr_bg == 1)

X0 <- colMeans(model.matrix(m0))
X1 <- colMeans(model.matrix(m1))

model.matrix(m0)

#forms row and column sums and means

beta0 <- coef(m0)
beta1 <- coef(m1)

explained <- (X1 - X0) %*% beta0
unexplained <- X1 %*% (beta1 - beta0)

explained
unexplained
