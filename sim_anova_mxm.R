#-------------------------------------------------------------------------
# Simulation Exercise 3 by DVM Bishop
# 20th March 2017
# Simulating data for multiway ANOVA
#-------------------------------------------------------------------------
## setwd('/Users/au194693/projects/repro17/anova_mxm_sim') # Here is command to set working directory to Dropbox/BBSCR_STARS/Bishop on my mac

# Remember! you will need to download these packages if they aren't already downloaded
library(reshape2) # for converting data from wide to long form (see below)
library(gridExtra) #for plotting output in a grid
library(grid) #https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html
require(lme4)

options(scipen=999) #disable scientific notation.
set.seed(183)
#-------------------------------------------------------------------------
# We are going to simulate data for a 3 way mixed anova: between/within/within
# Our dependent variable is Score. 
# We have two equal groups (0 and 1), as between subjects variables
# Each subject is measured on two occasions (1 and 2), on easy and hard task (1 and 2)
# This gives 4 scores (easytime1, hardtime1, easytime2, hardtime2) for each subject
# We are just going to generate random numbers - there is no real effect
# We will look for effects of factors group, occasion, and difficulty, 
# as well as interactions between these

# Be warned: as this website notes: http://www.statmethods.net/stats/anova.html
# 'If you have been analyzing ANOVA designs in traditional statistical packages, 
#  you are likely to find R's approach less coherent and user-friendly'
# Our aim here, though, is not to teach you to do ANOVA in R, so much as to 
# reveal some of the pitfalls in interpreting multiway ANOVA
#-------------------------------------------------------------------------

myM <- 0 # Mean score for all variables in the sample - we're using z scores for simplicity
mySD <- 1 #
myN <-30 #set sample size per group (You can vary this to see the effect)
n_sims <- 100000 # Specify number of simulated datasets
# We'll start with simulating 20 datasets, but can later update this number
ptable=matrix(rep(NA, (n_sims*10)), nrow=n_sims) #initialising a matrix that will hold p values in each run
table_names <- c("A", "B", "C", "AB", "AC", "BC", "ABC", "anysig", "Bonfsig", "fdr")
## data.frame for Mixed models
mxm_dt <- as.data.frame(matrix(0, ncol = length(table_names) , nrow = n_sims))
names(mxm_dt) <- table_names

# There are 7 p-values: 3 main effects, 3 2-way interactions, and 1 3-way interaction
# The last two columns will denote if any p-value in a run is <.05, or < .007 (Bonferroni corrected)
colnames(ptable) <- table_names

for (i in 1:n_sims){
    #----------------------------------------------------------------------------------------
    # We'll generate a dataset in 'wide' format, i.e. one row per subject, as this is easy to read
    #----------------------------------------------------------------------------------------
    mydata <-data.frame(matrix(rnorm(n = myN*2*5, mean = myM, sd = mySD),nrow=myN*2))
    mydata[,1]=1 #overwrite column 1 with ones
    mydata[1:myN,1]<-0  #overwrite first half of column 1 with zeros
    colnames(mydata)<-c('Group','Easy1','Hard1','Easy2','Hard2')
    mydata$Group <-as.factor(mydata$Group)  #this is important!
    # You can look at mydata by clicking on its label in the Environment tab


    #----------------------------------------------------------------------------------------
    mylongdata <- melt(mydata) #converts to long form, with one column with all scores
    mylongdata$Subject <- c(seq(1:(myN*2)),seq(1:(myN*2)),seq(1:(myN*2)),seq(1:(myN*2)))
    # Add a column giving subject ID: this is repeated for each combination of conditions
    
    mylongdata$time <-c(rep(1,myN*2*2),rep(2,myN*2*2)) #adds column giving time
    mylongdata$difficulty <-c(rep(1,myN*2),rep(2,myN*2),rep(1,myN*2),rep(2,myN*2))#adds column giving difficulty

    # Note that there are more general ways of reshaping data: this is just simple for this design
    # Look now at longdata by clicking on its name in the Environment tab
    #----------------------------------------------------------------------------------------
    # Save your simulated data as tab-separated text
    #----------------------------------------------------------------------------------------
    write.table(mylongdata,
                paste("/users/mje/anova_sim/repro17/mixed_model_eeg_simulation/data/sim_data/simulated_data_",
                      i, ".csv", sep =""),
                sep=",") 
    #----------------------------------------------------------------------------------------
    # Now run an ANOVA
    #----------------------------------------------------------------------------------------
    
    myaov<-summary(aov(value~(time*difficulty*Group)+Error(Subject/(time*difficulty)),data=mylongdata))
    myaovbit<-unlist(myaov$`Error: Within`) # tortuous way to extract the relevant bit of output
    ptable[i,1:7] <-myaovbit[33:39] #extract the p-values which happen to be values 33-39 in this output
    ptable[i,8:9]<-0 #initialise col 8-9 which will categorise each ANOVA in terms of whether *any* sig effects
    sigp<-which(ptable[i,1:7]<.05) #find whether any p-values are < .05
    if(length(sigp)>0)
    {ptable[i,8]<-1} # if so, assign a 1 to column 8
    sigp<-which(ptable[i,1:7]<.007) #now do the same with p < .007, i.e. .05/7
    if(length(sigp)>0)
    {ptable[i,9]<-1} #result is stored in column 9
    # fdr correction
    fdr_tmp <- which(p.adjust(ptable[i, 1:7], method = "fdr")<.05) #find whether any p-values are < .05
    ifelse(length(fdr_tmp)>0, ptable[i, 10]<-1, ptable[i, 10]<-0)
    
    ## Change to factor for the MM
    mylongdata$Group <- factor(mylongdata$Group)
    mylongdata$time <- factor(mylongdata$time)
    mylongdata$Subject <- factor(mylongdata$Subject)
    ## Fit Mixed model
    m1 <- lmer(value ~ -1 + ( 1 | Subject ), data = mylongdata, REML= FALSE)
    m2 <- update(m1, .~. + time)
    m3 <- update(m2, .~. + difficulty)
    m4 <- update(m3, .~. + Group)
    m5 <- update(m4, .~. + time:difficulty)
    m6 <- update(m5, .~. + time:Group)
    m7 <- update(m6, .~. + difficulty:Group)
    m8 <- update(m7, .~. + time:difficulty:Group)
    
    mxm_dt[i, 1:7] <- anova(m1, m2, m3, m4, m5, m6, m7, m8)$`Pr(>Chisq)`[2:8]
    # fdr correct 
    if (length(sigp)>0)
    {mxm_dt[i,8]<-1} # if so, assign a 1 to column 8
    sigp<-which(mxm_dt[i,1:7]<.007) #now do the same with p < .007, i.e. .05/7
    if(length(sigp)>0)
    {mxm_dt[i,9]<-1} #result is stored in column 9
    
    # fdr correction
    fdr_tmp <- which(p.adjust(mxm_dt[i, 1:7], method = "fdr")<.05) #find whether any p-values are < .05
    if (length(fdr_tmp)>0)
    {mxm_dt[i,8]<-1} # if so, assign a 1 to column 8
    
    
    } #repeat for next simulation
    
    mypf <- round(ptable,digits=2) #get rid of extraneous digits; This doesn't always work- not sure why!
    if (n_sims<21){
        grid.table(mypf) #doing gridtable gets v slow with more than 20 rows
}

# For an explanation of the issues raised by this exercise see:
# http://deevybee.blogspot.co.uk/2013/06/interpreting-unexpected-significant.html

write.csv(ptable, "/users/mje/anova_sim/repro17/mixed_model_eeg_simulation/data/anova_dt.csv")
write.csv(mxm_dt, "/users/mje/anova_sim/repro17/mixed_model_eeg_simulation/data/mxm_dt.csv")

# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. https://creativecommons.org/licenses/by-nc-sa/4.0/
#   

#   
