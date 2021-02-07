#if we're in the dashboard, use the information from there, otherwise upload locally for development purposes

if(!is.element("SHINY",ls())){
  #Define the settings matrix
  Settings<-data.frame("Description"=
                   c("Assay validation level",
                     "Total number of cycles used in qPCR run",
                     "Minimum number of positive qPCR replicates for DNA sample to be considered positive (default: 2)",
                     "Cq value for Limit of Detection (default: determined from standard curve [0])",
                     "Inhibition Threshold, cycle delay relative to pcrncs (default: 1)",
                     "Override Inhibition Threshold, just use minimum cycle number that IPC should have attained (default:0 [OFF])",
                     "IPC threshold: minimum number of replicates that must have been tested for DNA sample to be considered conclusively negative (default:12)",
                     "IPC threshold2: maximum number of replicates that exhibited inhibition for DNA sample to be considered conclusively negative (default:0)"
                   ),
                 "Setting"=c(
                   2,
                   45,
                   1,
                   30,
                   1,
                   37,
                   12,
                   4
                 ),
                 "shortname"=c(
                   "AVL",
                   "N_cycles",
                   "Pos_minreps",
                   "LOD_Cq",
                   "IPC_delay",
                   "IPC_Cq",
                   "IPC_minreps",
                   "IPC_Inh_maxreps"
                 )
                 )
  #Read the input data (which would normally be uploaded on the dashboard)
  DF<-read.csv("Data.csv")
}

#Column names for the table in the pdf - note I have added the last column header, which will contain the results which are to be displayed
#Important that the name "colNames" is not changed as the shiny code uses it
colNames<-c("Plate","Well","sample type","DNA Sample","Replicate","Target Cq","IPC Cq","std conc","Sampling point","Extraction Batch","Volume Water Processed","Interpretation")

#Some random data manipulation for illustration purposes. The vector Int contains the results, and will be appended to the dataframe and output in the pdf table
#We did some googling and it looks like the ifelse function in R can be used to copy excel if statements (more or less) directly
#I'm not sure if you're familiar with ifelse, but it looks like it could save you some time!
# Int<-array(NA,dim=dim(DF)[1])
# if(Override==0){
#   w<-which(!is.na(DF$IPC_Cq)&DF$IPC_Cq>35)
#   Int[w]<-"FAILED"
#   w<-which(!is.na(DF$IPC_Cq)&DF$IPC_Cq<=35)
#   Int[w]<-"PASSED"
# }else if(Override==1){
#   w<-which(!is.na(DF$IPC_Cq)&DF$IPC_Cq>35)
#   Int[w]<-"RANDOM"
#   w<-which(!is.na(DF$IPC_Cq)&DF$IPC_Cq<=35)
#   Int[w]<-"TEST"
# }

#TODO
#add setting on whether to use LOD at all
#add option to use concentration rather than cq value for LOD
#add fields to write validation of non-default settings
#if stds on plate, add plot stds, and add all other samples to plot
#possibly to be interactive so that can hover over samples to show ids? (guess not for pdf)
#add readme to github

##############################################

#step1 - formatting cols
DF$sample_type<-as.character(DF$sample_type)
DF$DNA_Sample<-as.character(DF$DNA_Sample)
DF$Replicate<-as.character(DF$Replicate)
DF$Target_Cq<-as.numeric(DF$Target_Cq)
DF$IPC_Cq<-as.numeric(DF$IPC_Cq)
DF$std_conc<-as.numeric(DF$std_conc)
#assume all NAs in Target_Cq and IPC_Cq are 0
DF$Target_Cq[is.na(DF$Target_Cq)]<-0
DF$IPC_Cq[is.na(DF$IPC_Cq)]<-0

#remove NA rows
DF<-DF[DF$DNA_Sample!="",]

#only allow "std","pcrnc","extnc","unkn", "fieldnc","pc" in DF$sample_type
if(sum(
      DF$sample_type == "std" | DF$sample_type =="pcrnc" | DF$sample_type == "extnc" |
      DF$sample_type == "unkn" | DF$sample_type == "fieldnc" | DF$sample_type == "pc"
      )!=nrow(DF)
   ) stop("sample_type must be one of: std, pcrnc, extnc, unkn, fieldnc, pc")

#split DF by sample type for later
DFlist<-split(DF,DF$sample_type)

#small settings function to get setting from Settings, based on shortname
getset<-function(Settingsdf,setting_shortname){
  Settingsdf[match(setting_shortname,table = Settingsdf$shortname),"Setting"]
}

##############################################
##qPCR replicate level report

#IPC filter
DF$IPC_threshold<-"NOT_DONE"
DF$IPC_threshold[which(DF$IPC_Cq==0)]<-"FAILED" #fail IPC if IPC_Cq was 0
DF$IPC_threshold[which(DF$IPC_Cq==-1)]<-"IPC NOT USED" #IPC not used if IPC_Cq was 0

if(getset(Settings,"IPC_Cq")==0){ #using the IPC delay approach
  meanNegIPC<-mean(DFlist$pcrnc$IPC_Cq[which(DFlist$pcrnc$IPC_Cq>-1)]) #get average of IPC for pcr negs (where used)
  #where not so far categorised, check that the IPC delay was not exceeded
  DF$IPC_threshold[which(DF$IPC_Cq-meanNegIPC>getset(Settings,"IPC_delay") & DF$IPC_threshold=="NOT_DONE")]<-"FAILED"
} else {#using the override
  DF$IPC_threshold[which(DF$IPC_Cq>getset(Settings,"IPC_Cq") & DF$IPC_threshold=="NOT_DONE")]<-"FAILED"
}
#everything else is passed
DF$IPC_threshold[DF$IPC_threshold=="NOT_DONE"]<-"PASSED"

#LOD filter
DF$LOD_threshold<-"NOT_DONE"
DF$LOD_threshold[which(DF$Target_Cq==0)]<-"NO DETECTABLE CURVE"

#getting LOD
if(getset(Settings,"LOD_Cq")==0){
  #get LOD from stds
  stds<-DFlist$std
  stdspos<-stds[stds$Target_Cq>0,]
  stdspos_agg<-aggregate(stdspos$Target_Cq,by = list(stdspos$std_conc),FUN=length)
  #Choose lowest conc that had n-1 reps working - why? Using a 95% rule, 
  #but leniently 1/2, 2/3, 3/4, 5/6, 6/7, 7/8, 8/9,9/10,10/11... 
  #note this will break down at >19 replicates...not very likely!
  nstdreps_m1<-max(stdspos_agg$x)-1
  stdspos_agg$threshold<-"NOT_DONE"
  stdspos_agg$threshold[which(stdspos_agg$x<nstdreps_m1)]<-"FAIL"
  stdspos_agg$threshold[stdspos_agg$threshold=="NOT_DONE"]<-"PASS"
  min_conc_pass<-min(stdspos_agg$Group.1[which(stdspos_agg$threshold=="PASS")])
  
  #if lowest concs passed, cannot use this LOD, revert to cycle number
  if(min(stds$std_conc)==min_conc_pass) {
    message("Limit of detection not reached (all stds had >=95% detection success, 
            not applying LOD threshold, suggest applying manual LOD value")
    LOD_Cq<-getset(Settings,"N_cycles")
  } else LOD_Cq<-max(stdspos[stdspos$std_conc==min_conc_pass,"Target_Cq"])
} else {
  LOD_Cq<-getset(Settings,"LOD_Cq") #absolute value from input
}

DF$LOD_threshold[which(DF$Target_Cq>LOD_Cq)]<-"FAILED" 
DF$LOD_threshold[which(DF$LOD_threshold=="NOT_DONE")]<-"PASSED"

#interpretations


#=IF(M90="PASSED",IF(L90="PASSED","positive",IF(L90="IPC NOT USED","positive","positive")),IF(M90="NO DETECTABLE CURVE",IF(L90="PASSED","negative",IF(L90="IPC NOT USED","not detected, inhibition not tested","not detected, PCR inhibited")),IF(L90="PASSED","weak signal below LOD, no inhibition observed",IF(L90="IPC NOT USED","weak signal below LOD, inhibition not tested","weak signal below LOD, inhibition observed"))))

##############################################
#DNA sample level report

##############################################
#This DF definition is important - this is what shiny will take for the pdf output, along with the colNames vector for column headers
#Also important that the name "DF" is not changed
DF<-cbind(DF,data.frame("Interpretation"=Int))
##############################################


