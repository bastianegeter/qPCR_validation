#if we're in the dashboard, use the information from there, otherwise upload locally for development purposes

if(!is.element("SHINY",ls())){
  #Define the settings matrix
  Settings<-data.frame("shortname"=c(
    "N_cycles",
    "Pos_minreps",
    "LOD_Cq",
    "IPC_delay",
    "IPC_Cq",
    "IPC_minreps",
    "IPC_Inh_maxreps"
  ),"Description"=c("Total number of cycles used in qPCR run",
                    "Minimum number of positive qPCR replicates for DNA sample to be considered positive (default: 2)",
                    "Cq value for Limit of Detection (default: determined from standard curve [0])",
                    "Inhibition cycle delay relative to pcrncs (default: 1)",
                    "Override cycle delay, just use minimum cycle number that IPC should have attained (default:0 [OFF])",
                    "Minimum number of IPC replicates that must have been tested for a DNA sample to be considered conclusively negative (default:12)",
                    "Maximum number of IPC replicates that may have been inhibited for a DNA sample to be considered conclusively negative (default:0)"
  ), "Setting"=c(
    55,
    1,
    30,
    1,
    37,
    12,
    4
  ))
  #Read the input data (which would normally be uploaded on the dashboard)
  DF<-read.csv("Data.csv")
  
  #Apply validation scale settings (should come from dashboard also)
  validation<-data.frame(Metric=c(
    "Was in silico testing conducted and primers shown to amplify the target species?",
    "Were the primers tested on tissue from the target species?",
    "Was in silico testing conducted and potential crossamplification of non-target species shown to be low?",
    "Were primers tested on non-target tissue of closely-related potentially co-occuring species?",
    "Did the assay successfully detect the target species at a site of known presence?",
    "Did the assay successfully detect the target species at multiple site of known presence?",
    "Did the assay return negative results for the target species at multiple sites of known absence?",
    "Has assay sensitivity (Limit of Detection and/or Limit of Quantification) been assessed?",
    "Has site occupancy modelling (or equivalent) been conducted?",
    "Has the probability of detecting a target species at a site been calculated?",
    "Has the number of water samples needed to achieve reliable detection from a site been calculated?",
    "Has the number of water samples needed to estimate probability of species absence given negative results from a site been calculated?",
    "Has the number of qPCR replicates needed to achieve reliable detection in an eDNA sample been calculated?"
  ), Level_of_Confidence=c(
    "Low",
    "Low",
    "Low",
    "Low",
    "Low",
    "Medium",
    "Medium",
    "Medium",
    "High",
    "High",
    "High",
    "High",
    "High"
  ), user_input=c(1,1,1,1,1,1,1,1,1,1,0,0,0)
  )
}

#Column names for the table in the pdf - note I have added the last column header, which will contain the results which are to be displayed
#Important that the name "colNames" is not changed as the shiny code uses it
colNames<-c("Plate","Well","Sample_Type","DNA_Sample","Replicate","Target_Cq","IPC_Cq",
            "Std_Conc","Sampling_Point","Extraction_Batch","Volume_Water_Processed","Interpretation")

#TODO
#add setting on whether to use LOD at all
#add option to use concentration rather than cq value for LOD
#add fields to write validation of non-default settings
#if stds on plate, add plot of stds, and add all other samples to plot
#possibly to be interactive so that can hover over samples to show ids? (guess not for pdf)
#add readme to github
#add fieldncs to tests
#add instruction for Sampling_Point for non-field controls


##############################################

#step1 - formatting cols
DF$Sample_Type<-as.character(DF$Sample_Type)
DF$Sampling_Point<-as.character(DF$Sampling_Point)
DF$Extraction_Batch<-as.character(DF$Extraction_Batch)
DF$DNA_Sample<-as.character(DF$DNA_Sample)
DF$Replicate<-as.character(DF$Replicate)
DF$Target_Cq<-as.numeric(DF$Target_Cq)
DF$IPC_Cq<-as.numeric(DF$IPC_Cq)
DF$Std_Conc<-as.numeric(DF$Std_Conc)
#assume all NAs in Target_Cq and IPC_Cq are 0
DF$Target_Cq[is.na(DF$Target_Cq)]<-0
DF$IPC_Cq[is.na(DF$IPC_Cq)]<-0

#remove NA rows
DF<-DF[DF$DNA_Sample!="",]

#some columns must be filled
if(nrow(DF[DF$Sampling_Point=="",])!=nrow(DF)) stop("All replicates must have a Sampling_Point") #must all have something
if(nrow(DF[DF$Replicate=="",])!=nrow(DF)) stop("All replicates must be assigned Replicate id") #must all have something
#DNA_sample*replicate must be unique
ids<-paste0(DF$DNA_Sample,DF$Replicate)
if(length(unique(ids))!=length(ids)) stop("Duplicate DNA_Sample:Replicate found")

#only allow "std","pcrnc","extnc","unkn", "fieldnc","pc" in DF$Sample_Type
if(sum(
  DF$Sample_Type == "std" | DF$Sample_Type =="pcrnc" | DF$Sample_Type == "extnc" |
  DF$Sample_Type == "unkn" | DF$Sample_Type == "fieldnc" | DF$Sample_Type == "pc"
)!=nrow(DF)
) stop("Sample_Type must be one of: std, pcrnc, extnc, unkn, fieldnc, pc")

#split DF by sample type for later
DFlist<-split(DF,DF$Sample_Type)

#small settings function to get setting from Settings, based on shortname
getset<-function(Settingsdf,setting_shortname){
  Settingsdf[match(setting_shortname,table = Settingsdf$shortname),"Setting"]
}

##############################################
##Put assay on validation scale (low,medium,high)
if(sum(validation[validation$Level_of_Confidence=="Low","user_input"])<5) {
  val_scale<-"Not usable"
} else if(sum(validation[validation$Level_of_Confidence=="Medium","user_input"])<3) {
  val_scale<-"Low"
} else if(sum(validation[validation$Level_of_Confidence=="High","user_input"])<5) {
  val_scale<-"Medium"
} else val_scale<-"High"

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
DF$interp1<-"NOT_DONE"

DF$interp1[which(DF$LOD_threshold=="PASSED" & DF$IPC_threshold=="PASSED")] <-"positive"
DF$interp1[which(DF$LOD_threshold=="PASSED" & DF$IPC_threshold=="IPC NOT USED")] <-"positive"
DF$interp1[which(DF$LOD_threshold=="PASSED" & DF$IPC_threshold=="FAILED")] <-"positive"

DF$interp1[which(DF$LOD_threshold=="NO DETECTABLE CURVE" & DF$IPC_threshold=="PASSED")] <-"negative"
DF$interp1[which(DF$LOD_threshold=="NO DETECTABLE CURVE" & 
                   DF$IPC_threshold=="IPC NOT USED")] <-"not detected, inhibition not tested"

DF$interp1[which(DF$LOD_threshold=="NO DETECTABLE CURVE" & 
                   DF$IPC_threshold!="PASSED" & 
                   DF$IPC_threshold!="IPC NOT USED")] <-"not detected, PCR inhibited"

DF$interp1[which(DF$LOD_threshold!="NO DETECTABLE CURVE" & 
                   DF$LOD_threshold!="PASSED" &
                   DF$IPC_threshold=="PASSED")] <-"weak signal below LOD, no inhibition observed"

DF$interp1[which(DF$LOD_threshold!="NO DETECTABLE CURVE" & 
                   DF$LOD_threshold!="PASSED" &
                   DF$IPC_threshold!="PASSED" &
                   DF$IPC_threshold=="IPC NOT USED")] <-"weak signal below LOD, inhibition not tested"

DF$interp1[which(DF$LOD_threshold!="NO DETECTABLE CURVE" & 
                   DF$LOD_threshold!="PASSED" &
                   DF$IPC_threshold!="PASSED" &
                   DF$IPC_threshold!="IPC NOT USED")] <-"weak signal below LOD, inhibition observed"


#INTERP2

#check on all the controls
#extnc
DFextnc<-DF[DF$Sample_Type=="extnc" & DF$interp1=="positive",]
batches_cont<-unique(DFextnc$Extraction_Batch)
#pcrnc
DFpcrnc<-DF[DF$Sample_Type=="pcrnc" & DF$interp1=="positive",]
if(nrow(DFpcrnc)>0) pcrnc_cont=1 else pcrnc_cont=0
#fieldnc
DFfieldnc<-DF[DF$Sample_Type=="fieldnc" & DF$interp1=="positive",]
Sampling_Point_batches_cont<-unique(DFfieldnc$Sampling_Point)
#pc
DFpc<-DF[which(DF$Sample_Type=="pc" & DF$interp1!="positive"),]
if(nrow(DFpc)>0) pc_amp=1 else pc_amp=0

#add control outcomes to interp2
#extnc
DF$extncs<-"NOT_DONE"
DF$extncs[which(is.na(DF$Extraction_Batch))]<-"NO BATCH SPECIFIED"
DF$extncs[which(DF$Extraction_Batch==batches_cont)]<-"CONTAMINATED"
DF$extncs[which(DF$Extraction_Batch!=batches_cont & DF$extncs!="NO BATCH SPECIFIED")]<-"CLEAN"
#fieldnc
DF$fieldncs<-"NOT_DONE"
DF$fieldncs[which(DF$Sampling_Point==Sampling_Point_batches_cont)]<-"CONTAMINATED"
DF$fieldncs[which(DF$Sampling_Point!=Sampling_Point_batches_cont)]<-"CLEAN"
#pcrnc
DF$pcrncs<-"NOT_DONE"
if(pcrnc_cont==1) DF$pcrncs="CONTAMINATED" else DF$pcrncs="CLEAN"
#pcs
DF$pcs<-"NOT_DONE"
if(pc_amp==1) DF$pcs="NOT ALL PCS AMPLIFIED" else DF$pcs="ALL PCS AMPLIFIED"

#INTERP 2 COL
DF$interp2<-"NOT_DONE"
DF$interp2[which((DF$interp1=="positive" | grepl("weak*",DF$interp1)) 
                 & DF$pcrncs=="CONTAMINATED")]<-"Inconclusive: PCR negative contaminated"

DF$interp2[which((DF$interp1=="positive" | grepl("weak*",DF$interp1)) 
                 & DF$pcrncs!="CONTAMINATED" 
                 & DF$extncs=="CONTAMINATED")]<-"Inconclusive: extraction batch contaminated"

DF$interp2[which((DF$interp1=="negative" | grepl("not detected*",DF$interp1)) 
                 & DF$pcs=="NOT ALL PCS AMPLIFIED")]<-"Inconclusive: not all positive controls amplified"

DF$interp2[which(DF$interp2=="NOT_DONE")]<-DF$interp1[which(DF$interp2=="NOT_DONE")]

DF$Interpretation<-DF$interp2

#final table 1 to report
DF1<-DF[,colNames]

##############################################
#DNA sample level report

#create dummy samples to ensure all categories are present
dummydf<-data.frame(DNA_Sample=paste("dummygg43266tahte6ytau",1:11),
                    Interpretation=c("negative",
                               "Inconclusive: PCR negative contaminated",
                               "Inconclusive: negative, PCR inhibited",
                               "not detected, PCR inhibited",
                               "not detected, inhibition not tested",
                               "weak signal below LOD, inhibition not tested",
                               "weak signal below LOD, inhibition observed",
                               "weak signal below LOD, no inhibition observed",
                               "positive",
                               "Inconclusive: extraction batch contaminated",
                               "Inconclusive: not all positive controls amplified"))

DFtemp<-rbind(DF[,c("DNA_Sample","Interpretation")],dummydf)

DF2<-as.data.frame.matrix(table(DFtemp$DNA_Sample,DFtemp$Interpretation))
DF2$DNA_Sample<-row.names(DF2)
DF2<-DF2[-grep("dummygg43266tahte6ytau",DF2$DNA_Sample),]

#add back sample type
DF2$Sample_Type<-DF[match(DF2$DNA_Sample,DF$DNA_Sample),"Sample_Type"]
#add n reps
DF2$nreps<-rowSums(DF2[,1:11])

#interp positives
DF2$interp_pos<-"NOT_DONE"

for(i in 1:nrow(DF2)){
  if(DF2$`Inconclusive: PCR negative contaminated`[i]>0) DF2$interp_pos[i]<-"Inconclusive: PCR negative contaminated"
  
  if(DF2$`Inconclusive: PCR negative contaminated`[i]==0
     & DF2$`Inconclusive: extraction batch contaminated`[i]>0) DF2$interp_pos[i]<-"Inconclusive: extraction batch contaminated"
  
  if(DF2$`Inconclusive: PCR negative contaminated`[i]==0
     & DF2$`Inconclusive: extraction batch contaminated`[i]==0
     & DF2$positive[i]>=getset(Settings,"Pos_minreps")) {
    DF2$interp_pos[i]<-
    paste0("Conclusive: ",DF2$positive[i], " of ",DF2$nreps[i],
          if(DF2$nreps[i]==1) " replicate" else " replicates",
          " positive and passing all controls (min=", getset(Settings,"Pos_minreps"),")")
  }
  
  if(DF2$`Inconclusive: PCR negative contaminated`[i]==0
     & DF2$`Inconclusive: extraction batch contaminated`[i]==0
     & DF2$positive[i]<getset(Settings,"Pos_minreps")
     & DF2$positive[i]>0) {
    DF2$interp_pos[i]<-
      paste0("Inconclusive: less than ", getset(Settings,"Pos_minreps"),
             " replicates positive and passing all controls, suggest further analysis")
  }
}

####
for(i in 1:nrow(DF2)){
  
  if(DF2$interp_pos[i]=="NOT_DONE" 
                       & DF2$negative[i]==DF2$nreps[i]) {
    DF2$interp_pos[i]<-"Conclusive: all replicates negative without inhibition"
  }
}

####
for(i in 1:nrow(DF2)){
  
  sumweak<-sum(DF2[i,grep("weak*",colnames(DF2))])
  
  if(DF2$interp_pos[i]=="NOT_DONE" & sumweak>0) {
    DF2$interp_pos[i]<-paste0("Inconclusive, weak signal in ", sumweak, 
                              if(sumweak==1) " replicate" else " replicates")
  }
    
    if(DF2$interp_pos[i]=="NOT_DONE" & sumweak==0
       & DF2$`not detected, PCR inhibited`[i]==0
       & DF2$nreps[i]-DF2$`not detected, inhibition not tested`[i]<=getset(Settings,"IPC_minreps")) {
      DF2$interp_pos[i]<-paste0("Inconclusive, not detected but inhibition was not tested in ", 
                                DF2$`not detected, inhibition not tested`[i],
                                if(DF2$`not detected, inhibition not tested`[i]==1) " replicate" else " replicates",
      " (min=",getset(Settings,"IPC_minreps"),")")
    }
      
      if(DF2$interp_pos[i]=="NOT_DONE" & sumweak==0
         & DF2$`not detected, PCR inhibited`[i]==0
         & DF2$nreps[i]-DF2$`not detected, inhibition not tested`[i]>getset(Settings,"IPC_minreps")) {
           DF2$interp_pos[i]<-paste0("Conclusive: Note that ", 
                                     DF2$negative[i]," replicates were negative, ",
                                     DF2$`not detected, PCR inhibited`[i]," showed inhibition and ",
                                     DF2$`not detected, inhibition not tested`[i],
                                     if(DF2$`not detected, inhibition not tested`[i]==1) " was" else " were",
                                     " not tested for inhibition")
      }
}

####
for(i in 1:nrow(DF2)){
  
  if(DF2$interp_pos[i]=="NOT_DONE" & DF2$`not detected, PCR inhibited`[i]<=getset(Settings,"IPC_Inh_maxreps")){
    DF2$interp_pos[i]<-paste0("Conclusive Negative: Note that ",DF2$negative[i], 
                              if(DF2$negative[i]==1) " replicate" else " replicates",
                              " were negative, ", DF2$`not detected, PCR inhibited`[i]," showed inhibition (min=",
                              getset(Settings,"IPC_minreps"),") and ",DF2$`not detected, inhibition not tested`[i],
                              if(DF2$`not detected, inhibition not tested`[i]==1) " was" else " were",
                              " not tested for inhibition")
  }
  
  if(DF2$interp_pos[i]=="NOT_DONE" & DF2$`not detected, PCR inhibited`[i]<=getset(Settings,"IPC_Inh_maxreps")){
    DF2$interp_pos[i]<-paste0("Conclusive Negative: Note that ",DF2$negative[i], 
                              if(DF2$negative[i]==1) " replicate" else " replicates",
                              " were negative, ", DF2$`not detected, PCR inhibited`[i]," showed inhibition (max=",
                              getset(Settings,"IPC_Inh_maxreps"),") and ",DF2$`not detected, inhibition not tested`[i],
                              if(DF2$`not detected, inhibition not tested`[i]==1) " was" else " were",
                              " not tested for inhibition (min=",getset(Settings,"IPC_minreps"),")")
  }
}

DF2$Interpretation<-DF2$interp_pos

#final table 2 to report
DF2<-DF2[,c("DNA_Sample","Sample_Type","nreps","Interpretation")]

##############################################
#Sampling point level report

DF3<-DF2
#add back sampling points
DF3$Sampling_Point<-DF1[match(DF3$DNA_Sample,DF1$DNA_Sample),"Sampling_Point"]
DF3<-as.data.frame.matrix(table(DF3$Sampling_Point,DF3$Interpretation))

#dev notes.
#validation 
#low: false pos possible and false neg possible
  #cannot say anything really
#medium: no false pos, possible false neg
  #can talk about 
#high: no false pos, no false neg


##############################################
#This DF definition is important - this is what shiny will take for the pdf output, along with the colNames vector for column headers
#Also important that the name "DF" is not changed
DF<-cbind(DF,data.frame("Interpretation"=Int))
##############################################


