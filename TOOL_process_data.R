#if we're in the dashboard, use the information from there, otherwise upload locally for development purposes

if(!is.element("SHINY",ls())){
  #Define the settings matrix
  Settings<-data.frame("Description"=
                   c("Assay validation level",
                     "Minimum number of positive qPCR replicates for DNA sample to be considered Conclusive: positive (default: 2)",
                     "Cq value for Limit of Detection (default: determined from standard curve [0])",
                     "Inhibition Threshold, cycle delay relative to pcrncs (default: 1)",
                     "Override Inhibition Threshold, just use minimum cycle number that IPC should have attained (default:0 [OFF])",
                     "IPC threshold: minimum number of replicates that must have been tested for DNA sample to be considered conclusively negative (default:12)",
                     "IPC threshold2: maximum number of replicates that exhibited inhibition for DNA sample to be considered conclusively negative (default:0)"
                   ),
                 "Setting"=c(
                   2,
                   1,
                   60,
                   1,
                   0,
                   12,
                   4
                 ),
                 "shortname"=c(
                   "AVL",
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

#Override<-Settings$Setting[5] #Extract the "Override Inhibition Threshold.." value

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

#checks to add



#step1 - formatting cols
DF$sample_type
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

#split DF by sample type for later
DFlist<-split(DF,DF$sample_type)

#small settings function to get setting from Settings, based on shortname
getset<-function(Settingsdf,setting_shortname){
  Settingsdf[match(setting_shortname,table = Settingsdf$shortname),"Setting"]
}



#IPC filter
DF$IPC_threshold<-"NOT_DONE"

if(getset(Settings,"IPC_Cq")=0){
  DF$IPC_threshold[which(DF$IPC_Cq==0)]<-"FAILED" #fail IPC if IPC_Cq was 0
  DF$IPC_threshold[which(DF$IPC_Cq==-1)]<-"IPC NOT USED" #fail IPC if IPC_Cq was 0
  if(DF$IPC_threshold>-1){
    meanNegIPC<-mean(DFlist$pcrnc$IPC_Cq[which(DFlist$pcrnc$IPC_Cq>-1)]) #get average of IPC for pcr negs (where used)
  }
  
  DF$IPC_threshold[which(DF$IPC_Cq==0)]<-"FAILED"
  if(DF$IPC_Cq)
}
#IF(README!$B$12=0,     IF(G2=0,"FAILED", 
IF(G2>-1,IF(G2-AVERAGEIFS(G:G, C:C,"pcrnc",G:G,">-1")>README!$B$11,"FAILED","PASSED"),"IPC NOT USED")),     IF(G2=0,"FAILED", IF(G2>-1,IF(G2>README!$B$12,"FAILED","PASSED"),"IPC NOT USED")))

##############################################
#This DF definition is important - this is what shiny will take for the pdf output, along with the colNames vector for column headers
#Also important that the name "DF" is not changed
DF<-cbind(DF,data.frame("Interpretation"=Int))
##############################################


