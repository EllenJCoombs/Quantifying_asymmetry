
#IMPORTANT NOTE - this code should be run after Step 4: Informing the curve protocol – which curves to manually place (see notes below)

#There are several general steps (i.e., not related to quantifying asymmetry specifically) that should be taken after 
#Step 4: Informing the curve protocol – which curves to manually place, before mirroring curves
#These additional steps have been highlighted here as a side note and are available in the detail 
#in the Supporting Information: Section 1: Additional steps before running geometric morphometric analyses 

#The example shown here is run on the asymmetrical data set (.pts) - run all the same steps for symmetrical data set (.pts)
#CAUTION - you will need a different curve_table for asymm and symm specimens as asymm specimens will have more manually placed curves (see below) 

library(tidyverse)
library(broom) #needed with tidyverse
library(latticeExtra) #qgraph
library(readr)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph) #load in packages bit 
library(ape)
library(geiger)
library(abind)
library(tibble)
library(dplyr)
library(tidyr)
library(rgl)
library("devtools")
devtools::install_github("rnfelice/SURGE")

install.packages("remotes")

library(SURGE)
library(RColorBrewer) # for color palettes
library(magick)

#import table defining curves for the asymmetrical specimens 
#different curves table for the asymmetrical specimens 
curve_table <- read_csv('new curves_asymm.csv')

#different curve_table for symmetrical specimens (only 64 curves)
curve_table <- read_csv('new curves_symm.csv')

#identify the folder where your pts files are
ptsfolder <- "D:/PTS WHOLE SKULL/asymm specs/pts"


#import the pts file names
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)

my_curves <- create_curve_info(curve_table, n_fixed = 123)
setwd(ptsfolder)
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE)

#capture the output if it's too big to read in the console 
capture.output(import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE), file = "odonts.txt", append = T)
#if you have any missing points, Checkpoint will set them to x,y,z=9999
#this makes for bad plots in checkLM, so switch them to NA
subsampled.lm[subsampled.lm == 9999] <- NA

#TEST TO SEE HOW THEY LOOK - check for any missing or incorrect curves 
#SET WORKING DIRECTORY TO ASCII PLY
#check to make sure your curves look okay on each specimen
checkLM(subsampled.lm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)

newpts <- subsampled.lm

#Create missing list 
misslist<-createMissingList(dim(newpts)[3])
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 
newpts2<-fixLMtps(newpts)
#check that your pts and ply match
plyfolder <- "D:/Ply ASCII/ply ASCII/Longirostrine/Schizodelphis/ply"
ptslist2<-gsub(pattern="\\.pts$","",ptslist)
plylist <-  dir(plyfolder, pattern='.ply', recursive=F)
plylist2<-gsub(pattern="\\.ply$","",plylist)
setdiff(plylist2,ptslist2) #should be zero if all matching up OK
#######NOTE
#####my plys are awkwardly in a folder called  "C:/Data/whale_ply/ply/ply"
##### I am setting the working directory ONE FOLDER UP from where the meshes are

#mc.cores = 1 for WINDOWS. mc.cores = 3 for MAC
####so for you it should be 

#################################################
#                                               #
#        RUN THE SLIDER 3D_2 CODE               #
#                                               #
#################################################

##################################
#                                #
#    STEP SIZE - 0.2 or 0.002    #
#                                #
##################################


setwd("D:/Ply ASCII/ply ASCII/Archs")
{
  slided4.all <- slider3d_2(newpts2$out, SMvector= my_curves$Sliding.LMs,
                          outlines = my_curves$Curves, sur.path = "./ply", sur.name = NULL, 
                          meshlist = paste("./ply/",dimnames(newpts2$out)[[3]],".ply",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = FALSE,stepsize=0.02,
                          missingList=misslist)
  dimnames(slided4.all[["dataslide"]])[3]<-dimnames(newpts2$out)[3]
  #save(slided4.all,file="~/Google Drive/NHM/crocs/data/slid.crocs.all.apr25.R")
  #load("~/Google Drive/NHM/crocs/data/slid.crocs.sept23.R")
}


#re-estimate missing post sliding
slided4.all$dataslide[which(is.na(newpts))]<-NA
#Fix missing landmarks
slid.lms<-fixLMtps(slided4.all$dataslide)
#the landmarks ready for mirroring are in an object called slid.lms$out


slidedlms <- slid.lms$out

checkLM(slidedlms,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)


slidedlms <- slidedlms
save(slidedlms, file = 'slidedlmsSCHIZO.R')


#####################
#                   #
#   ABSENT BONES    # 
#                   #
#####################

########################################################
#                                                      #
#  Dealing with specimens with variably present bones  #
#                                                      #
#                                                      #
########################################################



#add the data for absent bones for specific species 
#Jugal example 

absent<-read.csv("absent.csv")

lm_nasal_l <- my_curves$Curves[which(curve_table$Bone%in%c("nasal_l"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_nasal_r <- my_curves$Curves[which(curve_table$Bone%in%c("nasal_r"))]%>%unlist(.)%>%unique(.)%>%sort(.)

lm_max_v <- my_curves$Curves[which(curve_table$Bone%in%c("maxilla_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)

absentslid<-slid.lms

#testslid$out[nasal,,1] #specimen number on the end 
#testslid$out[lm_max_v,,207] #specimen number on the end 



#ventral maxilla - computer mirrored so only one needed 
for (i in 1:nrow(absent)){
  if( !is.na(absent$maxilla_v[i]))
    absentslid$out[lm_max_v,c(1:3),i] <- matrix(absentslid$out[61,c(1:3),i], nrow = length(lm_max_v), ncol=3,byrow=TRUE)
}



#LHS nasal - asymmetrical so pull out each side individually 
for (i in 1:nrow(absent)){
  if( !is.na(absent$nasal_l[i]))
    absentslid$out[lm_nasal_l,c(1:3),i] <- matrix(absentslid$out[1,c(1:3),i], nrow = length(lm_nasal_l), ncol=3,byrow=TRUE)
}



#LHS nasal 
for (i in 1:nrow(absent)){
  if( !is.na(absent$nasal_r[i]))
    absentslid$out[lm_nasal_r,c(1:3),i] <- matrix(absentslid$out[67,c(1:3),i], nrow = length(lm_nasal_r), ncol=3,byrow=TRUE)
}


slid.lms.neg<-absentslid #make this your new dataset 
slid.lms.neg$out[lm_nasal_l,,51] #test on the specimen with the bone missing 
slid.lms.neg$out[lm_max_v,,95] 

slidedlmsODONTS2 = slid.lms.neg$out 


#This step comes AFTER dealing with variably present bones, and BEFORE mirroring
##################################
#                                #
#                                #
#   ADD ADDITIONAL MIDLINE LMS   #
#   TO HELP ALIGNMENT            #
#                                #
##################################

LM1_bilat=slidedlmsZARHI[c(5,69),,] ## LM1 on left (5) and corresponding LM on the right side (69)
LM2_bilat=slidedlmsZARHI[c(6,70),,] ## LM2 on left (6) and corresponding LM on the right side (70)

#LM3_bilat=slidedlmsARCHS[c(58, 112),,] ##

#Check the dimensions 
#now want to find the midpoint of these pairs, 
#which is like finding the mean of each pair's position. 
#So to do this I summed each specimens x, then y then z coords, 
LM1_midline=colSums(LM1_bilat)/2
LM2_midline=colSums(LM2_bilat)/2
#LM3_midline=colSums(LM3_bilat)/2

#Then just visually check that you're happy the code worked (IMPORTANT!):
## read in a specimen ply and plot these landmarks on it:
Pipa=ply2mesh(file="D:/Ply ASCII/ply ASCII/Odonts/ply/Orcinus orca USNM 11980.ply") #This is a test specimen 
shade3d(Pipa, col="white")
spheres3d(LM1_midline[,1]) #bracketed is specimen number
spheres3d(LM2_midline[,1])
spheres3d(slidedlmsODONTS[c(5,69),, 1], col = 'green') #check how these look with 5, 69, 6, 70
spheres3d(slidedlmsODONTS[c(6,70),,1], col = 'red')
#spheres3d(slidedlmsARCHS[c(58,112),,1], col = 'red')
#Then you add these two new midline landmarks to your landmark set:
Shape_data_with_bilats=abind::abind(slidedlmsZARHI, LM1_midline, LM2_midline, along=1) 
## so here just add them to whatever dataset you were previously using 
#(that still includes the original 4 landmarks that you used for creating the midline points)

#Then include these two points in your midline definition
#symmetrical specimens 
midline<-as.integer(c(38,40,48,49,51,54,55,56,61,1114,1115)) 

#asymmetrical specimens
midline<-as.integer(c(38,40,48,49,51,54,55,56,61,1449,1450)) 
#check
#The additional LMs are higher up in the asymm specimens because they additional manually placed curves  
#Check symm specimens 
spheres3d(Shape_data_with_bilats[c(1114:1115),,1], radius = 1, col='red')
spheres3d(Shape_data_with_bilats[midline,,1], radius = 1, col = 'black')
spheres3d(Shape_data_with_bilats[c(1:123),,1], radius = 1, col = 'yellow')
spheres3d(Shape_data_with_bilats[c(124:1113),,1], radius = 1, col = 'green')

#Asymmetrical specimens 
midline<-as.integer(c(38,40,48,49,51,54,55,56,61,1449,1450)) 
#check
#The fake landmarks are 1449 and 1450 in the asymmetrical specimens 
spheres3d(Shape_data_with_bilats[c(1449:1450),,1], radius = 1, col='red')
spheres3d(Shape_data_with_bilats[midline,,1], radius = 1, col = 'black')
spheres3d(Shape_data_with_bilats[c(1:123),,1], radius = 1, col = 'yellow')
spheres3d(Shape_data_with_bilats[c(1114:1448),,1], radius = 1, col = 'green')
spheres3d(Shape_data_with_bilats[c(124:1113),,1], radius = 1, col = 'blue')

#MIRROR - the two data sets (asymm and symm) require different mirroring because of the difference in manually and computer mirrored LMs
#SEE 
