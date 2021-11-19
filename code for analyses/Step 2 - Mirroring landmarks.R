

#Code adapted from Guillerme, T and Weisbecker V (2019). 
#landvR: Tools for measuring landmark position variation. Zenodo. doi:10.5281/zenodo.2620785

###########################################
#                                         #
#      Step 2 - Coombs and Felice         #
#                                         #
###########################################


#Looking at asymmetry in odontocetes and mysticetes
#Run for all specimens (157) 
#Run for a subset of bilaterally symmetrical specimens (mysticetes in this example)

##### PLEASE NOTE ########
#The code shows you how to read in your own .pts files to make them an array
#If you would like to replicate the results shown in this study rather than your own .pts, 
#please use the data from the data file: full skull = 'manual_skull_LMs.R' (LMs 1:123) and mirrored skull = 'mirrored_skull_LMs.R' (LMs 1:66)

#Code to run to pull out radii for assessing asymmetry in cetaceans

library(Rvcg)
library(rgl)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)


#========================================#
#      1. READ IN THE MANUAL LMS         #  #Or read in our data set 'manual_skull_LMs.R' if you are not using your own data and skip this part
#========================================#


###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa <-157 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(123,3,157)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Need the .plys for this 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))
arraylm<-ptsarray #this is your array


##### MISSING LANDMARKS #########
#do this first and then rearrange landmarks 
#Landmarks need to be 'reshuffled' because we added 4 extra LMs on the nasal to be reflect that one midline landmark would 
#not work in the asymmetrical odontocetes 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#careful not to have missing data 
#check here and also check all of the LMs are numbered correctly 
text3d(arraylm[,,2], text=1:dim(arraylm)[1])

#This would be to change the names according to the species (as above with .ply)
#dimnames(arraylm)[3]=species

#let's call arraylm 'manual skull' to differentiate from the mirrored skull LMs
manual_skull <- arraylm

#############################
#                           #
#   PROCRUSTES THE DATA     #  Or read in our data set 'mirror_skull_LMs.R' if you are not using your own data and skip this part 
#                           #
#############################

manual_skull=gpagen(manual_skull) #Remove non-shape aspects 
Manual_skull_AB=manual_skull$coords #Subset out the coords 
size=manual_skull$Csize #centroid size
#PlotTangentSpace if you want to see a quick morphospace of the skulls 

PCA <- plotTangentSpace(manual_skull$coords, legend=TRUE) 

#Manual_skull_AB is the Procrusted coords dataset of manually placed landmarks
#AB denotes manually LMed skull, AC denotes computer mirrored skull

############################################################
#                                                          #
#     2. Load second dataset - half LM skull to mirror     #  Or read in our data set 'mirror_skull_LMs.R' if you are not using your own data and skip to line 126
#                                                          #                          
############################################################

setwd("X:xxxxxxx/xxxxxx") 
#NB the plys have been transformed to ASCII (from binary) so that landmarks can be visualised on the mesh 
#See code 'Binary_ASCII_ply.R" if needed for this step

ntaxa <- 157 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarrayAC<-array(dim=c(66,3,157)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarrayAC[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Set to the whole project directory for this - the code wants to look for 'ply' folder 
dimnames(ptsarrayAC)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))
arraylmAC<-ptsarrayAC


##### CHANGING MISSING LANDMARKS BEFORE MIRRORING ######### - NB, use 'mirrored_skull_LMs.R' if you are not using your own data 
#Otherwise you get weird numbers that aren't 9999 mirroring
#estimate.missing from Geomorph is used for estimating missing landmarks 

arraylmAC[which(arraylmAC==9999)] <- NA
arraylmAC <- estimate.missing(arraylmAC,method="TPS")

#let's call arraylmAC ' skull' to differentiate from the mirrored skull LMs
manual_skull <- arraylm


#MIRROR THESE LANDMARKS over the central line plane of the skull 

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################


midline<-as.integer(c(38,40,48,49,51,54,55,56,61)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(arraylmAC)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,157)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(arraylmAC)[1],,]=arraylmAC

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(arraylmAC[left.lm,,10],radius=4) #left LMs
spheres3d(arraylmAC[midline,,10],radius=5,col='red') #midline

right.lm <- c(67:123) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(left.lm, right.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC
#deformGrid3d(MirroredAC[67:123,,2], Matrice[,,2], ngrid=0) #This shows you the new mirroed landmarks 

#These visualisations are done before Procrusted data 
#This shows the original landmarks
spheres3d(MirroredAC[c(1:66),,1],col=2,radius=4)
#check dimensions

#############################
#                           #
#   PROCRUSTES THE DATA     #
#                           #
#############################

MirroredAC=gpagen(MirroredAC) #Remove non-shape aspects 
Mirrored_skull_AC=Mirrored_data$coords #Subset out the coords 

#Checking if landmark ordering is correct 
#Comparing landmarks on both sides oft he skull 

col=rainbow(length(1:dim(Mirrored_data)[1]))
shapes3d(Y.gpa$coords[,,1], color=col)

#We now have two data sets with which to quantify asymmetry: the manually landmarked skull (Manual_skull_AB) and the computer mirrored data (Mirrored_skull_AC)
#Use these in Step 3


