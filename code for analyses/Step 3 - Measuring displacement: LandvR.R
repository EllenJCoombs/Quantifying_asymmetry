
#Qauntifying landmark displacement 

####################################
#                                  #
#   Step 3 - Coombs and Felice     #
#                                  #
####################################

if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
Mirrored_skull_AC[1,,2] == Full_skull_AB[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(Mirrored_skull_AC[1,,2], digits = 9) == round(Manual_skull_AB[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = Mirrored_skull_AC[,,65],
                                                  reference = Manual_skull_AB[,,65],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
differences_between_lms[[1]][1:66, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(Manual_skull_AB[,,65], Mirrored_skull_AC[,,65], col.val = differences_between_lms[[1]][,1], col = colfunc)


test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######

library(landvR)
N=123 #number of total landmarks 
nfixed = 66 #number of fixed landmarks 
specs = 157 #number of specimens 
k = 3 #number of dimensions in the matrix
colfunc <- colorRampPalette(c("red", "yellow", "white")) #create colour function for visualising landmarks
colfunc(10) #choose number of increments for colour scale 

#make the array for analyses: 
all_combined = array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)
#manual_data is the fully landmarked skull (reference data)
#mirrored_data is the half-landmarked skull that has been mirrored 
#calculate the coordinates.differences between these data sets (i.e. how much the landmarks move between the manually placed landmarks and the mirrored landmarks) 

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = Mirrored_skull_AC[,,i],
                                            reference = Manual_skull_AB[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#landmarks 1:66 (nfixed) in this example are fixed and therefore have the value of zero 
all_combined[1:nfixed, 1:k, 1:specs] <- c(0.000000, 0.000000, 0.000000)

#save output if desired using: write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #radius per landmark for each specimen (second column of whole dataset with just the radii [,1,])

radii_mean=apply(radii, c(1), mean) #c(1) look at the first column - the radii 
#radii_mean is a mean radius value per landmark 

#save radii and radii_mean as .csv files for further analyses 

#example
#looking at the average radii compared to specimen 21 (or choose an average specimen if preferred)
get.col.spectrum <- landvR::procrustes.var.plot(Manual_skull_AB[,,21], Mirrored_skull_AC[,,21], col.val = radii_mean, col = colfunc)


