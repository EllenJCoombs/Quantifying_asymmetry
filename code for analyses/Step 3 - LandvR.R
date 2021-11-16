

#==========

if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
Mirrored_skull_AC[1,,2] == Full_skull_AB[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(Mirrored_skull_AC[1,,2], digits = 9) == round(FUll_skull_AB[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = Mirrored_skull_AC[,,65],
                                                  reference = Full_skull_AB[,,65],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
differences_between_lms[[1]][1:66, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(Full_skull_AB[,,65], Mirrored_skull_AC[,,65], col.val = differences_between_lms[[1]][,1], col = colfunc)


test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######
N=123 #number of landmarks 
specs=162 #number of specimens 
all_combined=array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = Mirrored_skull_AC[,,i],
                                            reference = Full_skull_AB[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#55, 56, 57, 59, 60 are all missing data and should be zero 
all_combined[1:66, 1:3, 1:162] <- c(0.000000, 0.000000, 0.000000)
#write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #looking at the second column (usually x,y,z) but here it is the radii, aziumuth, and polar 
radii_mean=apply(radii, c(1), mean) #c(1) looking at the first column which is the radii 
#test=all_combined[[1]][,,1] #this is a test for specimen 1 to look at the differences between lms 
#test

radii=all_combined[,1,] #second column of whole dataset with just the radii [,1,]


#Looking at the average radii compared to specimen 21 (or an average specimen)
get.col.spectrum <- landvR::procrustes.var.plot(Full_skull_AB[,,9], Mirrored_skull_AC[,,9], col.val = radii_mean, col = colfunc)
#datcol2<-c(rep("black",66),get.col.spectrum)
#open3d()
