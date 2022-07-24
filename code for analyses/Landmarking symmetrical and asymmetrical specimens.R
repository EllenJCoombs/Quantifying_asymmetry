
library(rgl)
library(paleomorph)
library(SURGE)
library(Morpho)

#This step comes after Step 4 - Informing the curve protocol – which curves to manually place - see paper 
#You will now have specimens that are either landmarked and semi-landmarked on one side (i.e., bialterally symmmetrical specimens) 
#Or specimens that have had had more manual placement of semi-landmarks because they are asymmetrical (see paper for Step 3: Locating asymmetry in the skull or structure)

##### IMPORTANT NOTE ########
#There are several general steps (i.e., not related to quantifying asymmetry specifically) that should be taken after Step 4: Informing the curve protocol – which curves to manually place, before mirroring curves. 
#These additional steps have been highlighted as a side note and are available in the paper 
#See 'Additional steps before running geometric morphometric analyses' both in the 'code for analyses' and the Supporting information. 

##########################################
#                                        #
#   landmarking symmetrical specimens    #
#                                        #
##########################################

#import .csv defining curves 
curve_table <- read_csv('new curves.csv')
my_curves <- create_curve_info(curve_table, n_fixed = 123) #define fixed curves 

#Symm_LMs is the data set of bilaterally symmetrical specimens - LMs and semi-LMs can be mirrored for these specimens 

#These data are slid with the fake landmarks 
symm_LMs <- Shape_data_with_bilats
symm_LMs[c(67:123),,]<-NA

open3d();spheres3d(symm_LMs[,,1])
left.curves<-c(1:64)
left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
right.lm <- c(67:123)
left.curve.list<-unlist(my_curves$Curve.in[left.curves])
leftside<-c(left.lm,left.curve.list) # RHS LMs+ RHS curves+all patch points
num.missing<-(length(leftside)-length(right.lm)) # number of LMs to create= total RHS-current LHS LMs
blanks<-c((dim(symm_LMs)[1]+1):(dim(symm_LMs)[1]+num.missing))
# to fill in blanks from one row past the last current point, for the number of rows needed (num.missing)
rightside<-c(right.lm,blanks)
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]),
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)),
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}
specimens2 <- add_col_or_row(symm_LMs, n=num.missing, add_col=FALSE,fill=NA)
dimnames(specimens2)[3]<-dimnames(symm_LMs)[3]
bilats<-cbind(leftside,rightside)
newarray<-mirrorfill(specimens2,l1=midline,l2=bilats)
dimnames(newarray)[3]<-dimnames(symm_LMs)[3]

symmLMs_final <- newarray
              
##########################################
#                                        #
#   landmarking asymmetrical specimens   #
#                                        #
##########################################

#import .csv defining curves 
curve_table <- read_csv('new curves.csv')
my_curves <- create_curve_info(curve_table, n_fixed = 123) #define fixed curves 

#asymmLMs is resampled, slid landmarks with missing landmarks and variably present bones corrected (see github for full code and ‘side notes’ section of methods)

#midline landmarks (anchor points via which landmarks are mirrored) in this example:
midline <-as.integer(c(38,40,48,49,51,54,55,56,61,1114,1115))

asymmLMs <- Shape_data_with_bilats #see above
#define the curves and landmarks on each side 
left.curves<-c(1:64) 
left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66) 
right.lm <- c(120,121,67:82,122,123,83:119)
right.curves <- c(65:85)
left.curve.list<-unlist(my_curves$Curve.in[left.curves])
right.curve.list<-unlist(my_curves$Curve.in[right.curves])
leftside<-c(left.lm,left.curve.list) 
rightside<-c(right.lm, right.curve.list)
num.missing<-(length(leftside)-length(rightside)) 
blanks<-c((dim(asymmLMs)[1]+1):(dim(asymmLMs)[1]+num.missing)) #to fill in blanks from one row past the last current point, for the number of rows needed (num.missing)
rightside<-c(rightside,blanks)
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]),
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)),
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}
specimens<-add_col_or_row(asymmLMs,n=num.missing,add_col=FALSE,fill=NA)
dimnames(specimens)[3]<-dimnames(asymmLMs)[3] #make sure the specimens match up 
bilats<-cbind(leftside,rightside) #bind the left and the right side 
newarray<-mirrorfill(specimens,l1=midline,l2=bilats) #newarray = final, correctly mirrored landmarks 
dimnames(newarray)[3]<-dimnames(asymmLMs)[3] #make sure the specimens match up

#plot how the mirrored landmarks look
open3d();
spheres3d(newarray[,,3],radius=1.5) #plot whole skull to check asymmetric curve and symmetric curve placement 
spheres3d(newarray[bilats[,1],,1],col='red',radius = 1.5) #plot left side
spheres3d(newarray[bilats[,2],,1],col='blue',radius = 1.5) #plot right side 
spheres3d(newarray[midline,,1], col = 'yellow', radius = 1.5) #plot midline

 asymmLMs_final <- newarray
             
########################################################
#                                                      #
#   Now combine the symmetrical specimens +            #
#   asymmetrical specimens into one universal data     #
#   set to start your geometric morphometric           #
#   analyses                                           #
#                                                      #
########################################################
                        
final_dataset_undescribed=abind::abind(symmLMs_final, asymmLMs_final along = 3) #column 3 is specimens - this binds the coords data in alphabetical order         
              
#Now run your usual GM analyses on a full data set that has all specimens (both symmetrical and asymmetrical) fully LMed and curved in a universal manner           
              
