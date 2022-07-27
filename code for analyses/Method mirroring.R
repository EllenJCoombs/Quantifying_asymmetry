
#Code for mirroring symmetrical landmarks only 
#Asymmetrical landmarks are already manually placed on the skull 

midline <-as.integer(c(38,40,48,49,51,54,55,56,61))

asymmLMs <- Shape_data_with_bilats #see above
#define the curves and landmarks on each side 
left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
right.lm <- c(67:81, 85:91,95, 104, 109, 119:123) 

leftside<-c(left.lm) 
rightside<-c(right.lm)
num.missing<-(length(leftside)-length(rightside)) 
blanks<-c((dim(manual_skull_LMs)[1]+1):(dim(manual_skull_LMs)[1]+num.missing)) #to fill in blanks from one row past the last current point, for the number of rows needed (num.missing)
rightside<-c(rightside,blanks)
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]),
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)),
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}
specimens<-add_col_or_row(manual_skull_LMs,n=num.missing,add_col=FALSE,fill=NA)
dimnames(specimens)[3]<-dimnames(manual_skull_LMs)[3] #make sure the specimens match up 
bilats<-cbind(leftside,rightside) #bind the left and the right side 
newarray<-mirrorfill(specimens,l1=midline,l2=bilats) #newarray = final, correctly mirrored landmarks 
#dimnames(newarray)[3]<-dimnames(asymmLMs)[3] #make sure the specimens match up

newarray2 <- newarray[-c(92:123),,] #delete the duplicated mirror LMs 
                                    #97:123 on the non-conservative skull because several rostrum landmarks were manually placed 

#plot how the mirrored landmarks look
open3d();
spheres3d(newarray[,,1],radius=1.5) #plot whole skull to check asymmetric curve and symmetric curve placement 
spheres3d(newarray[bilats[,1],,1],col='red',radius = 1.5) #plot left side
spheres3d(newarray[bilats[,2],,1],col='blue',radius = 5) #plot right side 
spheres3d(newarray[midline,,1], col = 'yellow', radius = 1.5) #plot midline

#Check the specimens
atarfa=ply2mesh(file="E:/Ply ASCII/ply ASCII/Odonts/ply/Delphinapterus leucas USNM 305071.ply")
shade3d(atarfa,col='white')
spheres3d(method_skull_LMs[c(67:81, 85:91,95, 104, 109),,21], col = 'blue', radius = 4) #add
spheres3d(method_skull_LMs[c(1:66),,21], col = 'red', radius = 4)
spheres3d(method_skull_LMs[c(82, 83, 84, 92:94, 96, 97, 98:103, 105:108, 110:123),, 21], col = 'green', radius = 4)
text3d(newarray2[,,21], text = 1:123)


atarfa=ply2mesh(file="E:/Ply ASCII/ply ASCII/Odonts/ply/Delphinapterus leucas USNM 305071.ply")
shade3d(atarfa,col='white')
#spheres3d(final_manual_odonts[,,21], col = 'red', radius = 4)
#spheres3d(final_dataset[c(1:123),,21], col = 'blue', radius = 4)
spheres3d(final_dataset[c(1:123),,178], col = 'red', radius = 4)
spheres3d(final_dataset[c(1:123),,335], col = 'green', radius = 4)



plotRefToTarget(manual_skull_LMs[,,21], newarray2[,,21], method = "vector", radius = 1)
texts3d(newarray2[,,21], text = 1:123)
