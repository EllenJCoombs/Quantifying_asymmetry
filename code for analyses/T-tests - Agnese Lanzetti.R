
#T-test to see whether there is significant variation between the manually landmarked side and the computer mirrored side of the skull

#Matrix with just the radius column from landVR output
radii <- symm_diff_radii[,1,]

#Create array to be filled with radii only
radii_mirror_real <- array(dim = c(nrow(radii), 1, ncol(radii))) #1 is the columns of data we need  - see output of previous line

#Loop linear distances calculations
for (b in 1:ncol(radii)){
  radii_mirror_real[,,b] <- radii[,b]
}

#Set dimnames as specimens
dimnames(radii_mirror_real)[[3]] <- dimnames(coords_LMs)[[3]]

#Overall variance in the dataset - max specimen - min specimen
var_range_radii <- rowSums(apply(radii_mirror_real, c(1,2), var))

#Pull out the fixed landmarks from the LHS (same in manual and mirrored data set)
landmarks_fix <- Mirrored_skull_AC[c(1:18),,]

#Test difference between sides in entire dataset - is the variance on the right higher than on the left+midline ?
side_test <- t.test(var_range_radii[landmarks_fix],
  var_range_radii[-landmarks_fix])

#Summary with p-value
side_test
