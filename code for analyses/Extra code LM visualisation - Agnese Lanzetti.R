
#Code from Agnese Lanzetti - visualising LMs on the mesh 

#Pairs of landmarks, real and symmed
landpairs_DA <- cbind(symm = fixed_LMs, DA = fixed_LMs+64)

#Create matrix for linear distances
distance_matrix_DA <-matrix(nrow = length(Ids), ncol = length(fixed_LMs))

for (r in 1:length(landpairs_DA[,1])){
  distance_matrix_DA[,r] <- linear.dist(coords_DA_symm, landpairs_DA[r,])
  
}

distance_matrix_DA <- t(distance_matrix_DA)

#Create array to be filled with distances
distance_symm_DA <- array(dim = c(nrow(distance_matrix_DA), 1, ncol(distance_matrix_DA))) #1 is the columns of data we need  - see output of previous line

#Loop linear distances calculations
for (w in 1:ncol(distance_matrix_DA)){
  distance_symm_DA[,,w] <- distance_matrix_DA[,w]
  
}
#Set dimnames as specimens
dimnames(distance_symm_DA)[[3]] <- dimnames(coords_DA_symm)[[3]]

#Calculate mean all data
distance_symm_DA_R_mean <- apply(distance_symm_DA_R, c(1), mean)

#Colors shpres 3D plot
colfunc_flip <- colorRampPalette(c("lightyellow1", "yellow", "orange" , "darkorange", "red", "darkred"))

#Plot landmarks colored by distance on mean shape
#4. Create colour gradient of distance values based on colour palette created earlier
DA_col <- colfunc_flip(30)[as.numeric(cut(distance_symm_DA_mean, breaks = 30))]

#5. Check output by plotting points as per-point DA distance
plot(distance_symm_DA_mean, col = DA_col, cex = 2, pch = 19)

#6. Create 3D plot of per-landmark DA distance, based on single reference specimen
spheres3d(coords[fixed_LMs,,23], col = DA_col, type = "s",
          radius = 0.001, aspect = T, main = "mean",axes = F, main = F, fov = 0)

#DA symm distances
shade3d(template_mesh, col = "white", alpha = 0.5, fastTransparency = T)
spheres3d(shape_array[fixed_LMs,,23], col =  DA_col, type = "s",
          radius = 1.5, aspect = T, main = "mean",axes = F, main = F, fov = 0)
