
library(tibble)
library(geomorph)
library(abind)
library(hablar)
library(dplyr)
library(ggplot2)


#load data 
load('manual_skull_LMs.R') #this is the manually LMd skull (123 LMs, 157 specimens)
load('mirror_skull_LMs.R') #this is the mirrored skull (57 LHS LMs, 9 midline LMs, mirrored to make 123 LMs total)

species_data <- read.csv('ODONTS_dataALL.csv') #this data set has info on whether the specimen was manually LMd or computer LMd (for plotting later)

#Note: all missing data has been dealt with in these data sets using code in 'Mirroing landmarks.R'
#Note: specimen names in the mirrored data set have _50 after them as only 50% of the skull was manually landmarked

#Bind the data set and Procrustes as one
final_dataset=abind::abind(manual_skull_LMs, mirror_skull_LMs.R, along = 3) 

#Add species names to data set (first 157 are manual LMd, second 157 are computer LMd (mirrored))
#This data set has the landmark names (computer or manual) for the plot 


#Procrustes the data 
final_data_PROC=geomorph::gpagen(final_dataset) #Remove non-shape aspects 
final_data_COORDS=final_data_PROC$coords #Subset out the coords 

#PCA
pca_res <- geomorph::gm.prcomp(final_data_COORDS)

#mutate the species data to add to the matrix
species_data <-  mutate(species_data, label = species_data$species)

pca_scores <- as_tibble(pca_res$x)
pca_scores <- pca_scores %>% mutate(., species=rownames(pca_res$x))

#join the coords and species data
pca_scores2 <- left_join(species_data, pca_scores, by = c("species"))
pca_scores2 <- pca_scores2 %>% convert(fct(landmarks))

#Plot - this is ugly compared to the final one but I will sort that later
g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=landmarks))+ 
  geom_point(size=2,aes(shape=landmarks, colour=landmarks))+
  scale_shape_manual(values=c(24, 22)) +
  scale_alpha(range=c(0.4, 1),guide=FALSE)+
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "57.2 % of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"8.9 % of Total Variance)")) +
  scale_colour_manual(values=c("#440154FF", "#2A788EFF")) +
  theme_bw()
  #geom_text(aes(label= number), hjust=1, vjust=2) # add this if you want it with numbers for each specimen

#plot
g

#Basic plot with specimen numbers on if needed
plot(pca_res, axis1 = 1, axis2 = 2)
#to look at specimen position 
text(pca_res$x[,1],pca_res$x[,2], pos = 4)




