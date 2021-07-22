## analysis of ACD staining : part 2
## july 21th, 2021

# CORRE Guillaume @ GENETHON, Evry, France

##
##
## Load data from ACD_{date}.rdata or create the object running 00-ACD.r
##
# Load libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(RImageJROI);
library(sp);
library(spatstat);




# Load the dataset --------------------------------------------------------

load("ACD_20210721.rdata")

# the environement contains 4 objects:
# - fibre : contains fiber information
# - spots : contains spots information and parent fiber id
# - nuclei : contains nuclei information and parent fiber id
# - polygon_rois : contains fiber boundaries as spatial objects, for ploting.







# filter spots ------------------------------------------------------------

      # remove spots that have a very low fluorescence compared to its parent fiber background

      ## add parent fiber information to each child
spots_parentFiber <- spots %>%
  filter(!is.na(ROI.fibre)) %>% 
  left_join(fibre %>% 
              select(sample,ROI.fibre,Mean_ACD.fibre), by = c("sample","ROI.fibre"))


ggplot(spots_parentFiber, aes((Mean_ACD.spots/Mean_ACD.fibre))) +
  geom_density() +
  scale_x_log10() +
  geom_vline(xintercept = 3,lty = 2)


spots_filtered = spots_parentFiber %>%
  filter((Mean_ACD.spots/Mean_ACD.fibre)>=5)
rm(spots_parentFiber)






# Count child object per fiber --------------------------------------------

spots_per_fibre <- fibre%>% 
  select(sample,ROI.fibre) %>%
  left_join(spots_filtered %>% 
              select(sample,ROI.fibre,ROI.spots), by = c("sample", "ROI.fibre")) %>%
  group_by(sample,ROI.fibre) %>%
  summarise(spot_count = length(which(!is.na(ROI.spots))))


nuclei_per_fibre <- fibre %>% 
  select(sample,ROI.fibre) %>%
  left_join(nuclei %>% 
              select(sample,ROI.fibre,ROI.nuclei), by = c("sample", "ROI.fibre")) %>%
  group_by(sample,ROI.fibre) %>%
  summarise(nuclei_count = length(which(!is.na(ROI.nuclei))))



fibre <- fibre %>% 
  left_join(spots_per_fibre, by = c("sample","ROI.fibre"))%>%
  left_join(nuclei_per_fibre, by = c("sample","ROI.fibre"))

rm(nuclei_per_fibre)
rm(spots_per_fibre)




# Filter fibers -----------------------------------------------------------

# we may want to remove some artefacts based on their shape , fluorescence intensity, position, children count.

ggplot(fibre, aes(Solidity.fibre)) + geom_density() 

ggplot(fibre, aes(AR.fibre)) + geom_density() 
ggplot(fibre, aes(log10(Mean_DAPI.fibre+1))) + geom_freqpoly() 
ggplot(fibre, aes(Mean_Dystrophin.fibre)) + geom_freqpoly()
ggplot(fibre, aes(MinFeret.fibre)) + geom_freqpoly()
ggplot(fibre, aes(nuclei_count)) + geom_bar()
ggplot(fibre, aes(spot_count)) + geom_bar()



fibres_filtered <- fibre %>% 
  filter(Solidity.fibre > .7, 
         AR.fibre < 4, 
         MinFeret.fibre < 500,
         log10(Mean_Dystrophin.fibre+1)>0.75,
         log10(Mean_DAPI.fibre+1)< 1.3,
         nuclei_count < 15,
         spot_count < 20)



write.table(x = fibres_filtered %>% select(sample,animal,slice,ROI.fibre,status,nuclei_count,spot_count), file = "fibers_childCount.tsv", sep = "\t", quote = F, row.names = F)



# Analyse child per fiber count -------------------------------------------

ggplot(fibres_filtered, aes(nuclei_count)) + geom_bar() + facet_wrap(~animal)

ggplot(fibres_filtered, aes(spot_count)) + geom_bar() + facet_wrap(~animal)


fibres_filtered %>% 
  count(sample,animal, slice, nuclei_count)
  
fibres_filtered %>% 
  count(sample,animal, slice, spot_count)



# Plotting functions -------------------------------------------------------


plot_slice <- function(slice_id, col = "status", fibre_df,spots_df,nuclei_df){
  
  fib = fibre_df %>% filter(sample == slice_id)
  fib_poly = fortify(polygon_rois[[slice_id]]) %>% inner_join(fib,by = c("id" = "ROI.fibre"))
  nuc = nuclei_df %>% filter(sample == slice_id)
  spo = spots_df %>% filter(sample == slice_id)
  

  ggp <- ggplot(fib_poly, aes(long,lat, group = id)) + 
      geom_polygon(lwd = .5, col = "white",fill="grey") + 
      coord_fixed() + 
      scale_y_reverse() + 
      theme_bw() +
      labs(fill = paste(col)) +
      geom_point(data = nuc, aes(X.nuclei,Y.nuclei), inherit.aes = F) +
      ggnewscale::new_scale_color()+
      geom_point(data = spo, aes(X.spots,Y.spots),col="red",inherit.aes = F)+
      scale_color_viridis_c() +
    ggtitle(label = slice_id)
print(ggp)
}


plot_slice(slice_id = "17-099_5e12-1e13-5e13_78-2", 
           fibre_df = fibres_filtered,
           spots_df = spots_filtered, 
           nuclei_df = nuclei,
           col = "Mean_ACD.fibre")




## save plots in pdf
pdf("plot_slice.pdf")



for(f in names(polygon_rois)){
  cat("...",f,"\n")
  plot_slice(slice_id = f, 
             fibre_df = fibres_filtered,
             spots_df = spots_filtered, 
             nuclei_df = nuclei)
}
dev.off()




