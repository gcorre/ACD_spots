### analyse images marquage ACD

### Juillet 2021
### guillaume Corre @ genethon



# libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(RImageJROI);library(spatstat); library(sp)





# datasets ----------------------------------------------------------------

files = list.files(pattern = "fibre.csv")
files <- str_remove(files, pattern = "_fibre.csv")



csa <- list()


for(f in files){
  
  cat("processing....",f,"\n")
  
  cat("\t\tReading fibers\n")
  
  fibers <- fread(paste(f,"fibre.csv", sep ="_"))
  
  fibers <- fibers %>% select(-V1,-XM,-YM,-Slice)%>%
    separate(Label, into = c("image","ROI","channel"), sep = ":") %>%
    pivot_wider(names_from = "channel", values_from = c("Mean","Median","Min","Max","Mode","StdDev"))
  
  
  ####
  cat("\t\tReading fibers ROIs\n")
  fibers_roi <- read.ijzip(file = paste(f,"fibre.roi.zip", sep ="_"),list.files = F)
  
  p = lapply(fibers_roi, function(y){
    Polygon(y$coords)
  })
  
  p <- lapply(seq_along(p), function(i) Polygons(list(p[[i]]), ID = names(p)[i]))
  p <- SpatialPolygons(p)
  
  
  ####
  cat("\t\tReading nuclei\n")
  nuclei <- fread(paste(f,"nuclei.csv", sep ="_"))
  
  nuclei <- nuclei %>% select(-V1,-XM,-YM,-Slice)%>%
    separate(Label, into = c("image","ROI","channel"), sep = ":") %>%
    pivot_wider(names_from = "channel", values_from = c("Mean","Median","Min","Max","Mode","StdDev"))
  
  
  ####
  cat("\t\tReading nuclei ROIs\n")
  nuclei_roi <- read.ijzip(file = paste(f,"fibre.roi.zip", sep ="_"),list.files = F)
  
  p_nuclei = lapply(nuclei_roi, function(y){
    Polygon(y$coords)
  })
  
  p_nuclei <- lapply(seq_along(p_nuclei), function(i) Polygons(list(p_nuclei[[i]]), ID = names(p_nuclei)[i]))
  p_nuclei <- SpatialPolygons(p_nuclei)
  
  
  ####
  cat("\t\tReading spots\n")
  spots <- fread(paste(f,"spots.csv", sep ="_"))
  spots <- spots %>% select(-V1,-XM,-YM,-Slice)%>%
    separate(Label, into = c("image","ROI","channel"), sep = ":") %>%
    pivot_wider(names_from = "channel", values_from = c("Mean","Median","Min","Max","Mode","StdDev"))
  
  
  spots.pp <- SpatialPoints(coords = spots[,c("X","Y")])
  
  
  
  
  ## spots per fiber
  spot_in_fiber <- sp::over(x = spots.pp,y = p, returnList = F,)
  
  
  spot_in_fiber= data.frame(ROI_spot = spots$ROI, 
           fiber_id = names(p)[spot_in_fiber]) %>%
    filter(!is.na(fiber_id)) 
  
  fibers <- fibers %>% left_join(spot_in_fiber, by = c("ROI" = "fiber_id"))
  
  
  
  
  cat("\t\toverlapping spots and fibers\n") 
  csa[[f]] <- list("fibers"=fibers, "nuclei" = nuclei, "spots" = spots, "nuclei_rois" = p_nuclei, "fibers_rois" = p)
  
}



fibers <- lapply(csa, "[[", 1)
fibers = bind_rows(fibers, .id = "sample")

fibers = fibers %>% filter(Solidity >.7, AR < 4, Mean_1<5)


sample = files[2]

x = fortify(csa[[sample]]$fibers_rois)

x = x %>% inner_join(fibers %>% group_by_at(vars(1:last_col()-1)) %>% summarise(n = n_distinct(ROI_spot,na.rm=T)), by = c("id" = "ROI"))



ggplot(x, aes(long,lat,group = id, fill = n)) + 
  geom_polygon(col = "black", lwd = 1) + 
  geom_point(data = csa[[sample]]$spots, aes(X,Y),col = "red", inherit.aes = F) +
  scale_y_reverse() +
  theme_bw() + 
  scale_color_viridis_c() +
  ggtitle(sample)




count_spots = fibers %>%
  group_by(sample,ROI) %>%
  summarise(n = n_distinct(ROI_spot,na.rm = T)) %>% ungroup %>%
  count(sample,n,name="count") %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)* 100)

ggplot(count_spots, aes(n,prop,col=sample)) + geom_path() + geom_point()

