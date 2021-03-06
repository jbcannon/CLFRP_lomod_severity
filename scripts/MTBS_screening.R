#############################################################################################
#### Script for delineating low- to moderate-severity portions of wildfires from MTBS data
#### Author: Jeff Cannon (Jeffery.Cannon@colostate.edu)
#### Institution: Colorado Forest Restoration Institute, Colorado State University
#### Date Created: 07/25/2017
#### Last Modified: 08/10/2017
#############################################################################################
# The goal of this script is to input data on fire severity (boundary, 6-class severity) from
# MTBS and identify potential sampling areas composed of low- to mixed-severity fire. The
# script identifies all high severity patches, and eliminates those greater than a given
# threshold (e.g., 20 ha) from consideration. The script then delineates potential sampling
# areas and outputs maps and shapefiles illustrating potential sampling areas, and areas
# removed from consideration. The creation of this script was sponsored by the US Geological
# Survey  and the US Forest Service Collaborative Forest Landscape Restoraiton Program.
# The script works by taking a plot size (e.g., 400 ha square), overlays a sampling grid
# composed of smaller areas (e.g., 16-25 ha quadrats, grain size), and eliminates those areas
# that intersect with 'large' patches of high severity fire.
#############################################################################################

#######################################START SET UP##########################################
#---> Load libraries
packages <- c('raster', 'rgdal', 'rgeos', 'maptools', 'igraph')
for(package in packages){
  if(suppressMessages(!require(package,character.only=T))){
    install.packages(package,repos='https://cran.mtu.edu/')
    suppressMessages(library(package,character.only=T))
  }
}

#---> Choose parameters involving size of sampling units
plot_size <- 400 #ha
plot_grain <- plot_size / 4^2 #ha, default grain is 1/16th plot_size
edge_size <- sqrt(plot_grain * 10000)

#---> Maximum allowable size of high severity patch, 
max_hi_Sev <-  200000 #m2, 20 ha... must be same units (m) as raster resolution
neighbor_rule <- 4 # 4 or 8 neighboring cells used for clump detection

#---> Input link table to retreive filenames associated with each fire. This is an input fire list 
#     that corresponds to data included for each fire.
link_table <- read.csv('MTBS_FRONT_RANGE/filename_link.csv')

#---> Remove the following fires which have no suitable areas following analysis (i.e., throw an error)
#removal.list <- c('Fern Lake', 'Galena', 'Springer', 'Wetmore', 'Crystal Fire', 'Burning Tree', 'Fourmile Canyon')
#link_table <- subset(link_table, ! fire_name %in% removal.list)

########################################END SET UP###########################################

######################################START ANALYSIS#########################################
#---> Remove completed samples from link_table
complete <- list.files('potential_samples/', pattern = '*.pdf') #If pdf already exists...
complete <- substr(complete, start = 0, stop = nchar(complete) - nchar('_SAMPLING.pdf'))
link_table <- subset(link_table, ! folder_name %in% complete) #then remove from link table

#---> Loop through all fires in fire table. For each loop, perform analysis, and output shapefile and map 
for(fire_nm in link_table$fire_name)
{
  completed = list.files('potential_samples/')

  
  print(paste(fire_nm, ': process initiated', sep = ''))
  t0 <- Sys.time()
  
  #---> Get filenames for selected fire in fire_nm
  link_attrib <- subset(link_table, fire_name == fire_nm) 
  if(link_attrib$folder_name %in% completed) next
  severity_filename <- paste('MTBS_FRONT_RANGE/', link_attrib$folder_name, '/', link_attrib$file_base, '_dnbr6.tif', sep = '')
  boundary_filename <- paste('MTBS_FRONT_RANGE/', link_attrib$folder_name, '/', link_attrib$file_base, '_burn_bndy.shp', sep = '')
  
  #---> Use filenames to load severity raster and boundary
  print(paste('Loading ', fire_nm, ' data.', sep = ''))
  fire_sev <- raster(severity_filename)
  fire_bnd <- readOGR(dsn = boundary_filename, layer = paste(link_attrib$file_base, '_burn_bndy', sep = ''), verbose = FALSE)
  
  #---> create raster of continous patches of hi severity fire
  print('Identifying hi severity areas')
  fire_hi <- fire_sev == 4 #subset high severity only
  fire_hi_clump <- clump(fire_hi, directions = neighbor_rule) #clump detection of hi sev.
  fire_hi_clump_list <- as.data.frame(freq(fire_hi_clump)) #create list of patch numbers and pixel count
  fire_hi_clump_list <- subset(fire_hi_clump_list, value != 'NA') #drop NA from list
  fire_hi_clump_list$area_m2 <- fire_hi_clump_list$count * prod(res(fire_sev)) # measure size of each patch from px count
  fire_hi_clump_list <- subset(fire_hi_clump_list, area_m2 >= max_hi_Sev)  #keep only large hi severity patches from list
  fire_largeHiSev <- fire_hi_clump %in% fire_hi_clump_list$value
  
  #---> create exclusion raster, eliminated NA and high severity patches
  fire_NA <- fire_sev == 6 #identify areas where severity is not known (cloud/SLC failure)
  sample_exclude <- fire_NA | fire_largeHiSev #get areas with NA or large severe patch
  sample_exclude <- sample_exclude * 1 #create raster grid, 1 = exclusion area
  
  sample_exclude = crop(sample_exclude, fire_bnd)
  sample_exclude = rasterToPolygons(sample_exclude, dissolve=TRUE)
  sample_exclude = subset(sample_exclude, layer == 1)
  
  if(length(sample_exclude)> 0) {
    #---> create fishnet on top of fire perimeter using spatial grain
    grid_r <- raster(extent(fire_bnd), resolution = edge_size, crs = proj4string(fire_bnd))
    grid_r[] <- 1:ncell(grid_r)
    grid_shp <- rasterToPolygons(grid_r)
    grid_cells_out = grid_shp[sample_exclude,]$layer
    grid_out = subset(grid_shp, ! layer %in% grid_cells_out)
    
    grid_out = gBuffer(grid_out,width=0); fire_bnd = gBuffer(fire_bnd, width=0)
    grid_sel = intersect(grid_out, fire_bnd)
    grid_sel <- gUnaryUnion(grid_sel)
    grid_sel <- sp::disaggregate(grid_sel)
    grid_sel$area_ha <- gArea(grid_sel, byid = TRUE) / 10000
    grid_sel <- subset(grid_sel, area_ha >= plot_size)
  } else grid_sel = fire_bnd
  
  #######################################END ANALYSIS##########################################
  
  ######################################START OUTPUTS##########################################
  #---> output shapefile of sampling areas
  shp_path <- paste('potential_samples/', link_attrib$folder_name, '_SAMPLING.shp', sep = '')
  print(c('Exporting sampling area shapefile to:', shp_path))
  layer_nm <- paste(link_attrib$folder_name, '_SAMPLING', sep = '')
  writeOGR(grid_sel, shp_path, layer = layer_nm, driver = 'ESRI Shapefile', overwrite = TRUE)
  ######################################END OUTPUTS ##########################################
  
  ######################################START GRAPHICS#########################################
  #--> Output plot map containing exclusion areas (hatched) and sample areas (outlined)
  print('Exporting map')
  img_path <- paste('potential_samples/', link_attrib$folder_name, '_SAMPLING.pdf', sep = '')
  #create exclusion areas
  pdf(img_path, width = 7, height = 7)
  if(length(grid_sel) > 0) {sampling_ha <- round(gArea(grid_sel) / 10000)} else {sampling_ha = 0}
  fire_ha <- round(gArea(fire_bnd) / 10000)
  plot(fire_bnd, main = c(paste(fire_nm, ' (', fire_ha, ' ha)' , sep = ''), paste('sample = ', sampling_ha, ' ha', sep = '')))
  plot(fire_sev, add = TRUE)
  #plot(out, add = TRUE, density = 30)
  plot(fire_bnd, add = TRUE, lwd = 3)
  plot(grid_sel, border = 'white', add = TRUE, lwd = 5)
  plot(grid_sel, border = 'black', lty = 2, lwd = 1, add = TRUE)
  axis(1); axis(2); box()
  scalebar(d = 5000, type = 'line', below = 'meters', xy=extent(fire_bnd)[c(1,3)], col = c('white', 'red'))
  dev.off()
  
  t1 <- Sys.time()
  print(paste(fire_nm, ' complete', sep = ''))
  print(t1-t0)
}

  #---> Merge all files into one pdf using pdf-toolkit
  dir <- paste(getwd(), 'potential_samples', sep = '/')
  pdfs <- list.files('potential_samples', pattern = '*.pdf')
  pdfs <- sort(pdfs)
  pdfs <- paste(dir, "/", pdfs, sep = '')
  pdfs <- paste(as.vector(pdfs), sep = ' ', collapse = ' ')
  cmd <- paste('pdftk ', pdfs, ' cat output ', dir, '/all_combined.pdf', sep = '')
  unlink(paste(dir, '/all_combined.pdf', sep = '')) # delete older combined PDF
  system(cmd)
