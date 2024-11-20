.BaBA <- function (animal, barrier, d, interval = NULL, b_time = 4, p_time = 36, 
          w = 168, tolerance = 0, units = "hours", max_cross = 0, sd_multiplier = 1, 
          exclude_buffer = F, export_images = F, img_path = "event_imgs", 
          img_suffix = NULL) 
{
  
  # b time
  b <- b_time/interval
  if (b < 1) 
    stop("interval needs to be set no bigger than b_time")
  if (round(b) != b) 
    stop("b_time must be divisible by interval")
  
  # p time
  p <- p_time/interval
  if (round(p) != p) 
    stop("p_time must be divisible by interval")
  
  # Assign point IDs
  animal <- animal[order(animal$Animal.ID, animal$date), ]
  ids <- unique(animal$Animal.ID)
  animal$ptsID <- NA
  for (id in ids) {
    inx <- animal$Animal.ID == id
    mov.seg.i <- animal[inx, ]
    animal$ptsID[inx] <- seq(nrow(mov.seg.i))
  }
  
  # create buffer 
  print("locating encounter events...")
  barrier_buffer <- st_buffer(barrier, dist=d)
  encounter <- st_intersection(animal, barrier_buffer)
  if (nrow(encounter) == 0)  { stop("no barrier encounter detected.")
  } else { message("barrier encounters detected!")}
  
  # loop over animals which encountered the fenceline
  ids <- unique(encounter$Animal.ID)
  for (id in ids) {
    message(id)
    encounter_i <- encounter[encounter$Animal.ID == id, ]
    if (!nrow(encounter_i)) {
      warning(paste0("Individual ", id, " has no locations overlapped with the barrier buffer and is eliminated from analysis."))
      (next)()
    }
    encounter_i$timediff <- c(interval, as.numeric(diff(encounter_i$date), units = units))
    encounter_i$timediff2 <- round(encounter_i$timediff - interval, digits = 1)
    
    if ( any(encounter_i$timediff2 > interval & 
            encounter_i$timediff2 <= tolerance, na.rm = T) ) {
      idx_pts_of_interest <- which(encounter_i$timediff2 > interval & 
                                   encounter_i$timediff2 <= tolerance)
      message('points of interest')
      for (pt in idx_pts_of_interest) {
        ptsID_of_interest_B <- encounter_i$ptsID[pt]
        ptsID_of_interest_A <- encounter_i$ptsID[pt - 1]
        inx.id <- animal$Animal.ID == id
        inx.A <- animal$ptsID > ptsID_of_interest_A
        inx.B <- animal$ptsID < ptsID_of_interest_B
        fetched_pts <- NULL
        fetched_pt <- animal[inx.id & inx.A & inx.B, ]
        
        if (!nrow(fetched_pt)) {
          encounter_i$timediff2[pt] <- 0
          (next)()
        } else {
          fetched_pt$timediff <- NA
          fetched_pt$timediff2 <- 0
          encounter_i$timediff2[pt] <- 0
          if (pt == idx_pts_of_interest[1])  fetched_pts <- fetched_pt
          if (is.null(fetched_pts)) fetched_pts <- fetched_pt
          else fetched_pts <- rbind(fetched_pts, fetched_pt)
        }
      }
      if (!is.null(fetched_pts)) {
        encounter_i <- rbind(encounter_i, fetched_pts)
      } else {
        message('-- no points to fetch!')
      }
      encounter_i <- encounter_i[order(encounter_i$ptsID), ]
    }
    encounter_i$burstID <- paste(id, cumsum(encounter_i$timediff2), sep = "_")
    if (id == unique(encounter$Animal.ID[1])) 
      encounter_complete <- encounter_i
    else encounter_complete <- rbind(encounter_complete, 
                                     encounter_i)
  }
  encounter <- encounter_complete
  
  
  
  # now to classify behaviors into categories
  print("classifying behaviors...")
  pb <- txtProgressBar(style = 3)
  event_df <- NULL
  barrier_sp <- barrier %>% as("Spatial")
  p4s <- CRS("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
  barrier_sp <- spTransform(barrier_sp, p4s)
  ids = unique(encounter$burstID)
  for (i in ids) {
    setTxtProgressBar(pb, which(ids == i)/length(ids))
    encounter_i <- encounter[encounter$burstID == i, ]
    encounter_i_sp <- spTransform(as(encounter_i, "Spatial"), p4s)
    animal_i <- animal[animal$Animal.ID == encounter_i$Animal.ID[1], ] 
    animal_i_sp <- spTransform(as(animal_i, "Spatial"), p4s)
    start_time <- encounter_i$date[1]
    end_time <- encounter_i$date[nrow(encounter_i)]
    duration <- difftime(end_time, start_time, units = units)
    straightness_i <- strtns(encounter_i_sp)
    coords_i <- st_coordinates(encounter_i)
    
    if (duration <= b_time) {
      pt.first <- encounter_i$ptsID[1]
      pt.last <- encounter_i$ptsID[nrow(encounter_i)]
      mov_seg_i <- movement.segment.b(animal_i_sp, pt.first, pt.last)
      int.num <- length(rgeos::gIntersection(mov_seg_i, barrier_sp))
      # if (int.num == 0 & nrow(coords_i) != (nrow(encounter_i) + 2)) {
      #   classification <- "unknown"
      # }
      # else {
        classification <- ifelse(int.num == 0, "Bounce", "Cross")
      # }
      if (export_images) {
        mov_seg_i_sf <- mov_seg_i %>% st_as_sf() %>% st_transform(st_crs(barrier))
        bb <- st_bbox(mov_seg_i_sf)
        p <- ggplot() + 
          geom_sf(data=barrier_buffer, fill='lightpink') + 
          geom_sf(data=barrier, color='red') +
          geom_sf(data=mov_seg_i_sf, color='black') +
          geom_sf(data=encounter_i, pch = 20, color = "cyan3") +
          coord_sf(xlim = c(bb[['xmin']], bb[['xmax']]),
                   ylim = c(bb[['ymin']], bb[['ymax']]))
        
        fname <- paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png")
        ggsave(file=fname, plot=p, width = 6, height = 6)
      }
    }
    if (duration > b_time) {
      mov_seg_i <- SpatialLines(list(Lines(Line(coords_i), 
                     ID = encounter_i$date[1])), proj4string = p4s)
      int.num <- length(rgeos::gIntersection(mov_seg_i, barrier_sp))
      if (duration > p_time) {
        classification <- "Trapped"
      }
      else {
        classification <- "TBD"
      }
      if (export_images & !classification %in% "TBD") {
        mov_seg_i_sf <- mov_seg_i %>% st_as_sf() %>% st_transform(st_crs(barrier))
        bb <- st_bbox(mov_seg_i_sf)
        p <- ggplot() + 
          geom_sf(data=barrier_buffer, fill='lightpink') + 
          geom_sf(data=barrier, color='red') +
          geom_sf(data=mov_seg_i_sf, color='black') +
          geom_sf(data=encounter_i, pch = 20, color = "cyan3") +
          coord_sf(xlim = c(bb[['xmin']], bb[['xmax']]),
                   ylim = c(bb[['ymin']], bb[['ymax']]))
        
        fname <- paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png")
        ggsave(file=fname, plot=p, width = 6, height = 6)
      }
    }
    event_df <- rbind(event_df, data.frame(AnimalID = encounter_i$Animal.ID[1], 
                                           burstID = i, 
                                           easting = coords_i[1, 1], 
                                           northing = coords_i[1, 2], 
                                           start_time, end_time, duration, cross = int.num, 
                                           straightness = ifelse(classification %in% c("Bounce", "Cross"), NA, straightness_i), 
                                           eventTYPE = classification, 
                                           stringsAsFactors = F))
  }
  close(pb)
  
  print('dealing with unknowns/TBD...')
  for (i in 1:nrow(event_df)) {
    # print(i)
    
    if (event_df[i, ]$eventTYPE == "TBD") {
      event_i <- event_df[i, ]
      duration_i <- event_i$duration
      straightness_i <- event_i$straightness
      animal_i <- animal[animal$Animal.ID == event_i$AnimalID, ]
      encounter_i <- encounter[encounter$burstID %in% event_i$burstID, ]
      if (exclude_buffer) {
        animal_i <- animal_i[!animal_i$ptsID %in% encounter$ptsID[encounter$Animal.ID == event_i$AnimalID], ]
      }
      animal_i <- animal_i[animal_i$date >= event_i$start_time - 
                             as.difftime(w/2, units = units) & animal_i$date <= 
                             event_i$end_time + as.difftime(w/2, units = units), ]
      animal_i$continuousID <- cumsum(abs(c(interval, 
                                            round(diff(animal_i$date, units = units), digits = 1)) - interval))
      straightnesses_i <- NULL
      for (ii in unique(animal_i$continuousID)) {
        animal_ii <- animal_i[animal_i$continuousID == ii, ]
        animal_ii_sp <- spTransform(as(animal_ii, "Spatial"), p4s)
        duration_ii <- difftime(animal_ii$date[nrow(animal_ii)], animal_ii$date[1], units = units)
        if (duration_ii >= duration_i) {
          for (iii in c(1:(which(animal_ii$date > (animal_ii$date[nrow(animal_ii)] - 
                                                   as.difftime(duration_i, units = units)))[1] - 
                           1))) {
            mov_seg <- animal_ii_sp[iii:(iii + duration_i/interval), ]
            straightnesses_i <- c(straightnesses_i, strtns(mov_seg))
          }
        }
      }
      LS <- length(straightnesses_i)
      if (LS >= (w/interval + 1)/4 && LS > 1) {
        upper <- mean(straightnesses_i) + sd_multiplier *  sd(straightnesses_i)
        lower <- mean(straightnesses_i) - sd_multiplier *  sd(straightnesses_i)
        if (straightness_i < lower) {
          event_df[i, ]$eventTYPE <- ifelse(event_i$cross < max_cross, "Back_n_forth", "unknown")
        }
        if (straightness_i > upper) {
          event_df[i, ]$eventTYPE <- ifelse(event_i$cross <  max_cross, "Trace", "unknown")
        }
        if (straightness_i >= lower & event_i$straightness <= upper) 
          event_df[i, ]$eventTYPE <- "Average"
      }
      else {
        event_df[i, ]$eventTYPE = "Trace_OR_Back_n_Forth"
        if (is.null(straightnesses_i)) {
          straightnesses_i <- NA
        }
      }
      # if (export_images) {
      #   mov_seg_i_sf <- mov_seg_i %>% st_as_sf() %>% st_transform(st_crs(barrier))
      #   bb <- st_bbox(mov_seg_i_sf)
      #   p <- ggplot() + 
      #     geom_sf(data=barrier_buffer, fill='lightpink') + 
      #     geom_sf(data=barrier, color='red') +
      #     geom_sf(data=mov_seg_i_sf, color='black') +
      #     geom_sf(data=encounter_i, pch = 20, color = "cyan3") +
      #     coord_sf(xlim = c(bb[['xmin']], bb[['xmax']]),
      #              ylim = c(bb[['ymin']], bb[['ymax']]))
      #   
      #   fname <- paste0(img_path, "/", classification, "_", i, "_", img_suffix, ".png")
      #   ggsave(file=fname, plot=p, width = 6, height = 6)
      # 
      #       width = 6, height = 6, res = 96, units = "in")
      #   A = by(as.data.frame(st_coordinates(animal_i)), 
      #          animal_i$continuousID, Line, simplify = T)
      #   A = SpatialLines(mapply(sp::Lines, A, ID = names(A), 
      #                           SIMPLIFY = F))
      #   plot(A, main = event_df[i, ]$eventTYPE, 
      #        sub = paste("cross = ", event_df[i, ]$cross, 
      #                    ", duration =", event_df[i, ]$duration, 
      #                    ", stri =", round(straightness_i, 2), 
      #                    ", str_mean = ", round(mean(straightnesses_i), 
      #                                                                                                                                         2)))
      #   plot(barrier_buffer, border=scales::alpha("red", 0.5), lty = "dashed", add = T)
      #   lines(barrier, col = "red", lwd = 2)
      #   points(encounter_i, pch = 20, col = "cyan3", type = "o", lwd = 2)
      #   dev.off()
      # }
    }
  }
  print("creating dataframe...")
  encounter <- encounter[!duplicated(encounter$burstID), ]
  encounter <- encounter[, c("Animal.ID", "burstID", "date")]
  encounter <- merge(encounter, event_df[, c("burstID", "eventTYPE")])
  return(list(encounters = encounter, classification = event_df))
}
