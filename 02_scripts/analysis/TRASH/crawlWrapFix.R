.crawlWrap <- function (obsData, timeStep = 1, ncores = 1, retryFits = 0, retrySD = 1, 
          retryParallel = FALSE, mov.model = ~1, err.model = NULL, 
          activity = NULL, drift = NULL, coord = c("x", "y"), proj = NULL, 
          Time.name = "time", time.scale = "hours", theta, fixPar, 
          method = "L-BFGS-B", control = NULL, constr = NULL, prior = NULL, 
          need.hess = TRUE, initialSANN = list(maxit = 200), attempts = 1, 
          predTime = NULL, fillCols = FALSE, coordLevel = NULL, ...) 
{
  if (is.data.frame(obsData)) {
    if (!inherits(obsData, "sf")) {
      if (!is.character(coord) | length(coord) != 2) 
        stop("coord must be character vector of length 2")
      if (any(!(c("ID", Time.name, coord) %in% names(obsData)))) 
        stop("obsData is missing ", paste(c("ID", Time.name, 
                                            coord)[!(c("ID", Time.name, coord) %in% names(obsData))], 
                                          collapse = ","))
      if (any(coord %in% c("mu.x", "nu.x", "mu.y", "nu.y", 
                           "se.mu.x", "se.nu.x", "se.mu.y", "se.nu.y", "speed"))) 
        stop("coordinates cannot include the following names: mu.x, nu.x, mu.y, nu.y, se.mu.x, se.nu.x, se.mu.y, se.nu.y, or speed \n   please choose different coord names")
    }
    else {
      if (any(!(c("ID", Time.name) %in% names(obsData)))) 
        stop("obsData is missing ", paste(c("ID", Time.name)[!(c("ID", 
                                                                 Time.name) %in% names(obsData))], collapse = ","))
    }
  }
  else if (inherits(obsData, "SpatialPoints")) {
    if (any(!(c("ID", Time.name) %in% names(obsData)))) 
      stop("obsData is missing ", paste(c("ID", Time.name)[!(c("ID", 
                                                               Time.name) %in% names(obsData))], collapse = ","))
    coord <- colnames(sp::coordinates(obsData))
  }
  else stop("obsData must be a data frame or a SpatialPointsDataFrame")
  if (retryFits < 0) 
    stop("retryFits must be non-negative")
  if (attempts < 1) 
    stop("attempts must be >=1")
  ids = as.character(unique(obsData$ID))
  ind_data <- list()
  hierInd <- FALSE
  for (i in ids) {
    ind_data[[i]] = obsData[which(obsData$ID == i), ]
    if (any(is.na(ind_data[[i]][[Time.name]]))) 
      stop("obsData$", Time.name, " cannot contain missing values")
    if (!is.null(suppressWarnings(ind_data[[i]]$level))) {
      if (is.factor(ind_data[[i]]$level) & is.null(coordLevel)) 
        stop("'level' field cannot be a factor unless 'coordLevel' is specified")
      if (!is.factor(ind_data[[i]]$level) & !is.null(coordLevel)) 
        stop("'level' field must be a factor when 'coordLevel' is specified")
      if (!is.character(coordLevel) | length(coordLevel) != 
          1) 
        stop("coordLevel must be a character string")
      if (!(coordLevel %in% levels(ind_data[[i]]$level))) 
        stop("'coordLevel' not found in 'level' field")
      ind_data[[i]] <- obsData[which(obsData$ID == i & 
                                       obsData$level == coordLevel), ]
      hierInd <- TRUE
    }
    else if (!is.null(coordLevel)) 
      stop("coordLevel can not be specified unless obsData includes a 'level' field")
  }
  nbAnimals <- length(ids)
  crawlArgs <- checkCrawlArgs(ids, nbAnimals, ind_data, obsData, 
                              Time.name, retryFits, attempts, timeStep, mov.model, 
                              err.model, activity, drift, theta, fixPar, retrySD, constr, 
                              prior, proj, predTime)
  id <- NULL
  cat("Fitting", nbAnimals, "track(s) using crawl::crwMLE...\n")
  if (ncores > 1) {
    for (pkg in c("doFuture", "future")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"", pkg, "\" needed for parallel processing to work. Please install it.", 
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    future::plan(future::multisession, workers = ncores)
    on.exit(with(oldDoPar, foreach::setDoPar(fun = fun, data = data, 
                                             info = info)), add = TRUE)
  }
  else {
    doParallel::registerDoParallel(cores = ncores)
  }
  if (missing(theta)) 
    theta <- NULL
  if (missing(fixPar)) 
    fixPar <- NULL
  withCallingHandlers(model_fits <- foreach(id = ind_data, 
                                            i = ids, .errorhandling = "pass", .packages = "crawl", 
                                            .final = function(x) stats::setNames(x, ids)) %dorng% 
                        {
                          cat("Individual ", i, "...\n", sep = "")
                          fit <- .crwMLE(data = id, mov.model = crawlArgs$mov.model[[i]], 
                                               err.model = crawlArgs$err.model[[i]], activity = crawlArgs$activity[[i]], 
                                               drift = crawlArgs$drift[[i]], coord = coord, 
                                               proj = crawlArgs$proj[[i]], Time.name = Time.name, 
                                               time.scale = time.scale, theta = crawlArgs$theta[[i]], 
                                               fixPar = crawlArgs$fixPar[[i]], method = method, 
                                               control = control, constr = crawlArgs$constr[[i]], 
                                               prior = crawlArgs$prior[[i]], need.hess = need.hess, 
                                               initialSANN = initialSANN, attempts = attempts, 
                                               ... = ...)
                        }, warning = muffleRNGwarning)
  if (ncores == 1) 
    doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  rm(ind_data)
  convFits <- ids[which(unlist(lapply(model_fits, function(x) inherits(x, "crwFit"))))]
  if (!length(convFits)) 
    stop("crawl::crwMLE failed for all individuals.  Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
  cat("DONE\n")
  if ((ncores == 1 & !retryFits) | retryFits | ncores > 1) 
    cat("\n")
  if (retryParallel & ncores > 1) {
    for (i in ids) {
      if (!inherits(model_fits[[i]], "crwFit")) 
        warning("crawl::crwMLE for individual ", i, " failed;\n", 
                model_fits[[i]], "   Check crawl::crwMLE arguments and/or consult crawl documentation.")
    }
  }
  if (ncores > 1 & nbAnimals > 1 & retryFits & retryParallel) {
    doFuture::registerDoFuture()
    future::plan(future::multisession, workers = ncores)
    if (retryParallel) 
      cat("Attempting to achieve convergence and valid variance estimates for each individual in parallel.\n    Press 'esc' to force exit from 'crawlWrap'... \n", 
          sep = "")
  }
  else {
    doParallel::registerDoParallel(cores = 1)
  }
  withCallingHandlers(model_fits <- foreach(mf = model_fits, 
                                            i = ids, .export = c("quietCrawl"), .errorhandling = "pass", 
                                            .packages = "crawl", .final = function(x) stats::setNames(x, 
                                                                                                      ids)) %dorng% {
                                                                                                        if (inherits(mf, "crwFit")) {
                                                                                                          if ((mf$convergence | any(is.na(mf$se[which(is.na(crawlArgs$fixPar[[i]]))]))) | 
                                                                                                              retryFits) {
                                                                                                            if (retryFits) {
                                                                                                              fitCount <- 0
                                                                                                              fitPar <- mf$estPar
                                                                                                              if (mf$convergence | any(is.na(mf$se[which(is.na(crawlArgs$fixPar[[i]]))]))) {
                                                                                                                if (mf$convergence) 
                                                                                                                  cat("\ncrawl::crwMLE for individual", i, 
                                                                                                                      "has suspect convergence: ", mf$message, 
                                                                                                                      "\n")
                                                                                                                if (any(is.na(mf$se[which(is.na(crawlArgs$fixPar[[i]]))]))) 
                                                                                                                  cat("\ncrawl::crwMLE for individual", i, 
                                                                                                                      "has NaN variance estimate(s)\n")
                                                                                                                cat("Attempting to achieve convergence and valid variance estimates for individual ", 
                                                                                                                    i, ". Press 'esc' to force exit from 'crawlWrap'\n", 
                                                                                                                    sep = "")
                                                                                                              }
                                                                                                              else {
                                                                                                                cat("Attempting to improve fit for individual ", 
                                                                                                                    i, ". Press 'esc' to force exit from 'crawlWrap'\n", 
                                                                                                                    sep = "")
                                                                                                              }
                                                                                                              while (fitCount < retryFits) {
                                                                                                                cat("\r    Attempt ", fitCount + 1, " of ", 
                                                                                                                    retryFits, " -- current log-likelihood value: ", 
                                                                                                                    mf$loglik, "  ...", sep = "")
                                                                                                                tmpFun <- function() {
                                                                                                                  tryCatch(suppressWarnings(suppressMessages(crawl::crwMLE(data = mf$data, 
                                                                                                                                                                           mov.model = crawlArgs$mov.model[[i]], 
                                                                                                                                                                           err.model = crawlArgs$err.model[[i]], 
                                                                                                                                                                           activity = crawlArgs$activity[[i]], drift = crawlArgs$drift[[i]], 
                                                                                                                                                                           coord = coord, proj = crawlArgs$proj[[i]], 
                                                                                                                                                                           Time.name = Time.name, time.scale = time.scale, 
                                                                                                                                                                           theta = fitPar + rnorm(length(fitPar), 
                                                                                                                                                                                                  0, crawlArgs$retrySD[[i]]), fixPar = crawlArgs$fixPar[[i]], 
                                                                                                                                                                           method = method, control = control, constr = crawlArgs$constr[[i]], 
                                                                                                                                                                           prior = crawlArgs$prior[[i]], need.hess = need.hess, 
                                                                                                                                                                           initialSANN = list(maxit = 0, trace = 0), 
                                                                                                                                                                           attempts = 1, ... = ...))), error = function(e) {
                                                                                                                                                                             e
                                                                                                                                                                           })
                                                                                                                }
                                                                                                                tmp <- NULL
                                                                                                                if (retryParallel) {
                                                                                                                  tmp <- tmpFun()
                                                                                                                }
                                                                                                                else {
                                                                                                                  tmp <- quietCrawl(tmpFun())
                                                                                                                }
                                                                                                                if (inherits(tmp, "crwFit")) {
                                                                                                                  if (tmp$convergence == 0) {
                                                                                                                    if (tmp$loglik > mf$loglik | all(!is.na(tmp$se[which(is.na(crawlArgs$fixPar[[i]]))]))) 
                                                                                                                      fitPar <- tmp$estPar
                                                                                                                    if (tmp$loglik >= mf$loglik & all(!is.na(tmp$se[which(is.na(crawlArgs$fixPar[[i]]))]))) 
                                                                                                                      mf <- tmp
                                                                                                                  }
                                                                                                                  rm(tmp)
                                                                                                                }
                                                                                                                fitCount <- fitCount + 1
                                                                                                              }
                                                                                                              if (mf$convergence | any(is.na(mf$se[which(is.na(crawlArgs$fixPar[[i]]))]))) {
                                                                                                                message("FAILED\n")
                                                                                                              }
                                                                                                              else {
                                                                                                                cat("DONE\n")
                                                                                                              }
                                                                                                            }
                                                                                                            else {
                                                                                                              if (mf$convergence) 
                                                                                                                warning("crawl::crwMLE for individual ", 
                                                                                                                        i, " has suspect convergence: ", mf$message)
                                                                                                              if (any(is.na(mf$se[which(is.na(crawlArgs$fixPar[[i]]))]))) 
                                                                                                                warning("crawl::crwMLE for individual ", 
                                                                                                                        i, " has NaN variance estimate(s)")
                                                                                                            }
                                                                                                          }
                                                                                                        }
                                                                                                        else {
                                                                                                          warning("\ncrawl::crwMLE for individual ", i, " failed;\n", 
                                                                                                                  mf, "   Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
                                                                                                        }
                                                                                                        mf
                                                                                                      }, warning = muffleRNGwarning)
  if (!(ncores > 1 & nbAnimals > 1 & retryFits & retryParallel)) {
    doParallel::stopImplicitCluster()
  }
  else future::plan(future::sequential)
  if (ncores > 1) {
    doFuture::registerDoFuture()
    future::plan(future::multisession, workers = ncores)
  }
  else {
    doParallel::registerDoParallel(cores = ncores)
  }
  convFits <- ids[which(unlist(lapply(model_fits, function(x) inherits(x, 
                                                                       "crwFit"))))]
  if (!length(convFits)) 
    stop("crawl::crwMLE failed for all individuals.  Check crawl::crwMLE arguments and/or consult crawl documentation.\n")
  model_fits <- model_fits[convFits]
  txt <- NULL
  if (inherits(obsData[[Time.name]], "POSIXct")) {
    td <- list()
    for (i in convFits) {
      td[[i]] <- crawlArgs$predTime[[i]]
      if (length(crawlArgs$predTime[[i]]) > 1) {
        td[[i]] <- utils::capture.output(difftime(crawlArgs$predTime[[i]][2], 
                                                  crawlArgs$predTime[[i]][1], units = "auto"))
        td[[i]] <- substr(td[[i]], 20, nchar(td[[i]]))
      }
    }
    if (length(unique(td)) == 1) 
      txt <- paste("at", td[[1]], "time steps")
  }
  if (retryFits > 0) 
    cat("\n")
  cat("Predicting locations (and uncertainty)", txt, "for", 
      length(convFits), "track(s) using crawl::crwPredict... ")
  if (time.scale %in% c("hours", "hour")) {
    ts <- 60 * 60
  }
  else if (time.scale %in% c("days", "day")) {
    ts <- 60 * 60 * 24
  }
  else if (time.scale %in% c("sec", "secs", "second", "seconds")) {
    ts <- 1
  }
  else if (time.scale %in% c("min", "mins", "minute", "minutes")) {
    ts <- 60
  }
  withCallingHandlers(predData <- foreach(mf = model_fits[convFits], 
                                          i = convFits, .combine = rbind, .errorhandling = "remove", 
                                          .packages = "crawl") %dorng% {
                                            pD <- crawl::crwPredict(mf, predTime = crawlArgs$predTime[[i]][crawlArgs$predTime[[i]] >= 
                                                                                                             min(mf$data[[Time.name]])], return.type = "flat")
                                            if (inherits(mf$data[[Time.name]], "POSIXct") && attributes(pD[[Time.name]])$tzone != 
                                                attributes(mf$data[[Time.name]])$tzone) {
                                              if (!requireNamespace("lubridate", quietly = TRUE)) {
                                                stop("Package \"lubridate\" needed for this function to work. Please install it.", 
                                                     call. = FALSE)
                                              }
                                              pD[[Time.name]] <- lubridate::with_tz(pD[[Time.name]], 
                                                                                    tz = attributes(mf$data[[Time.name]])$tzone)
                                            }
                                            if (length(crawlArgs$predTime[[i]]) == 1 && inherits(mf$data[[Time.name]], 
                                                                                                 "POSIXct")) {
                                              pD[[Time.name]] <- as.POSIXct(pD$TimeNum * ts, origin = "1970-01-01 00:00:00", 
                                                                            tz = attributes(pD[[Time.name]])$tzone)
                                            }
                                            if (!fillCols) {
                                              dups <- duplicated(mf$data[[Time.name]])
                                              tmpind_data <- as.data.frame(mf$data[!dups, , drop = FALSE])
                                              for (j in names(pD)[names(pD) %in% names(mf$data)]) {
                                                if (!(j %in% c(Time.name, "ID", coord))) {
                                                  if (!isTRUE(all.equal(pD[[j]], tmpind_data[[j]]))) {
                                                    pD[[j]][pD[[Time.name]] %in% tmpind_data[[Time.name]]] <- tmpind_data[[j]]
                                                    pD[[j]][!(pD[[Time.name]] %in% tmpind_data[[Time.name]])] <- NA
                                                  }
                                                }
                                              }
                                            }
                                            if (!is.null(coordLevel)) 
                                              pD$level <- coordLevel
                                            pD
                                          }, warning = muffleRNGwarning)
  if (ncores == 1) 
    doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  if (hierInd) {
    pData <- predData
    predData <- data.frame()
    for (i in convFits) {
      ipData <- pData[which(pData$ID == i), ]
      ipData$level <- factor(ipData$level, levels = levels(obsData$level))
      tmpData <- obsData[which(obsData$ID == i), ]
      tmpData <- merge(tmpData, ipData[, c(Time.name, "level"), 
                                       drop = FALSE], all = TRUE, by = c(Time.name, 
                                                                         "level"))
      tmpData$level[is.na(tmpData$level)] <- coordLevel
      tmpData$ID[is.na(tmpData$ID)] <- i
      for (jj in names(tmpData)[!(names(tmpData) %in% names(ipData))]) {
        ipData[[jj]] <- rep(NA, nrow(ipData))
      }
      for (jj in names(ipData)[!(names(ipData) %in% names(tmpData))]) {
        tmpData[[jj]] <- rep(NA, nrow(tmpData))
      }
      ipData <- ipData[, names(tmpData)]
      tmpInd <- which(tmpData$time %in% ipData$time & tmpData$level == 
                        coordLevel)
      tmpData[tmpInd, ] <- Map(function(x, y) {
        x[is.na(x)] <- y[is.na(x)]
        x
      }, tmpData[tmpInd, ], ipData[ipData$time %in% tmpData$time, 
      ])
      predData <- rbind(predData, tmpData)
    }
    predData <- predData[, names(pData)]
    attrNames <- names(attributes(pData))[!(names(attributes(pData)) %in% 
                                              names(attributes(predData)))]
    attributes(predData)[attrNames] <- attributes(pData)[attrNames]
    attr(predData, "coordLevel") <- coordLevel
    class(predData) <- append(c("crwPredict", "hierarchical"), 
                              class(predData))
    cat("DONE\n")
    return(crwHierData(list(crwFits = model_fits, crwPredict = predData)))
  }
  else {
    cat("DONE\n")
    return(crwData(list(crwFits = model_fits, crwPredict = predData)))
  }
}


.crwMLE <- function(
    data,
    mov.model = ~ 1,
    err.model = NULL,
    activity = NULL,
    drift = FALSE,
    coord = c("x", "y"),
    proj = NULL,
    Time.name = "time",
    time.scale = NULL,
    theta = NULL,
    fixPar = NULL,
    method = "Nelder-Mead",
    control = NULL,
    constr = list(lower = -Inf, upper = Inf),
    prior = NULL,
    need.hess = TRUE,
    initialSANN = list(maxit = 200),
    attempts = 1,
    retrySD = 1,
    skip_check = FALSE,
    ...
)

{
  
  if (!inherits(data, c("data.frame", "tbl_df"))) {
    stop("data must be an 'sf'/'sp' spatial object or a data.frame/tibble")
  }
  
  p4 <- NULL
  if (!is.null(proj) && !inherits(proj, "crs")) {
    if (inherits(proj, "numeric") || inherits(proj, "character")) {
      p4 <- sf::st_crs(proj)
    } else {
      stop(
        "provided projection is not valid. must be an EPSG integer code or proj4string character"
      )
    }
  }
  if (inherits(proj, "crs")) {
    p4 <- proj
  }
  
  if(!is.null(p4) && p4$IsGeographic) {
    stop(
      "Provided projection is geographic (e.g. longlat). Coordinates must be in
      a non-geographic, projected coordinate reference system."
    )
  }
  
  if (inherits(data[[Time.name]], "POSIXct")) {
    if (is.null(time.scale)) {
      time.scale = crawl::detect_timescale(data[[Time.name]])
    }
    if (time.scale %in% c("hours", "hour")) {
      ts = 60 * 60
    } else if (time.scale %in% c("days", "day")) {
      ts = 60 * 60 * 24
    } else if (time.scale %in% c("sec", "secs", "second", "seconds")) {
      ts = 1
    } else if (time.scale %in% c("min", "mins", "minute", "minutes")) {
      ts = 60
    } else
      stop("'time.scale' not specified correctly!")
    data$TimeNum <- as.numeric(data[[Time.name]]) / ts
  } else{
    data$TimeNum <- as.numeric(data[[Time.name]])
    ts = 1
  }
  st <- Sys.time()
  ## SET UP MODEL MATRICES AND PARAMETERS ##
  errMod <- !is.null(err.model)
  #if(!errMod) stop("Error model must be specified! (argument 'err.model' is currently set to NULL)")
  activeMod <- !is.null(activity)
  driftMod <- drift
  # if (
  #   length(initial.state$a) != 2*(driftMod+2) | all(dim(initial.state$P) != c(2*(driftMod+2), 2*(driftMod+2)))
  # ) stop("Dimentions of 'initial.state' argument are not correct for the specified model")
  mov.mf <-
    model.matrix(mov.model, model.frame(mov.model, data, na.action = na.pass))
  if (any(is.na(mov.mf)))
    stop("Missing values are not allowed in movement covariates!")
  n.mov <- ncol(mov.mf)
  if (errMod) {
    err.mfX <-
      model.matrix(err.model$x,
                   model.frame(err.model$x, data, na.action = na.pass))
    err.mfX <- ifelse(is.na(err.mfX), 0, err.mfX)
    n.errX <- ncol(err.mfX)
    if (!is.null(err.model$y)) {
      err.mfY <-
        model.matrix(err.model$y,
                     model.frame(err.model$y, data, na.action = na.pass))
      err.mfY <- ifelse(is.na(err.mfY), 0, err.mfY)
      n.errY <- ncol(err.mfY)
    } else {
      err.mfY <- NULL
      n.errY <- 0
    }
    if (!is.null(err.model$rho)) {
      rho = model.matrix(err.model$rho,
                         model.frame(err.model$rho, data, na.action = na.pass))[, -1]
      if (any(rho > 1 |
              rho < -1, na.rm = TRUE))
        stop("Error model correlation outside of the range (-1, 1).")
      rho <- ifelse(is.na(rho), 0, rho)
    } else
      rho = NULL
  } else {
    n.errY <- n.errX <- 0
    err.mfX <- err.mfY <- rho <- NULL
  }
  if (activeMod) {
    #stop.model
    activity <-
      model.matrix(activity, model.frame(activity, data, na.action = na.pass))
    if (ncol(activity) > 2)
      stop("There can only be one activity variable.")
    activity <- as.double(activity[, 2])
    if (any(activity < 0) |
        any(activity > 1))
      stop("'activity' variable must be >=0 and <=1.")
    if (any(is.na(activity)))
      stop("Missing values are not allowed in the activity variable.")
  } else
    activity <- NULL
  n.drift <- as.integer(driftMod)
  n.activ <- as.integer(activeMod)
  b.nms <- paste("ln beta ", colnames(mov.mf), sep = "")
  sig.nms <- paste("ln sigma ", colnames(mov.mf), sep = "")
  if (errMod) {
    if (!is.null(err.model$y)) {
      tau.nms <- c(paste("ln tau.x ", colnames(err.mfX), sep = ""),
                   paste("ln tau.y ", colnames(err.mfY), sep = ""))
    } else
      tau.nms <- paste("ln tau ", colnames(err.mfX), sep = "")
  } else
    tau.nms <- NULL
  if (activeMod) {
    active.nms <- "ln phi"
  } else
    active.nms <- NULL
  if (driftMod) {
    drift.nms <- c("ln sigma.drift/sigma", "ln psi-1")
  } else
    drift.nms <- NULL
  nms <- c(tau.nms, sig.nms, b.nms, active.nms, drift.nms)
  n.par <- length(nms)
  if (is.null(fixPar))
    fixPar <- rep(NA, n.par)
  n.theta = sum(is.na(fixPar))
  if (length(fixPar) != n.par)
    stop(
      "'fixPar' argument is not the right length! The number of parameters in the model is ",
      n.par,
      "\n"
    )
  if (!(length(constr$lower) == 1 |
        length(constr$lower) == sum(is.na(fixPar))))
    stop(
      "The number of lower contraints specified is not correctly! The number of free parameters is ",
      sum(is.na(fixPar)),
      "\n"
    )
  if (!(length(constr$upper) == 1 |
        length(constr$upper) == sum(is.na(fixPar))))
    stop(
      "The number of upper contraints specified is not correctly! The number of free parameters is ",
      sum(is.na(fixPar)),
      "\n"
    )
  # if(length(constr$upper)==1) constr$upper <- rep(constr$upper, sum(is.na(fixPar)))
  # if(length(constr$lower)==1) constr$lower <- rep(constr$lower, sum(is.na(fixPar)))
  if (is.null(theta))
    theta <- rep(0, n.theta)
  theta[theta < constr$lower] = constr$lower[theta < constr$lower] + 0.01
  theta[theta > constr$upper] = constr$upper[theta > constr$upper] - 0.01
  if (driftMod &
      is.na(fixPar[n.par]))
    theta[sum(is.na(fixPar))] <- log(diff(range(data$TimeNum)) / 9)
  if (length(theta) != n.theta) {
    stop("\nWrong number of parameters specified in start value.\n")
  }
  
  y = data[, c(coord[1], coord[2])]
  noObs <- as.numeric(is.na(y[, 1]) | is.na(y[, 2]))
  y[noObs == 1, ] = 0
  
  checkFit <- 1
  thetaAttempt <- theta
  if (!is.null(initialSANN) & method != 'SANN') {
    message("Beginning SANN initialization ...")
    init <-
      optim(
        thetaAttempt,
        crwN2ll,
        method = 'SANN',
        control = initialSANN,
        fixPar = fixPar,
        y = y,
        noObs = noObs,
        delta = c(diff(data$TimeNum), 1),
        #a=initial.state$a, P=initial.state$P,
        mov.mf = mov.mf,
        err.mfX = err.mfX,
        err.mfY = err.mfY,
        rho = rho,
        activity = activity,
        n.mov = n.mov,
        n.errX = n.errX,
        n.errY = n.errY,
        driftMod = driftMod,
        prior = prior,
        need.hess = FALSE,
        constr = constr
      )
    print(init)
    #thetaAttempt <- init$par
  } else
    init <- list(par = thetaAttempt)
  #if(any(init$par<lower)) init$par[init$par<lower] <- lower[init$par<lower] + 0.000001
  #if(any(init$par>upper)) init$par[init$par>upper] <- upper[init$par>upper] - 0.000001
  message("Beginning likelihood optimization ...")
  myPar <- init$par
  while (attempts > 0 & checkFit) {
      mle = tryCatch(
        {
          message("\n-----\nattempt #", attempts)
          optim(
              init$par,
              crwN2ll,
              method = method,
              hessian = need.hess,
              lower = constr$lower,
              upper = constr$upper,
              control = control,
              fixPar = fixPar,
              y = y,
              noObs = noObs,
              delta = c(diff(data$TimeNum), 1),
              mov.mf = mov.mf,
              err.mfX = err.mfX,
              err.mfY = err.mfY,
              activity = activity,
              n.mov = n.mov,
              n.errX = n.errX,
              n.errY = n.errY,
              rho = rho,
              driftMod = driftMod,
              prior = prior,
              need.hess = need.hess,
              constr = constr
            )
        },
        error = function(cond) {
          message("Here's the original error message:")
          message(conditionMessage(cond))
          # Choose a return value in case of error
          NA
        })
    attempts <- attempts - 1
    checkFit <- check_fit(mle)
    myPar <- myPar + rnorm(length(myPar), 0, retrySD)
    init$par <- myPar
  }
  
  if (checkFit & !skip_check) {
    return(simpleError(
      paste(
        "crwMLE failed. Try increasing attempts or changing user",
        "provided sd values. Other parameter values may need to be changed.",
        "Good Luck!"
      )
    ))
  }
  par <- fixPar
  par[is.na(fixPar)] <- mle$par
  Cmat <- matrix(NA, n.par, n.par)
  Cmat[is.na(fixPar), is.na(fixPar)] <- 2 * solve(mle$hessian)
  se <- sqrt(diag(Cmat))
  ci.l <- par - 1.96 * se
  ci.u <- par + 1.96 * se
  
  out <-
    list(
      par = par,
      estPar = mle$par,
      se = se,
      ci = cbind(L = ci.l, U = ci.u),
      Cmat = Cmat,
      loglik = -mle$value / 2,
      aic = mle$value + 2 * sum(is.na(fixPar)),
      coord = coord,
      fixPar = fixPar,
      convergence = mle$convergence,
      message = mle$message,
      activity = activity,
      random.drift = drift,
      mov.model = mov.model,
      err.model = err.model,
      n.par = n.par,
      nms = nms,
      n.mov = n.mov,
      n.errX = n.errX,
      n.errY = n.errY,
      mov.mf = mov.mf,
      err.mfX = err.mfX,
      err.mfY = err.mfY,
      rho = rho,
      Time.name = Time.name,
      init = init,
      data = data,
      lower = constr$lower,
      upper = constr$upper,
      prior = prior,
      need.hess = need.hess,
      runTime = difftime(Sys.time(), st)
    )
  attr(out, "epsg") <- p4$epsg
  attr(out, "proj4") <- p4$proj4string
  attr(out, "time.scale") = ts
  class(out) <- c("crwFit")
  if (drift) {
    class(out) <- c("crwFit_drift", "crwFit")
  }
  return(out)
}



