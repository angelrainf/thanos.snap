# THANOS - 25 body-size classes, 1000 iterations, four mortality scenarios.
# Draw 30 random, non-overlapping quadrats occupied before disturbance.
# Hold locations fixed across scenarios, including quadrats emptied afterward.

# Avoid nested BLAS/OpenMP threading inside parallel R workers.
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)

library(dplyr)
library(tidyr)
library(ggplot2)
library(mobsim)
library(vegan)
library(ecopart)
library(spatstat.random)
library(spatstat.geom)

# ---- 1. SETTINGS --------------------------------------------------------------
cfg <- list(
  seed=404,n_iter=1000,n_size_classes=25,mass_min_g=.1,mass_max_g=1e5,
  damuth_exponent=-.75,damuth_ref_mass_g=100,damuth_ref_N_per_species=500,
  max_range=2,n_quadrats=30,quadrat_area=4/90,disturb_frac=.5, # 33.333% coverage
  meteor_decay_fraction=.20,mother_points=2,
  example_size_classes=c(1,9,17,20,25),example_iteration=1,
  max_workers=200,
  max_parallel_n=1e10,
  sampling_max_attempts=150000
)

# max_parallel_n caps concurrent simulated individuals, not memory in bytes.
# Sampling may fail if 30 occupied, non-overlapping quadrats cannot be placed.

cfg$sampling_validator_version <- "occupied_random_direct_counts_v2"
cfg$analysis_quadrat_rule <- "exactly_30_preoccupied_fixed_locations_v1"

# Quadrat coordinates are lower-left corners; verify boundaries and non-overlap.
validate_layout <- function(xy) {
  xy <- as.data.frame(xy)
  names(xy)[1:2] <- c("x", "y")
  if (nrow(xy) != cfg$n_quadrats ||
      any(!is.finite(xy$x)) || any(!is.finite(xy$y)))
    stop("Invalid quadrat count or non-finite coordinates")

  side <- sqrt(cfg$quadrat_area)
  tol <- 1e-10
  if (any(xy$x < -tol | xy$y < -tol |
          xy$x + side > cfg$max_range + tol |
          xy$y + side > cfg$max_range + tol))
    stop("Quadrat outside the simulated landscape")

  dx <- abs(outer(xy$x, xy$x, "-"))
  dy <- abs(outer(xy$y, xy$y, "-"))
  overlap_pairs <- sum(upper.tri(dx) &
                       dx < side - tol & dy < side - tol)
  if (overlap_pairs > 0)
    stop("Overlapping quadrats detected: ", overlap_pairs, " pairs")
  invisible(as.integer(overlap_pairs))
}

out_dir <- Sys.getenv("THANOS_OUT",unset="output_30q_30_occupied_simple")
data_dir <- file.path(out_dir,"data")
fig_dir <- file.path(out_dir,"figures")
checkpoint_dir <- file.path(out_dir,paste0("checkpoints_n",cfg$n_iter))
config_path <- file.path(data_dir, "config.R")
if (file.exists(config_path)) {
  old_cfg <- dget(config_path)
  fixed_keys <- setdiff(names(cfg), c("max_workers", "max_parallel_n"))
  if (!isTRUE(all.equal(old_cfg[fixed_keys], cfg[fixed_keys],
                        check.attributes=FALSE)))
    stop("Output directory has a different sampling/model configuration: ", out_dir)
}
if (dir.exists(checkpoint_dir) &&
    length(list.files(checkpoint_dir, pattern="^class_[0-9][0-9]\\.rds$")) > 0 &&
    !file.exists(config_path))
  stop("Found checkpoints without a validated configuration in: ", out_dir)
dir.create(data_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(fig_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(checkpoint_dir,showWarnings=FALSE,recursive=TRUE)
dput(cfg,file=config_path)

# ---- 2. BODY-SIZE PARAMETERIZATION -------------------------------------------
anchors <- data.frame(
  body_mass_g=c(.1,100,1e5),
  s_pool=c(1000,100,30),
  sigma=c(.020,.080,.200)
)

mass <- 10^seq(log10(cfg$mass_min_g),log10(cfg$mass_max_g),length.out=cfg$n_size_classes)
interp_log <- function(x,y)
  10^splinefun(log10(anchors$body_mass_g),log10(y),method="monoH.FC")(log10(x))

regional_area <- cfg$max_range^2
s_pool <- as.integer(round(interp_log(mass,anchors$s_pool)))
sigma <- interp_log(mass,anchors$sigma)
mean_N_species <- cfg$damuth_ref_N_per_species*(mass/cfg$damuth_ref_mass_g)^cfg$damuth_exponent
n_sim_theoretical <- mean_N_species*s_pool

size_params <- data.frame(
  size_class=seq_along(mass),body_mass_g=mass,log10_body_mass=log10(mass),s_pool=s_pool,
  density_per_species=mean_N_species/regional_area,mean_regional_N_per_species=mean_N_species,
  n_sim_theoretical=n_sim_theoretical,n_sim=as.integer(round(n_sim_theoretical)),
  sigma=sigma,mother_points=cfg$mother_points
)
write.csv(size_params,file.path(data_dir,"size_class_parameters.csv"),row.names=FALSE)

# One sub-regional event center per iteration, shared among body-size classes.
set.seed(cfg$seed)
impact_table <- data.frame(
  iteration=seq_len(cfg$n_iter),
  impact_x=runif(cfg$n_iter,0,cfg$max_range),
  impact_y=runif(cfg$n_iter,0,cfg$max_range)
)
write.csv(impact_table,file.path(data_dir,"paired_impact_coordinates.csv"),row.names=FALSE)

# ---- 3. HELPERS ---------------------------------------------------------------
snap_n <- function(n,frac){x <- n*frac; floor(x)+rbinom(1,1,x-floor(x))}

remove_counts <- function(counts,n_remove){
  if(n_remove==0) return(counts)
  total <- sum(counts)
  idx <- findInterval(sample.int(total,n_remove,replace=FALSE)-1,c(0,cumsum(counts)))
  counts-tabulate(idx,nbins=length(counts))
}

mean_jaccard <- function(x){
  if(nrow(x)<2L || ncol(x)==0L) return(NA_real_)
  value <- mean(suppressWarnings(vegan::vegdist(x,method="jaccard",binary=TRUE)),na.rm=TRUE)
  if(is.nan(value)) NA_real_ else value
}

tatsumi_net <- function(pre,post){
  if(nrow(pre)<2L || sum(pre)==0 || sum(post)==0) return(NA_real_)
  ep <- ecopart::ecopart.multi(1*(pre>0),1*(post>0),index="whittaker",components="four")
  sum(ep[c("Extinction homogenization","Extinction differentiation")])
}

# Reject empty or overlapping random candidates until 30 are accepted.
sample_community <- function(sim){
  if(nrow(sim$census)<cfg$n_quadrats)
    stop("Fewer individuals than the 30 occupied quadrats required")
  side <- sqrt(cfg$quadrat_area)
  xy <- data.frame(x=numeric(cfg$n_quadrats),y=numeric(cfg$n_quadrats))
  counts <- vector("list",cfg$n_quadrats)
  accepted <- 0L
  candidates_drawn <- 0L
  rejected_overlap <- 0L
  rejected_empty <- 0L

  while(accepted < cfg$n_quadrats && candidates_drawn < cfg$sampling_max_attempts){
    candidates_drawn <- candidates_drawn+1L
    x <- runif(1L,0,cfg$max_range-side)
    y <- runif(1L,0,cfg$max_range-side)

    # Square-square overlap occurs only if both coordinate intervals overlap.
    if(accepted>0L && any(abs(xy$x[seq_len(accepted)]-x)<side &
                           abs(xy$y[seq_len(accepted)]-y)<side)){
      rejected_overlap <- rejected_overlap+1L
      next
    }

    # Exact local abundance from the pre-disturbance community.
    observed <- mobsim::abund_rect(x,y,side,side,sim)
    if(sum(observed)==0){
      rejected_empty <- rejected_empty+1L
      next
    }

    accepted <- accepted+1L
    xy$x[accepted] <- x
    xy$y[accepted] <- y
    counts[[accepted]] <- observed
  }

  if(accepted!=cfg$n_quadrats)
    stop("Only ",accepted," of ",cfg$n_quadrats,
         " occupied, non-overlapping quadrats found after ",candidates_drawn,
         " random draws. Increase sampling_max_attempts, or examine sparsity.")

  overlap_pairs <- validate_layout(xy)
  species <- unique(unlist(lapply(counts,names),use.names=FALSE))
  comm <- matrix(0L,nrow=cfg$n_quadrats,ncol=length(species),
                 dimnames=list(NULL,species))
  for(i in seq_len(cfg$n_quadrats))
    comm[i,names(counts[[i]])] <- as.integer(counts[[i]])
  comm <- comm[,colSums(comm)>0,drop=FALSE]
  if(ncol(comm)==0L || any(rowSums(comm)==0))
    stop("Sampling validation failed: an accepted quadrat was empty")

  rownames(comm) <- rownames(xy) <- paste0("Q",seq_len(cfg$n_quadrats))
  list(comm=comm,xy=xy,overlap_pairs=overlap_pairs,
       candidates_drawn=candidates_drawn,
       candidates_rejected_empty=rejected_empty,
       candidates_rejected_overlap=rejected_overlap)
}

# ---- 4. MORTALITY SCENARIOS --------------------------------------------------
regional_fixed <- function(comm,frac){
  out <- comm
  for(j in seq_len(ncol(out))){
    N <- sum(out[,j])
    if(N) out[,j] <- remove_counts(out[,j],snap_n(N,frac))
  }
  out
}

local_fixed <- function(comm,frac){
  out <- comm
  for(i in seq_len(nrow(out))) for(j in seq_len(ncol(out))){
    N <- out[i,j]
    if(N) out[i,j] <- N-snap_n(N,frac)
  }
  out
}

individual_mortality <- function(comm,frac)
  matrix(rbinom(length(comm),as.integer(comm),1-frac),
         nrow=nrow(comm),ncol=ncol(comm),dimnames=dimnames(comm))

allocate_gradient <- function(site_N,risk,target){
  if(target==0) return(rep(0L,length(site_N)))

  risk <- risk/max(risk)
  f <- function(z) sum(site_N*pmin(1,z*risk))

  upper <- 1
  while(f(upper)<target) upper <- upper*2

  # High precision is important here because site_N can contain millions of
  # individuals. The default uniroot tolerance can leave an integer rounding
  # remainder larger than the number of quadrats available for rounding.
  lambda <- uniroot(
    function(z) f(z)-target,
    c(0,upper),
    tol=1e-12
  )$root

  q <- site_N*pmin(1,lambda*risk)
  rem <- floor(q)
  left <- as.integer(target-sum(rem))

  if(left>0){
    residual <- q-rem
    ok <- which(residual>0 & rem<site_N)

    if(left>length(ok))
      stop("Sub-regional allocation rounding mismatch: left=",left,
           ", available quadrats=",length(ok))

    add <- ok[sample.int(
      length(ok),left,replace=FALSE,prob=residual[ok]
    )]
    rem[add] <- rem[add]+1L
  }

  as.integer(rem)
}

subregional_gradient <- function(comm,xy,cx,cy,frac){
  site_N <- rowSums(comm); target <- snap_n(sum(site_N),frac)
  d <- sqrt((xy$x-cx)^2+(xy$y-cy)^2)
  scale <- cfg$meteor_decay_fraction*cfg$max_range
  risk <- exp(-.5*(d/scale)^2)
  rem <- allocate_gradient(site_N,risk,target); out <- comm
  for(i in seq_len(nrow(out))) if(rem[i]) out[i,] <- remove_counts(out[i,],rem[i])
  out
}

# ---- 5. METRICS ---------------------------------------------------------------
# Every sampled quadrat is initially occupied; retain all 30 after mortality.
state_metrics <- function(x){
  a <- mean(specnumber(x)); g <- sum(colSums(x)>0)
  list(alpha=a,gamma=g,beta=if(a>0) g/a else NA_real_,jaccard=mean_jaccard(x))
}

baseline_metrics <- function(pre,s_pool_full){
  s <- state_metrics(pre)
  pre_empty <- rowSums(pre)==0
  data.frame(
    sampled_individuals=sum(pre),species_detection_fraction=ncol(pre)/s_pool_full,
    median_sampled_species_abundance=median(colSums(pre)),
    empty_quadrat_fraction=mean(pre_empty),
    n_quadrats_analysed=nrow(pre),
    n_quadrats_excluded_pre_empty=sum(pre_empty),
    mean_local_richness_pre=s$alpha,regional_richness_pre=s$gamma,
    whittaker_beta_pre=s$beta,mean_jaccard_pre=s$jaccard
  )
}

response_metrics <- function(pre,post){
  stopifnot(identical(dim(pre),dim(post)))
  s0 <- state_metrics(pre); s1 <- state_metrics(post)
  data.frame(
    mean_local_richness_post=s1$alpha,regional_richness_post=s1$gamma,
    whittaker_beta_post=s1$beta,mean_jaccard_post=s1$jaccard,
    n_quadrats_analysed=nrow(pre),
    empty_quadrat_fraction_post=mean(rowSums(post)==0),
    mortality_fraction_realized=1-sum(post)/sum(pre),
    alpha_change_pct=100*(s1$alpha-s0$alpha)/s0$alpha,
    gamma_change_pct=100*(s1$gamma-s0$gamma)/s0$gamma,
    tatsumi_extirpation_net=tatsumi_net(pre,post)
  )
}

# ---- 6. ONE CLASS x ITERATION ------------------------------------------------
run_one <- function(pars,iteration,impact,return_example=FALSE){
  run_id <- (pars$size_class-1)*cfg$n_iter+iteration
  seed <- cfg$seed+5*run_id

  set.seed(seed)
  sim <- sim_thomas_community(
    s_pool=pars$s_pool,n_sim=pars$n_sim,sad_type="lnorm",fix_s_sim=TRUE,
    sad_coef=list(cv_abund=1),sigma=pars$sigma,mother_points=cfg$mother_points,
    xrange=c(0,cfg$max_range),yrange=c(0,cfg$max_range)
  )
  s <- sample_community(sim); pre <- s$comm; xy <- s$xy
  stopifnot(nrow(pre)==cfg$n_quadrats,all(rowSums(pre)>0))

  scenario_seed <- setNames(
    seed+1:4,
    c("regional_fixed_50","local_fixed_50","probabilistic_individual_50","subregional_50")
  )
  set.seed(scenario_seed[1]); regional <- regional_fixed(pre,cfg$disturb_frac)
  set.seed(scenario_seed[2]); local <- local_fixed(pre,cfg$disturb_frac)
  set.seed(scenario_seed[3]); individual <- individual_mortality(pre,cfg$disturb_frac)
  set.seed(scenario_seed[4]); subregional <- subregional_gradient(
    pre,xy,impact$impact_x,impact$impact_y,cfg$disturb_frac
  )

  disturbances <- list(
    regional_fixed_50=regional,
    local_fixed_50=local,
    probabilistic_individual_50=individual,
    subregional_50=subregional
  )

  common <- data.frame(
    size_class=pars$size_class,body_mass_g=pars$body_mass_g,
    log10_body_mass=pars$log10_body_mass,iteration=iteration,
    run_id=run_id,community_seed=seed
  )

  res <- bind_rows(lapply(names(disturbances),function(z)
    cbind(common,scenario=z,scenario_seed=unname(scenario_seed[z]),
          response_metrics(pre,disturbances[[z]]))))
  base <- cbind(common,baseline_metrics(pre,pars$s_pool),
                n_quadrats_sampled=nrow(pre),quadrat_overlap_pairs=s$overlap_pairs,
                quadrat_candidates_drawn=s$candidates_drawn,
                quadrat_rejected_empty=s$candidates_rejected_empty,
                quadrat_rejected_overlap=s$candidates_rejected_overlap)

  ex <- NULL
  if(return_example){
    states <- c(list(pre_disturbance=pre),disturbances); preN <- rowSums(pre)
    ex <- bind_rows(lapply(names(states),function(z){
      postN <- rowSums(states[[z]])
      data.frame(
        size_class=pars$size_class,body_mass_g=pars$body_mass_g,
        x=xy$x,y=xy$y,quadrat=rownames(pre),scenario=z,
        pre_N=preN,post_N=postN,
        mortality_fraction=1-postN/preN,
        occupied_pre=preN>0
      )
    }))
  }

  rm(sim,regional,local,individual,subregional,disturbances)
  gc(verbose=FALSE)
  list(results=res,baseline=base,example=ex)
}

# ---- 7. PARALLEL RUN + CHECKPOINTS -------------------------------------------
available_cores <- parallel::detectCores()
max_workers <- min(cfg$max_workers,cfg$n_iter,max(1L,available_cores-2L))
workers_for <- function(n_sim)
  max(1L,min(max_workers,as.integer(floor(cfg$max_parallel_n/n_sim))))

class_runs <- vector("list",nrow(size_params))

for(c in seq_len(nrow(size_params))){
  pars <- size_params[c,,drop=FALSE]
  checkpoint <- file.path(checkpoint_dir,sprintf("class_%02d.rds",c))

  if(file.exists(checkpoint)){
    cat("Class",c,"already complete; loading checkpoint\n")
    class_runs[[c]] <- readRDS(checkpoint)
    if(!all(class_runs[[c]]$results$scenario %in%
            c("regional_fixed_50","local_fixed_50","probabilistic_individual_50","subregional_50")) ||
       !all(c("n_quadrats_analysed","empty_quadrat_fraction_post") %in%
            names(class_runs[[c]]$results)) ||
       !all(c("quadrat_candidates_drawn","quadrat_rejected_empty",
                "quadrat_rejected_overlap") %in% names(class_runs[[c]]$baseline)) ||
       !all(class_runs[[c]]$baseline$n_quadrats_analysed==cfg$n_quadrats))
      stop("Checkpoint predates the 30-occupied-quadrat sampling rule: ",checkpoint)
    next
  }

  n_workers <- workers_for(pars$n_sim)
  cat("Class",c,"/",nrow(size_params),"|",signif(pars$body_mass_g,4),"g | n_sim",
      format(pars$n_sim,big.mark=",",scientific=FALSE),"| workers",n_workers,"\n")

  ans <- parallel::mclapply(
    seq_len(cfg$n_iter),
    function(i){
      make_ex <- c %in% cfg$example_size_classes && i==cfg$example_iteration
      tryCatch(
        list(
          ok=TRUE,
          iteration=i,
          value=run_one(pars,i,impact_table[i,,drop=FALSE],make_ex)
        ),
        error=function(e) list(
          ok=FALSE,
          iteration=i,
          error=conditionMessage(e),
          call=paste(deparse(conditionCall(e)),collapse=" ")
        )
      )
    },
    mc.cores=n_workers,
    mc.preschedule=FALSE,
    mc.set.seed=FALSE
  )

  failed <- !vapply(ans,function(x) isTRUE(x$ok),logical(1))
  if(any(failed)){
    cat("\nFailed iterations in class",c,":\n")
    for(z in ans[failed]){
      cat(" iteration",z$iteration,":",z$error,"\n",
          " call:",z$call,"\n")
    }
    stop("Class ",c," had ",sum(failed)," failed iteration(s).")
  }

  ans <- lapply(ans,`[[`,"value")

  class_runs[[c]] <- list(
    results=bind_rows(lapply(ans,`[[`,"results")),
    baseline=bind_rows(lapply(ans,`[[`,"baseline")),
    example=bind_rows(lapply(ans,`[[`,"example"))
  )
  saveRDS(class_runs[[c]],checkpoint,compress=FALSE)
  rm(ans); gc(verbose=FALSE)
}

results <- bind_rows(lapply(class_runs,`[[`,"results"))
baseline <- bind_rows(lapply(class_runs,`[[`,"baseline"))
example_data <- bind_rows(lapply(class_runs,`[[`,"example"))

# Auditable guarantee covering every completed iteration.
stopifnot(nrow(baseline) == cfg$n_size_classes * cfg$n_iter,
          all(baseline$n_quadrats_sampled == cfg$n_quadrats),
          all(baseline$quadrat_overlap_pairs == 0L),
          all(baseline$n_quadrats_analysed == cfg$n_quadrats),
          all(baseline$n_quadrats_excluded_pre_empty == 0L),
          nrow(results) == 4L*nrow(baseline),
          all(results$n_quadrats_analysed == cfg$n_quadrats))

main_metrics <- c("mortality_fraction_realized","alpha_change_pct","gamma_change_pct","tatsumi_extirpation_net")
summarise_95 <- function(dat,groups,metrics) dat %>%
  group_by(across(all_of(groups))) %>%
  summarise(across(all_of(metrics),list(
    mean=~mean(.x,na.rm=TRUE),
    q025=~quantile(.x,.025,na.rm=TRUE,names=FALSE),
    q975=~quantile(.x,.975,na.rm=TRUE,names=FALSE)),
    .names="{.col}_{.fn}"),.groups="drop")

summary_results <- summarise_95(
  results,c("size_class","body_mass_g","log10_body_mass","scenario"),main_metrics
)
baseline_plot_metrics <- c(
  "sampled_individuals","species_detection_fraction",
  "median_sampled_species_abundance","empty_quadrat_fraction"
)
baseline_summary <- summarise_95(
  baseline,c("size_class","body_mass_g","log10_body_mass"),baseline_plot_metrics
)

write.csv(results,file.path(data_dir,"iteration_results.csv"),row.names=FALSE)
write.csv(summary_results,file.path(data_dir,"scenario_summary.csv"),row.names=FALSE)
write.csv(baseline,file.path(data_dir,"baseline_iteration_results.csv"),row.names=FALSE)
write.csv(baseline_summary,file.path(data_dir,"baseline_summary.csv"),row.names=FALSE)
if(nrow(example_data)) write.csv(example_data,file.path(data_dir,"example_scenario_map_data.csv"),row.names=FALSE)

run_summary <- c(
  paste("Size classes:",nrow(size_params)),
  paste("Iterations:",cfg$n_iter),
  paste("Detected cores:",available_cores),
  paste("Maximum workers:",max_workers),
  paste("Parallel n_sim cap:",format(cfg$max_parallel_n,scientific=FALSE)),
  paste("Mother points/species:",cfg$mother_points),
  paste("Quadrats per iteration:", cfg$n_quadrats),
  paste("Analysis quadrat rule:",cfg$analysis_quadrat_rule),
  paste("Occupied quadrats in every iteration:",min(baseline$n_quadrats_analysed)),
  paste("Mean candidate quadrats drawn:",round(mean(baseline$quadrat_candidates_drawn),1)),
  paste("Max candidate quadrats drawn:",max(baseline$quadrat_candidates_drawn)),
  paste("Sampled fraction (%):",100*cfg$n_quadrats*cfg$quadrat_area/regional_area),
  paste("Layouts validated:", nrow(baseline)),
  paste("Max overlapping quadrat pairs:", max(baseline$quadrat_overlap_pairs)),
  paste("Largest-class expected N/species:",signif(tail(size_params$mean_regional_N_per_species,1),4)),
  paste("Largest-class total N:",tail(size_params$n_sim,1))
)
writeLines(run_summary,file.path(data_dir,"run_summary.txt"))
cat("\n",paste(run_summary,collapse="\n"),"\n",sep="")

# ---- 8. PLOT DATA + PUBLICATION STYLE ---------------------------------------
scenario_order <- c("regional_fixed_50","local_fixed_50","probabilistic_individual_50","subregional_50")
scenario_labels <- c(
  regional_fixed_50="Regional fixed 50%",
  local_fixed_50="Local fixed 50%",
  probabilistic_individual_50="Probabilistic individual 50%",
  subregional_50="Sub-regional event"
)
scenario_colours <- c(
  regional_fixed_50="#0072B2",local_fixed_50="#009E73",
  probabilistic_individual_50="#D55E00",subregional_50="#CC79A7"
)
scenario_lty <- c(
  "Regional fixed 50%"="solid","Local fixed 50%"="dashed",
  "Probabilistic individual 50%"="dotdash","Sub-regional event"="longdash"
)
plot_cols <- setNames(unname(scenario_colours[scenario_order]),unname(scenario_labels[scenario_order]))
state_cols <- c("Pre-disturbance"="#000000",plot_cols)
state_lty <- c("Pre-disturbance"="solid",scenario_lty)

format_mass <- function(g)
  if(g>=1000) paste0(format(signif(g/1000,3),trim=TRUE)," kg") else paste0(format(signif(g,3),trim=TRUE)," g")

mass_axis <- function() scale_x_continuous(
  breaks=-1:5,
  labels=c("0.1 g","1 g","10 g","100 g","1 kg","10 kg","100 kg"),
  expand=expansion(mult=c(.015,.02))
)

theme_pub <- function() theme_bw(base_size=10.5)+theme(
  panel.grid.minor=element_blank(),
  panel.grid.major=element_blank(),
  strip.background=element_rect(fill="grey92",colour="grey45",linewidth=.4),
  legend.position="bottom",legend.key.width=grid::unit(1.25,"cm"),
  legend.key.height=grid::unit(.45,"cm"),plot.margin=margin(6,8,6,6)
)

save_pub <- function(p,name,w,h){
  dev <- if(capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(fig_dir,paste0(name,".pdf")),p,width=w,height=h,units="in",device=dev,bg="white")
  ggsave(file.path(fig_dir,paste0(name,".png")),p,width=w,height=h,units="in",dpi=600,bg="white")
}

param_labels <- c(
  s_pool="Species pool",mean_regional_N_per_species="Mean regional N / species",
  n_sim="Total regional abundance",sigma="Spatial spread (sigma)"
)
param_df <- size_params %>%
  select(log10_body_mass,all_of(names(param_labels))) %>%
  pivot_longer(-log10_body_mass,names_to="parameter",values_to="value") %>%
  mutate(parameter=factor(unname(param_labels[parameter]),levels=unname(param_labels)))
anchor_df <- bind_rows(
  anchors %>% transmute(log10_body_mass=log10(body_mass_g),parameter="Species pool",value=s_pool),
  anchors %>% transmute(log10_body_mass=log10(body_mass_g),parameter="Spatial spread (sigma)",value=sigma)
) %>% mutate(parameter=factor(parameter,levels=unname(param_labels)))

base_labels <- c(
  sampled_individuals="Sampled individuals",
  species_detection_fraction="Species-pool detection fraction",
  median_sampled_species_abundance="Median sampled N / species",
  empty_quadrat_fraction="Empty before disturbance (always zero by design)"
)
base_df <- baseline_summary %>%
  pivot_longer(-c(size_class,body_mass_g,log10_body_mass),
               names_to=c("metric","stat"),names_pattern="(.+)_(mean|q025|q975)$",values_to="value") %>%
  pivot_wider(names_from=stat,values_from=value) %>%
  mutate(metric=factor(unname(base_labels[metric]),levels=unname(base_labels)))

metric_labels <- c(
  mortality_fraction_realized="Realized mortality (%)",
  alpha_change_pct="Local richness change (%)",
  gamma_change_pct="Regional richness change (%)",
  tatsumi_extirpation_net="Net Whittaker beta change (Tatsumi)"
)
response_df <- summary_results %>%
  pivot_longer(-c(size_class,body_mass_g,log10_body_mass,scenario),
               names_to=c("metric","stat"),names_pattern="(.+)_(mean|q025|q975)$",values_to="value") %>%
  pivot_wider(names_from=stat,values_from=value) %>%
  mutate(
    mean=ifelse(metric=="mortality_fraction_realized",mean*100,mean),
    q025=ifelse(metric=="mortality_fraction_realized",q025*100,q025),
    q975=ifelse(metric=="mortality_fraction_realized",q975*100,q975),
    metric=factor(unname(metric_labels[metric]),levels=unname(metric_labels[main_metrics])),
    scenario=factor(scenario,levels=scenario_order,labels=scenario_labels[scenario_order])
  )

# ---- 9. FIGURES --------------------------------------------------------------
p_parameters <- ggplot(param_df,aes(log10_body_mass,value))+
  geom_line(linewidth=.7)+geom_point(size=1)+
  geom_point(data=anchor_df,shape=21,fill="white",size=2.8,stroke=.9)+
  scale_y_log10()+mass_axis()+facet_wrap(~parameter,scales="free_y",ncol=2)+theme_pub()+
  labs(x="Representative body mass",y=NULL,
       caption=paste0("Mother points fixed at ",cfg$mother_points,
                      "; open circles mark species-pool and sigma anchors."))

p_baseline <- ggplot(base_df,aes(log10_body_mass,mean))+
  geom_ribbon(aes(ymin=q025,ymax=q975),alpha=.16,colour=NA)+
  geom_line(linewidth=.75)+geom_point(size=.9)+mass_axis()+
  facet_wrap(~metric,scales="free_y",ncol=2)+theme_pub()+
  labs(x="Representative body mass",y=NULL)

p_response <- ggplot(response_df,aes(log10_body_mass,mean,group=scenario,colour=scenario,linetype=scenario))+
  geom_hline(yintercept=0,linewidth=.25,colour="grey55")+
  geom_ribbon(aes(ymin=q025,ymax=q975,fill=scenario),alpha=.16,colour=NA,show.legend=FALSE,na.rm=TRUE)+
  geom_line(linewidth=.8,na.rm=TRUE)+geom_point(size=.85,na.rm=TRUE)+mass_axis()+
  scale_colour_manual(values=plot_cols,drop=FALSE)+scale_fill_manual(values=plot_cols,drop=FALSE)+
  scale_linetype_manual(values=scenario_lty,drop=FALSE)+
  facet_wrap(~metric,scales="free_y",ncol=2)+theme_bw()+
  theme(legend.position="bottom",panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))+
  labs(x="Representative body mass",y=NULL,colour="Scenario",linetype="Scenario",fill="Scenario")+
  scale_y_continuous(expand=expansion(mult=c(.03,.15)))

p_response_tagged <- egg::tag_facet(
  p_response,open="",close="",
  tag_pool=c("A) Realized mortality (%)","B) Local richness change (%)",
             "C) Regional richness change (%)","D) Net Whittaker beta-diversity change"),
  hjust=-0.05,size=3.5
)

raw_pre <- baseline %>%
  transmute(size_class,body_mass_g,log10_body_mass,iteration,series="Pre-disturbance",
            `Mean Jaccard dissimilarity`=mean_jaccard_pre,`Whittaker beta`=whittaker_beta_pre,
            `Mean local richness`=mean_local_richness_pre,`Regional richness`=regional_richness_pre) %>%
  pivot_longer(c(`Mean Jaccard dissimilarity`,`Whittaker beta`,`Mean local richness`,`Regional richness`),
               names_to="metric",values_to="value")
raw_post <- results %>%
  transmute(size_class,body_mass_g,log10_body_mass,iteration,series=unname(scenario_labels[scenario]),
            `Mean Jaccard dissimilarity`=mean_jaccard_post,`Whittaker beta`=whittaker_beta_post,
            `Mean local richness`=mean_local_richness_post,`Regional richness`=regional_richness_post) %>%
  pivot_longer(c(`Mean Jaccard dissimilarity`,`Whittaker beta`,`Mean local richness`,`Regional richness`),
               names_to="metric",values_to="value")
raw_diversity <- bind_rows(raw_pre,raw_post) %>%
  group_by(size_class,body_mass_g,log10_body_mass,series,metric) %>%
  summarise(mean=mean(value,na.rm=TRUE),q025=quantile(value,.025,na.rm=TRUE,names=FALSE),
            q975=quantile(value,.975,na.rm=TRUE,names=FALSE),.groups="drop")
raw_order <- c("Pre-disturbance",unname(scenario_labels[scenario_order]))
raw_diversity$series <- factor(raw_diversity$series,levels=raw_order)
write.csv(raw_diversity,file.path(data_dir,"raw_diversity_summary.csv"),row.names=FALSE)

p_raw_diversity <- ggplot(raw_diversity,aes(log10_body_mass,mean,group=series,colour=series,linetype=series))+
  geom_ribbon(aes(ymin=q025,ymax=q975,fill=series),alpha=.12,colour=NA,show.legend=FALSE,na.rm=TRUE)+
  geom_line(linewidth=.75,na.rm=TRUE)+geom_point(size=.8,na.rm=TRUE)+mass_axis()+
  scale_colour_manual(values=state_cols,drop=FALSE)+scale_fill_manual(values=state_cols,drop=FALSE)+
  scale_linetype_manual(values=state_lty,drop=FALSE)+facet_wrap(~metric,scales="free_y",ncol=2)+theme_pub()+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))+
  labs(x="Representative body mass",y=NULL,colour="State / scenario",linetype="State / scenario",fill="State / scenario")

p_examples <- list(); example_order <- c("pre_disturbance",scenario_order)
example_labels <- c(pre_disturbance="Pre-disturbance",scenario_labels)
if(nrow(example_data)){
  example_data$scenario <- factor(example_data$scenario,levels=example_order,labels=example_labels[example_order])
  for(sc_class in cfg$example_size_classes){
    ex <- example_data %>% filter(size_class==sc_class)
    if(!nrow(ex)) next
    event_center <- data.frame(
      x=impact_table$impact_x[cfg$example_iteration],y=impact_table$impact_y[cfg$example_iteration],
      scenario=factor(example_labels["subregional_50"],levels=example_labels[example_order])
    )
    p_examples[[as.character(sc_class)]] <- ggplot(ex,aes(x,y))+
      geom_point(aes(colour=mortality_fraction),size=3.5,alpha=.95)+
      geom_point(data=event_center,aes(x,y),inherit.aes=FALSE,shape=8,size=4,colour="black",stroke=1)+
      scale_colour_gradientn(colours=c("#0072B2","#F0E442","#D55E00"),limits=c(0,1),
                             breaks=c(0,.25,.5,.75,1),name="Local mortality\nfraction")+
      coord_equal(xlim=c(0,cfg$max_range),ylim=c(0,cfg$max_range),expand=FALSE)+
      facet_wrap(~scenario,ncol=3)+theme_pub()+theme(legend.position="right")+
      labs(x="Simulated x coordinate",y="Simulated y coordinate",
           title=paste0("Class ",sc_class," (",format_mass(unique(ex$body_mass_g)),")"))
  }
}

# ---- 10. SAVE + REPRODUCIBILITY ---------------------------------------------
save_pub(p_parameters,"figure_parameters",7.2,6.2)
save_pub(p_baseline,"figure_baseline_diagnostics",7.2,6.2)
save_pub(p_response_tagged,"figure_main_response",7.2,6.5)
save_pub(p_raw_diversity,"figure_supplement_raw_diversity",7.2,6.5)
if(length(p_examples)) for(nm in names(p_examples))
  save_pub(p_examples[[nm]],paste0("figure_spatial_example_class_",nm),8.2,5.7)

capture.output(sessionInfo(),file=file.path(data_dir,"sessionInfo.txt"))
cat("\nOutputs saved in:",normalizePath(out_dir,mustWork=FALSE),"\n")
