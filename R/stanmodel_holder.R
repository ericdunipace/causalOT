
# names of stan models
cot_stanfiles <- c("barycenter_projection", "gp_hyper_dose", "gp_hyper")

cot_stanfiles <- sapply(cot_stanfiles, function(sm){
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(sm, ".stan"))
  
  return(readLines(stan_file))
})

stanbuilder <- function(model_name) {
  # if (is.null(cot_stanmodels[[model_name]])) {
    
    stanfit <- rstan::stan_model(model_code = cot_stanfiles[[model_name]],
                                 model_name = model_name,
                                    allow_undefined = TRUE,
                                    obfuscate_model_name = FALSE)
    return(stanfit)
    # create stanmodel object
  #   cot_stanmodels[[model_name]] <<- stanfit
  # }
  # return(cot_stanmodels[[model_name]])
}
