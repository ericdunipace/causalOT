# Bounded replay of the saved precision diagnostic. No fits or simulation grid.
# Rscript dev/polish_precision_followup_validate.R [checkpoint-root] [new-output-dir]
args <- commandArgs(TRUE)
root <- if(length(args)) args[[1]] else '/Users/eifer/GoogleDrive/Research/Jose Zubizarreta Research/causalOT'
out <- if(length(args)>1L) args[[2]] else 'dev/polish_precision_followup_results'
if(dir.exists(out)) stop('Use a new output directory; do not mix corrected and old results')
dir.create(out,recursive=TRUE)
suppressMessages(pkgload::load_all('.',quiet=TRUE))
torch::torch_set_num_threads(1L)
phase <- readRDS(file.path(root,'dev/newton_warm_evaluation_results/small/cells',
  'small_n1000_d3_l0.05_r1_cpu_float64_mirror_one/phase1.rds'))
d <- readRDS(file.path(root,'dev/newton_warm_diagnostic_results/diagnostic.rds'))
T <- causalOT:::cot_polish_tensor
x <- T(phase$x); y <- T(phase$y)
stopifnot(phase$cost_source=='solver_rounded')
C <- lapply(list(list(x,y),list(x,x),list(y,y)),function(z)
  causalOT:::cot_polish_cost(causalOT:::costTensor$new(z[[1]],z[[2]],p=2)))
names(C) <- c('xy','xx','yy')
for(nm in names(C)) C[[nm]]$diameter <- causalOT:::cot_polish_diameter(C[[nm]])
stopifnot(is.null(C$xy$dense),is.null(C$xx$dense),is.null(C$yy$dense))
control <- causalOT:::cot_polish_control(sinkhorn_niter=256L)
cold_control <- causalOT:::cot_polish_control(sinkhorn_niter=1000L)
cat('Cold corrected candidate evaluation\n');flush.console()
cold <- causalOT:::cot_polish_eval(d$candidate,phase$b,C$xy,C$xx,phase$lambda,cold_control)
y_self <- causalOT:::cot_polish_self(phase$b,C$yy,phase$lambda,cold_control)
residual <- function(weights,ev) {
  a <- T(weights);a <- a/a$sum();b <- T(phase$b);b <- b/b$sum();pot <- ev$pot
  row <- causalOT:::cot_polish_softmin(C$xy,pot$g_yx+b$log()*phase$lambda,phase$lambda)
  col <- causalOT:::cot_polish_softmin(C$xy,pot$f_xy+a$log()*phase$lambda,phase$lambda,TRUE)
  target <- causalOT:::cot_polish_softmin(C$xx,pot$f_xx+a$log()*phase$lambda,phase$lambda)
  c(marginal=((a*((pot$f_xy-row)/phase$lambda)$exp()-a)$abs()$sum()+
    (b*((pot$g_yx-col)/phase$lambda)$exp()-b)$abs()$sum())$item(),
    self=(target-pot$f_xx)$abs()$max()$item())
}
rows <- list(); evaluations <- list()
starts <- list(previous_new_weights=list(a=d$candidate,pot=lapply(d$incoming,T)),
  old_cold_new_weights=list(a=d$candidate,pot=lapply(d$cold$pot,T)),
  corrected_cold_new_weights=list(a=d$candidate,pot=cold$pot),
  previous_unchanged_weights=list(a=phase$a,pot=lapply(d$incoming,T)))
for(name in names(starts)) {
  cat('Warm corrected evaluation:',name,'\n');flush.console()
  start <- starts[[name]]; t <- proc.time()[['elapsed']]
  ev <- causalOT:::cot_polish_warm(start$a,phase$b,C$xy,C$xx,phase$lambda,start$pot,control)
  r <- residual(start$a,ev)
  paired <- name!='previous_unchanged_weights'
  loss_difference <- if(paired) abs(ev$loss-cold$loss) else NA_real_
  grad_difference <- if(paired) max(abs(ev$grad-cold$grad)) else NA_real_
  stopifnot(r[['marginal']]<=1.2e-13,r[['self']]<=1.2e-13,ev$iterations<256L)
  if(paired) stopifnot(loss_difference<=1e-12,grad_difference<=1e-11)
  rows[[name]] <- data.frame(case=name,iterations=ev$iterations,time=proc.time()[['elapsed']]-t,
    loss=ev$loss,marginal=r[['marginal']],self=r[['self']],
    cold_loss_difference=loss_difference,cold_gradient_difference=grad_difference)
  evaluations[[name]] <- list(loss=ev$loss,grad=ev$grad,pot=lapply(ev$pot,as.numeric))
  print(rows[[name]])
  write.csv(do.call(rbind,rows),file.path(out,'summary.csv'),row.names=FALSE)
}
saveRDS(list(rows=rows,cold=list(loss=cold$loss,grad=cold$grad,pot=lapply(cold$pot,as.numeric)),
  y_self=y_self,evaluations=evaluations),file.path(out,'corrected.rds'))
cat('Saved bounded corrected results:',out,'\n')
