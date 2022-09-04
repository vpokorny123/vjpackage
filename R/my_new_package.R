make_demo_table<-function(data_frame){
for (j in colnames(data_frame)) {
  if (j == "grp" | j == "group" | j == "Group") {
    group_var_name = j
    freqs = table(data_frame[group_var_name])
    final_table = t(as.data.frame(paste0(names(freqs)," n=(",unname(freqs),")")))
    final_table = unname(final_table)
    final_table = cbind('Variable',final_table,'Statistic','Post Hoc Contrasts')
    }

  else if (j == "Gender" | j == "gender" | j == "gndr" | j== "Sex" | j == "sex") {
    print("calculating percent female: assuming lowest number is male and highest number is female (e.g. 0 = male, 1 = female). YOU NEED TO CHECK TO MAKE SURE THAT IS RIGHT. If it isn't right then recode so it is right lol")
    freqs = table(as.factor(data_frame[[group_var_name]]),as.factor(data_frame[[j]]))
    res = chisq.test(freqs)
    es = effectsize::effectsize(res)
    res_pretty = paste0(greekLetters::greeks("Chi^2"),"(",res$parameter,")=",round(res$statistic,digits = 3),
                        ', p=', round(res$p.value, digits = 3),
                        ', Cramer',"'",'s V=',round(es$Cramers_v, digits = 2))
    posthoc_res = chisq.posthoc.test::chisq.posthoc.test(freqs)
    #browser()
    posthoc_pvals = posthoc_res[posthoc_res["Value"] == 'p values',]
    if (any(posthoc_pvals["1"]<.05)){
      print("haven't coded this part yet")
    }
    else{
      posthoc_string = " "
    }
    perc = round(freqs[,2]/(freqs[,1]+freqs[,2])*100)
    final_row = unname(cbind("Percent Female",t(unname(perc)),res_pretty,posthoc_string))
    final_table = as.data.frame(rbind(final_table,final_row))
    rownames(final_table)<- NULL
  }
  else{
  means = aggregate(get(j) ~ get(group_var_name),data = data_frame, mean, drop = FALSE)
  means = round(means[,2],digits = 2)
  stds = aggregate(get(j) ~ get(group_var_name),data = data_frame, sd, drop = FALSE)
  stds = round(stds[,2],digits = 2)
  final_row = paste0(means," (",stds,")")
  aov_res = aov(get(j)~ get(group_var_name), data =data_frame)
  es = effectsize::effectsize(aov_res)
  #browser()
  aov_sum = unlist(summary(aov_res))
  pval = round(aov_sum["Pr(>F)1"],digits = 2)
  ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
  stats = paste0("F(",aov_sum[1], ",", aov_sum[2],")=",round(aov_sum["F value1"],digits = 2),
                 ", p", pval,
                 ", ",greekLetters::greeks("eta^2"),"=",round(es$Eta2,digits = 2))
  final_final_row = as.data.frame(unname(cbind(j,t(final_row),stats," ")))
  final_table = rbind(final_table,final_final_row)
}
}
return(final_table)
}

#######################################
#######################################
#######################################


pub_ready_stats<-function(x) {
  #this is for afex aov_ez
  #correction should either be none, hf or gg
  output = NULL
  if (class(x)[1]=="afex_aov") {
    anov_table<-x[["anova_table"]]
    pest <- afex::nice(x, 'pes')
    rownames(pest)<- pest$Effect
    for (j in rownames(anov_table)) {
    fstat = round(anov_table[j,"F"], digits = 2)
    pval = round(anov_table[j,"Pr(>F)"], digits = 2)
    ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
    df1 = round(anov_table[j,"num Df"],digits =2)
    df2 = round(anov_table[j,"den Df"],digits =2)
    p_eta <- pest[j,'pes']
    pub_ready = paste0('F(',df1,',',df2,')=',fstat,', p',pval,', ', greekLetters::greeks("eta^2"), '=', p_eta)
    pub_ready = unname(cbind(j,pub_ready))
    output = unname(rbind(output,pub_ready))
    }
  }
  if (all(class(x)==c("psych","fa"))){
    #x$chi is the emprically derived which is recommended when normal theory fails (e.g. non-positive definite matrix)
    #x$STATISTIC is what we want to report unless something funky's goin on
    chi<-round(x$STATISTIC,digits = 2)
    ifelse(round(x$PVAL,digits = 2)==0, pval<-"<.001",pval<-paste0('=',round(x$PVAL,digits = 2)))
    chi_rep<- paste0(greekLetters::greeks("chi^2"),'(',x$dof,')=',chi,', p',pval)
    x$CFI<-((x$null.chisq-x$null.dof)-
                    (x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
    other_fits<-paste0(chi_rep,', TLI=',round(x$TLI,digits=2),
                       ', CFI=', round(x$CFI,digits=2),
                       ', RMSEA=', round(x$RMSEA,digits=2)[1])
    output = other_fits
  }
  #browser()
  if ("method" %in% names(x)){
  if (x$method == "Pearson's product-moment correlation") {
    r = round(x$estimate,digits = 2)
    df = x$parameter
    ifelse(round(x$p.value,digits = 2)==0, pval<-"<.001",pval<-paste0('=',round(x$p.value,digits = 3)))
    output<-paste0("r(",df,")=",r,", p",pval)}
  if (grepl('t-test',x$method )) {
    #browser()
      t = round(x$statistic,digits = 2)
      df = round(x$parameter,digits = 2)
      ifelse(round(x$p.value,digits = 2)==0, pval<-"<.001",pval<-paste0('=',round(x$p.value,digits = 3)))
      output<-paste0("t(",df,")=",t,", p",pval)

  }
  }
  if (all(grepl("mediate",class(x)))) {
    library(greekLetters)
    acme_b<-round(x$d0,digits = 2)
    acme_p<-round(x$d0.p,digits = 2)
    acme_ci<-round(x$d0.ci,digits = 2)
    output<-paste0("ACME ", greeks("beta"),"=",acme_b,", 95% CI [",acme_ci[1],", ",
                     acme_ci[2],"], p=",acme_p)
  }
return(output)
}

#######################################
#######################################
#######################################

cor_conf_int<-function(raw_r,n,alpha){
  phishman_z = .5*log((1+raw_r)/(1-raw_r)) # log is ln in R; log10 is base 10 log in r
  thresh = qnorm(1-(alpha/2))
  se = 1/sqrt(n-3)
  upper = phishman_z + (thresh *se)
  lower = phishman_z - (thresh *se)
  backtrans_upper = (exp(2*upper) - 1) / (exp(2*upper) + 1)
  backtrans_lower = (exp(2*lower) - 1) / (exp(2*lower) + 1)
  conf_int = c(backtrans_lower,backtrans_upper)
  return(conf_int)
}

#######################################
#######################################
#######################################

vjpscatter<- function(x,y,x_label,y_label,SB,plot_title, save_title, alpha=NULL) {
  #manually listwise delete
  main_data = data.frame(x,y)
  main_data <- main_data[complete.cases(main_data), ]
  x = main_data[1]
  y = main_data[2]
  if (!hasArg(group)) {
  x = scale(x)
  y = scale(y)
  if (!hasArg(SB)){
    SB = FALSE
  }
  if (is.null(alpha)){
    alpha = .05
  }
  cor_res = round(cor(x, y,use="pairwise.complete.obs"),digits = 3)
  raw_conf_int = round(cor_conf_int(cor_res, length(x),alpha),digits = 3)
  #cal
  if (SB == TRUE) {
  SB_cor = round((2*cor_res)/ (1+cor_res),digits = 3)
  SB_conf_int = round(cor_conf_int(SB_cor, length(x),alpha),digits = 3)
  #browser()
  }
  scatter_df = data.frame(x,y)
  pretty_plot<-ggplot(scatter_df, aes(x,y)) +
    geom_point(alpha=.3) +
    theme_classic() +
    geom_smooth(method= "lm",color = "black",fill = "black",alpha = .3)+
    annotate(geom = "text",fontface = c('bold'),
             min(na.omit(x)),
             max(na.omit(y)),
             label = paste0("Raw r = ",cor_res,
                            " [",raw_conf_int[1],", ",raw_conf_int[2],"]",
                            "\nSB r = ",SB_cor,
                            " [",SB_conf_int[1],", ",SB_conf_int[2],"]"),
             hjust = 0, vjust = .83, size = 2.8) +
    xlab(x_label) +
    ylab(y_label) +
    theme(axis.title=element_text(size=10))

  if (hasArg(save_title)){
    ggsave(paste0(save_title,".png"),width = 5.5, height = 5)
    print(getwd())
  }
  }
  print(pretty_plot)
}

#######################################
#######################################
#######################################

vjplot_paired_ttest<- function(x, y, sub_col_name, x_lab, y_lab, cond_labs,
                               plot_title, main_data, alternative,save_out,
                               BF_out=TRUE,font_mltplyr=1, outlier_sd_thresh = NULL){
  library(ggplot2)
  library(tidyr)
  library(see)
  #remove outliers
  #browser()
  if (is.numeric(outlier_sd_thresh)) {
  main_data_old = as.data.frame(main_data)
  main_data = main_data_old[scale(main_data_old[[x]])<outlier_sd_thresh & scale(main_data_old[[y]])<outlier_sd_thresh, ]
  if (nrow(main_data) < nrow(main_data_old) ){
    print(paste0("removed", nrow(main_data) - nrow(main_data_old)," outlier(s)"))
  }
  }

  min_y_ax = min(na.omit(c(main_data[[x]],main_data[[y]])))
  max_y_ax = max(na.omit(c(main_data[[x]],main_data[[y]])))
  sd_y_ax = sd(na.omit(c(main_data[[x]],main_data[[y]])))
  data_long<-pivot_longer(main_data,
                          cols = c(x,y),
                          names_to = x_lab, values_to = y_lab)
  res<-t.test(main_data[[x]], main_data[[y]],  alternative = alternative, paired = TRUE)
  stats_text <- pub_ready_stats(res)
  if (BF_out==TRUE){
    library(BayesFactor)
    #annoyingly have to do listwise deletion by hand
    main_data <- main_data[complete.cases(main_data), ]
    BF = ttestBF(main_data[[x]], main_data[[y]], paired = TRUE)
    BF_val = round(as.vector(BF),digits = 2)
    stats_text = paste0(stats_text,', BF=',BF_val)
  }
  pretty_plot<-data_long %>%
    ggplot(aes_string(x = as.name(x_lab), y = as.name(y_lab))) +
    geom_violinhalf(flip = 1,width = 2,fill=alpha("black", 0.2)) +
    geom_point(alpha = .5)+
    geom_line(aes_string(group=sub_col_name),alpha = .3) +
    stat_summary(fun = "mean",geom = "crossbar", width = .15, colour = "red")+
    theme_classic()+
    theme(legend.position = "none",axis.title.x = element_blank()) +
    scale_x_discrete(expand = c(0, 1.1),
                     labels = cond_labs) +
    scale_y_continuous(limits = c(min_y_ax, max_y_ax+ (.5*sd_y_ax) )) +
           #+ .1*sd(c(na.omit(x),na.omit(y))))) +
    annotate(geom = "text",fontface = "bold",size = 2.7 * font_mltplyr, x =1.5, #1.5 is center
             y = Inf, # finding max of both
             label=stats_text, vjust = 1)+
    if (save_out == TRUE){
  ggsave(paste0(plot_title,".png"),width = 5.5, height = 5)
    print(getwd())
    }
  print(pretty_plot)
}

#######################################
#######################################
#######################################

mediate_bin_IV<-function(Y0, Y1, M0, M1, wide_df, sims){
  require(tidyr)
  colnames(wide_df)[colnames(wide_df)==Y0]<- "Y0"
  colnames(wide_df)[colnames(wide_df)==Y1]<- "Y1"
  colnames(wide_df)[colnames(wide_df)==M0]<- "M0"
  colnames(wide_df)[colnames(wide_df)==M1]<- "M1"
  long_df<-pivot_longer(wide_df,
                     cols = c(Y0,
                              Y1,
                              M0,
                              M1),
                     names_to = c(".value","X"),
                     names_pattern = "(.)(.)")
  a_to_b<-lm(M ~ X, data = long_df)
  a_to_c_cov_b<-lm(Y ~ X + M, data = long_df)
  res<-mediate(a_to_b,a_to_c_cov_b, treat = "X", mediator = "M", sims = sims)
  summary_res<-summary(res)
  return(summary_res)
}
