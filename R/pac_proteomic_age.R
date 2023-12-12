
# pac_proteomic_age: This function is used to calculate the PAC proteomic age 

pac_proteomic_age <- function(data){
  
  shape=0.148766187078631 #shape from gompertz model with age and 128 selected proteins   
  rate=0.000231560557596016 #rate from gompertz model with age and 128 selected proteins
  shape0=0.1200806 #shape from gompertz model with age only     
  rate0=0.000004885383 #rate from gompertz model with age only  
  beta_age=0.1087932 #coefficient associated with age in the gompertz model with age only                
  
  beta_age_protein=c(0.029353964,0.099491261,-0.128153586,0.097070584,
                     -0.005514185,0.191335791,-0.276924813,0.034399365,
                     -0.075457263,-0.005805245,-0.215725117,-0.064391821,
                     0.194400347,-0.014944780,0.099807081,0.050241198,
                     -0.029057558,-0.103814564,0.070480748,-0.034660872,
                     -0.155752380,-0.001621645,-0.015657010,0.119829451,
                     0.071034581,-0.050333113,-0.055156094,-0.078615313,
                     0.098018235,-0.006973108,-0.045531367,-0.227985467,
                     -0.137532908,0.090923928,-0.037373094,-0.108456609,
                     -0.062976739,-0.084424204,-0.048030320,0.167737797,
                     0.005878648,-0.007561524,-0.137882529,-0.016749283,
                     -0.018107792,-0.068519646,-0.027362969,-0.125079859,
                     0.026829529,0.138523985,-0.126913742,0.107668152,
                     0.209132471,0.096307256,-0.152706324,0.053938526,
                     0.110094224,-0.012193783,0.288576126,0.043023223,
                     -0.030544583,0.111604776,-0.001024753,0.028680301,
                     -0.013744572,-0.013033678,-0.139848216,0.006030185,
                     -0.159282896,-0.047343343,0.029428855,0.015070364,
                     0.068816751,-0.085180823,0.052069065,-0.208651159,
                     0.457163801,-0.136903943,-0.019177778,-0.083227652,
                     -0.080694831,-0.113440341,0.147377640,0.058012728,
                     -0.126449842,0.325071818,0.023575072,-0.132356723,
                     0.018121720,0.046941782,-0.048970791,-0.083726604,
                     -0.096547934,-0.071408718,0.111906716,-0.062176138,
                     0.027498932,0.116975605,0.123769062,0.068896784,
                     -0.137662829,-0.047219443,-0.056378508,0.265602762,
                     -0.115249658,-0.080153962,0.014768758,-0.045588135,
                     0.069474034,0.037816669,-0.055014204,0.067015239,
                     -0.104542859,-0.123590138,0.164663205,0.130347786,
                     0.123681499,0.118673810,0.031169845,0.063336979,
                     0.034805482,-0.098905845,-0.144461235,-0.127344390,
                     -0.020352303,0.079187131,-0.305091894,0.118800770,-0.089338440
  ) #coefficients associated with age and proteins in the gompertz model with age and 128 selected proteins
  
  names_beta=c("age_recruitment","ada2", "adamts13", "adamts16", "adgrg2",
               "adm", "ager", "agr2", "apoe", "areg", "art3", "bag3",
               "bcam", "bcan", "brk1", "c7", "ca14", "ca4", "calca", 
               "cblif", "cd248", "cd34", "cdcp1", "ceacam5", "ceacam6", 
               "cela2a", "chad", "cndp1", "col24a1", "cpa2", "crh", "crisp3", 
               "crtac1", "cthrc1", "ctrc", "dcxr", "ddc", "dkk4", "dsg4",
               "eda2r", "egfr", "epha4", "fap", "faslg", "fcrl1", "fgfbp1",
               "fgfbp2", "flt3", "fut3_fut5", "galnt5", "galnt7", "gdf15",
               "gfap", "gp2", "gzma", "havcr1", "hgf", "hpgds", "icam3",
               "icam5", "igfbp3", "igfbp7", "il1rl1", "il22", "il6", "itga2",
               "itgam", "itgav", "izumo1", "kit", "klk3", "klk4", "krt19", 
               "lcp1", "lrrn1", "lrtm2", "ltbp2", "lto1", "mep1b", "mepe", 
               "met", "mmp10", "mmp12", "mmp9", "ncan", "nefl", "nppb", "ntf4",
               "ntprobnp", "ocln", "odam", "optc", "pepd", "prg3", "prss8", 
               "ptgr1", "ptn", "rbfox3", "reg3g", "ren", "scg2", "scg3", 
               "scgb1a1", "scgb3a1", "sdc4", "selenop", "sell", "sez6l", 
               "sfrp1", "sftpd", "skap1", "slitrk1", "slitrk2", "sost", "spock1",
               "spon1", "spp1", "tff3", "tnc", "tnfrsf10b", "tnfrsf6b", 
               "tnn", "tnr", "tpk1", "ttr", "txndc15", "vgf", "wfdc2", "xg") # variable names corresponding to the coefficients in "beta_age_protein"
  
  # match the variable names in the input data name and "names_beta"
  beta_age_protein=beta_age_protein[match(rownames(data), names_beta)]
  
  # b(x)=b*exp(x*beta)
  b_x<- t(as.matrix(data))%*%as.matrix(beta_age_protein)
  rate_new <- rate*exp(b_x)
  
  # 10-year mortality risk 
  cdf_10_year=1-exp(-(rate_new/shape)*(exp(shape*10)-1))
  
  # PAC proteomic age
  proteomic_age=(1/beta_age)*log(shape0*log(1-cdf_10_year)/(rate0*(1-exp(10*shape0))))
  colnames(proteomic_age)="proteomic age"
  return(proteomic_age)
}




