#!/usr/bin/env Rscript
source("PCAplot.R")


library("pls")

#############
#  quality  #
#############


perftable = function (list_predict, list_real){
  nb_value = length (list_real)
  i = 1
  tp = 0
  fp = 0
  tn = 0
  fn = 0
  while(i <= nb_value ){
    if (as.character(list_predict[i])=="d"){
      if (list_predict[i] == list_real[i]){
        tp = tp + 1
      }else {
        fp = fp + 1
      }
    }else{
      if (list_predict[i] == list_real[i]){
        tn = tn + 1
      }else {
        fn = fn + 1
      }
    }
    i = i + 1
  }
  print (paste ("TP : ", tp, sep = ""))
  print (paste ("TN : ", tn, sep = ""))
  print (paste ("FP : ", fp, sep = ""))
  print (paste ("FN : ", fn, sep = ""))
  
  tableval = c(tp,tn,fp,fn)
  return (tableval)
}


#####################
# Perf for classes  #
#####################

accuracy = function (tp, tn, fp, fn){
  return ((tp + tn)/(tp + fp + tn +fn))
}

precision = function (tp, fp){
  return (tp/(tp + fp))
}

recall = function (tp, fn){
  return (tp/(tp + fn))
}

specificity = function (tn, fp){
  return (tn/(tn + fp))
}

sensibility = function (tp, fn){
  return (tp/(tp + fn))
}


BCR = function (tp, tn, fp, fn){
  return (0.5*(tp/(tp+fn) + tn/(tn+fp)))
  
}

MCC = function (tp, tn, fp, fn){
  numerator = tp*tn-fp*fn
  denumerator = (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
  return (numerator / sqrt(denumerator))
}


qualityPredict = function (predict, Y2){
  print (as.vector(predict)[[1]])
  print (as.vector(Y2)[[1]])
  v_predict = calculTaux (as.vector(predict)[[1]], as.vector(Y2)[[1]])
  print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
  #print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  return (v_predict)
}




qualityPredictList = function (test_vector, real_vector){
  v_predict = calculTaux (test_vector, real_vector)
  print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
  #print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  return (v_predict)
}


qualityShowModelSelection = function (list_des_model, coef, v_real_train, v_predict_train, v_real_test, v_predict_test, v_real_loo, v_predict_loo, l_out_CV){
  
  # loo
  criteria_loo = computedCriteria(v_real_loo, v_predict_loo)
  # CV
  criteria_CV = computedCriteriaCV(l_out_CV)
  # train
  criteria_train = computedCriteria(v_real_train, v_predict_train)
  # test
  criteria_test = computedCriteria(v_real_test, v_predict_test)
  
  # show
  print ("descriptor")
  print (list_des_model)
  print (as.vector(abs(coef[list_des_model])))
  print ("Acc_loo --- Acc_train --- Acc_test --- Acc_CV_train --- SD_CV_train --- Acc_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[1]], criteria_train[[1]], criteria_test[[1]], criteria_CV["acc_train"], criteria_CV["acc_train_SD"], criteria_CV["acc_test"], criteria_CV["acc_test_SD"], sep = "---"))
  print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test --- Se_CV_train --- SD_CV_train --- Se_CV_test --- SD_CV_test --- Sp_CV_train --- SD_CV_train --- Sp_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[2]], criteria_loo[[3]], criteria_train[[2]], criteria_train[[3]], criteria_test[[2]], criteria_test[[3]], criteria_CV["se_train"], criteria_CV["se_train_SD"], criteria_CV["se_test"],criteria_CV["se_test_SD"], criteria_CV["sp_train"], criteria_CV["sp_train_SD"], criteria_CV["sp_test"],criteria_CV["sp_test_SD"], sep = "---"))
  print ("MCC_loo --- MCC_train --- MCC_test --- MCC_CV_train --- SD_CV_train --- MCC_CV_test --- SD_CV_test")
  print (paste (criteria_loo[[4]], criteria_train[[4]], criteria_test[[4]], criteria_CV["mcc_train"], criteria_CV["mcc_train_SD"], criteria_CV["mcc_test"], criteria_CV["mcc_test_SD"], sep = "---"))
  print ("**********************************************************************")
}



computedCriteria = function (v_real, v_predict){
  
  rate = calculTaux2 (v_predict, v_real)
  acc = accuracy(rate[1], rate[2], rate[3], rate[4])
  se = sensibility(rate[1], rate[4])
  sp = sensibility(rate[2], rate[3])
  mcc = MCC(rate[1], rate[2], rate[3], rate[4])
  
  return (list (acc, se, sp, mcc))
  
}


computedCriteriaCV = function (l_out_CV){
  
  CV_train = l_out_CV[[1]]
  CV_test = l_out_CV[[2]]
  
  v_acc_train = NULL
  v_acc_test = NULL
  v_se_train = NULL
  v_se_test = NULL
  v_sp_train = NULL
  v_sp_test = NULL
  v_mcc_train = NULL
  v_mcc_test = NULL
  
  
  for (i in seq (1,length (CV_train))){
    v_acc_train = append (v_acc_train, CV_train[[i]][[1]])
    v_acc_test = append (v_acc_test, CV_test[[i]][[1]])
    v_se_train = append (v_se_train, CV_train[[i]][[2]])
    v_se_test = append (v_se_test, CV_test[[i]][[2]])
    v_sp_train = append (v_sp_train, CV_train[[i]][[3]])
    v_sp_test = append (v_sp_test, CV_test[[i]][[3]])
    v_mcc_train = append (v_mcc_train, CV_train[[i]][[4]])
    v_mcc_test = append (v_mcc_test, CV_test[[i]][[4]])
  }
  
  v_out = c(mean (v_acc_train), sd (v_acc_train), mean (v_acc_test), sd (v_acc_test), mean (v_se_train), sd (v_se_train), mean (v_se_test), sd (v_se_test), mean (v_sp_train), sd (v_sp_train), mean (v_sp_test), sd (v_sp_test),mean (v_mcc_train), sd (v_mcc_train), mean (v_mcc_test), sd (v_mcc_test) )
  names (v_out) = c("acc_train", "acc_train_SD","acc_test", "acc_test_SD","se_train", "se_train_SD", "se_test", "se_test_SD", "sp_train", "sp_train_SD", "sp_test", "sp_test_SD", "mcc_train", "mcc_train_SD", "mcc_test", "mcc_test_SD")
  return (v_out)
}



cumulTaux = function (taux1, taux2){
  
  tp = taux1[1] + taux2[1]
  tn = taux1[2] + taux2[2]
  fp = taux2[3]
  fn = taux2[4]
  
  print (paste ("accuracy : ", accuracy(tp, tn, fp, fn), sep = ""))
  print (paste ("precision : ",precision(tp, fp), sep = ""))
  #print (paste ("recall : ", recall(tp, fn), sep = ""))
  print (paste ("sensibility : ", sensibility(tp, fn), sep = ""))
  print (paste ("specificity : ", sensibility(tn, fp), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (tp, tn, fp, fn), sep = ""))	
  print (paste ("BER (balanced error rate) : ", 1 - BCR (tp, tn, fp, fn), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (tp, tn, fp, fn), sep = ""))	
  
}


# for ROC curve -> calcul vecteur prediction with probability (just for druggability)
generateVect = function(proba_out_predict, threshold){
  
  proba_class1 = proba_out_predict[,1]
  
  vect_out = NULL
  
  for (proba in proba_class1){
    if (proba > threshold){
      vect_out = c(vect_out, "d")
    }else{
      vect_out = c(vect_out, "nd")
    }
  }
  return (vect_out)
}


#########################
#    PERF regression    #
#########################


vrmsep = function(dreal, dpredict){
  
  #dpredict = dpredict[rownames(dreal),]
  
  i = 1
  imax = length(dreal)
  
  valout = 0
  while(i <= imax){
    valout = valout + ((dreal[i] - dpredict[i])^2)
    i = i + 1
  }
  return(sqrt(valout))
  
}

###############################
# divise the dataset in folds #
###############################

samplingDataNgroupClass = function (t_din, i_nb_group, s_nameclass){
  
  # divise two classes
  v_class = as.factor(t_din[,s_nameclass])
  t_dc0 = t_din [which(v_class == 0),]
  t_dc1 = t_din [which(v_class == 1),]
  
  
  # sample data
  v_sampledc0 = sample (dim (t_dc0)[1])
  v_sampledc1 = sample (dim (t_dc1)[1])
  
  # ind limit
  i_limitc0 = as.integer (dim(t_dc0)[1] / i_nb_group)
  i_limitc1 = as.integer (dim(t_dc1)[1] / i_nb_group)
  
  #print (i_limitc0)
  #print (i_limitc1)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = rbind (t_dc0[v_sampledc0[1:i_limitc0],], t_dc1[v_sampledc1[1:i_limitc1],])
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = rbind (t_dc0[v_sampledc0[((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))],], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(length (v_sampledc1))]),])
    }
    else{
      t_group = rbind (t_dc0[(v_sampledc0[(i_limitc0 * (g-1) + 1):(i_limitc0 * g)]),], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(i_limitc1 * g)]),])
    }
    # append list
    output[[g]] = t_group
  }
  
  return (output)
}




samplingDataFraction = function (t_din, fract){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = round((dim (t_din)[1]) * fract)
  
  dtrain = t_din[v_sample[(i_limitc + 1):(length (v_sample))],]
  dtest = t_din[v_sample[1:i_limitc],]
  
  return (list(dtrain, dtest))
  
}




samplingDataNgroup = function (t_din, i_nb_group){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = as.integer (dim(t_din)[1] / i_nb_group)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = t_din[v_sample[1:i_limitc],]
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = t_din[v_sample[((i_limitc * (g-1) + 1):(length (v_sample)))],]
    }
    else{
      t_group = t_din[(v_sample[(i_limitc * (g-1) + 1):(i_limitc * g)]),]
    }
    # append list
    output[[g]] = t_group
  }
  
  return (output)
}




#################################
#   Control dataset integrity   #
#################################

controlDatasets = function(ldataset, prin){
  
  nbsplit = length(ldataset)
  
  print (nbsplit)
  
  pdf(paste(prin,".pdf", sep = ""), width = 10, height = 10)
  
  colorrainbow = rainbow(nbsplit)
  i = 1
  colorpoint = NULL
  dPCA = NULL
  for (d in ldataset){
    hist(d[,(dim(d)[2])], col = "grey", main = paste("Fold ", i, "-Dim=", dim(d)[1], sep = ""))

    #points for PCA
    colorpoint = append(colorpoint, rep(colorrainbow[i], (dim(d)[1])))
    dPCA = rbind(dPCA, d[,-dim(d)[2]]) # remove affinity coulum
    i = i + 1
  }
  
  # plot PCA
  dplot = generatePCAcoords(dPCA)
  var_cap = dplot[[2]]
  data_plot = dplot[[1]]
  #par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = colorpoint, xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 2, cex.main = 2, cex.axis = 1.75, cex = 2)
  abline(h=0,v=0)
  dev.off()
  
}


########################
#    MODEL regression  #
########################


# PCR - cross validation
PCRCV = function(lfolds, prout){
  
  # performance in CV and number of cmp
  # return best number of cmp
  
  maxCp = dim(lfolds[[1]])[2]
  
  vcpRMSEP = NULL
  vcpR2 = NULL
  for (cp in seq(1,maxCp)){
    i = 1
    imax = length(lfolds)
    vpred = NULL
    vreal = NULL
    while(i <= imax){
      dtrain = NULL
      dtest = NULL
      for (j in seq(1:imax)){
        if (j == i){
          dtest = lfolds[[j]]
        }else{
          dtrain = rbind(dtrain, lfolds[[j]])
        }
      }
      modelpcr = pcr(Aff~., data=dtrain, ncomp = cp)
      predpcr = predict(modelpcr, ncomp = cp, newdata = dtest)
      
      vpred = append(vpred, predpcr)
      vreal = append(vreal, dtest[,"Aff"])
      i = i + 1
    }
    
    corpred = cor(vreal, vpred)
    rmsepcp = vrmsep(vreal, vpred)
    
    vcpR2 = append(vcpR2, corpred)
    vcpRMSEP = append(vcpRMSEP,rmsepcp)
  }
  
  #print(vcpR2)
  #print (vcpRMSEP)
  
  png(width = 960, height = 480, filename = paste(prout ,"CVpcr.png", sep = ""))
  par(mfrow = c(1,2))
  plot(seq(1,length(vcpR2)), vcpR2, type = "l", cex = 3, main = paste("R2 by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "R2")
  plot(seq(1,length(vcpRMSEP)), vcpRMSEP, type = "l", cex = 3, main = paste("RMSEP by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "RMSEP")
  dev.off()
  
  # optimal number of component
  outcp = cbind(seq(1, length(vcpRMSEP)), vcpRMSEP)
  dout = NULL
  for (i in seq(1,length(vcpRMSEP))){
    dout = append(dout, eucdist(0, outcp[i,1], 0, outcp[i,2]))
  }
  outcp = cbind(outcp, dout)
  nbCPoptimun = outcp[which(outcp[,3] == min(outcp[,3])),1]
  
  print("****PCR in CV****")
  print(paste("Optimal component: ", nbCPoptimun, sep = ""))
  print(paste("Perf optimal in CV: RMSEP=", vcpRMSEP[nbCPoptimun], " R2=", vcpR2[nbCPoptimun], sep = ""))
  
  return (nbCPoptimun)
  
}


# PCR - real
PCRTrainTest = function(dtrain, dtest, nbcp){
  
  modelpcr = pcr(Aff~., data=dtrain, ncomp = nbcp)
  predpcrtest = predict(modelpcr, ncomp = nbcp, newdata = dtest)
  predpcrtrain = predict(modelpcr, ncomp = nbcp, newdata = dtrain)
  
  r2train = cor(dtrain[,"Aff"], predpcrtrain)
  r2test = cor(dtest[,"Aff"], predpcrtest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predpcrtrain)
  rmseptest = vrmsep(dtest[,"Aff"], predpcrtest)
  
  
  print("****PCR model****")
  print(paste("NB components = ", nbcp, sep = ""))
  print(paste("Perf training (dim = ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("- R2=", r2train))
  print(paste("- RMSEP=", rmseptrain))
  
  print(paste("Perf test (dim = ", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("- R2=", r2test))
  print(paste("- RMSEP=", rmseptest))
  
  return(modelpcr)
  
}



# PLS in cross validation
PLSCV = function(lfolds, prout){
  
  # performance in CV and number of cmp
  # return best number of cmp
  
  maxCp = dim(lfolds[[1]])[2] -1 # remove col with affinity
  
  vcpRMSEP = NULL
  vcpR2 = NULL
  for (cp in seq(1,maxCp)){
    i = 1
    imax = length(lfolds)
    vpred = NULL
    vreal = NULL
    while(i <= imax){
      dtrain = NULL
      dtest = NULL
      for (j in seq(1:imax)){
        if (j == i){
          dtest = lfolds[[j]]
        }else{
          dtrain = rbind(dtrain, lfolds[[j]])
        }
      }
      modelpls = plsr(Aff~., data=dtrain, ncomp = cp)
      predpls = predict(modelpls, ncomp = cp, newdata = dtest)
      
      vpred = append(vpred, predpls)
      vreal = append(vreal, dtest[,"Aff"])
      i = i + 1
    }
    
    corpred = cor(vreal, vpred)
    rmsepcp = vrmsep(vreal, vpred)
    
    vcpR2 = append(vcpR2, corpred)
    vcpRMSEP = append(vcpRMSEP,rmsepcp)
  }
  
  #print(vcpR2)
  #print (vcpRMSEP)
  
  png(width = 960, height = 480, filename = paste(prout ,"CVpls.png", sep = ""))
  par(mfrow = c(1,2))
  plot(seq(1,length(vcpR2)), vcpR2, type = "l", cex = 3, main = paste("R2 by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "R2")
  plot(seq(1,length(vcpRMSEP)), vcpRMSEP, type = "l", cex = 3, main = paste("RMSEP by components in CV", length(lfolds), sep = ""), xlab = "Components", ylab = "RMSEP")
  dev.off()
  
  # optimal number of component
  outcp = cbind(seq(1, length(vcpRMSEP)), vcpRMSEP)
  dout = NULL
  for (i in seq(1,length(vcpRMSEP))){
    dout = append(dout, eucdist(0, outcp[i,1], 0, outcp[i,2]))
  }
  outcp = cbind(outcp, dout)
  nbCPoptimun = outcp[which(outcp[,3] == min(outcp[,3])),1]
  
  print("****PLS in CV****")
  print(paste("Optimal component: ", nbCPoptimun, sep = ""))
  print(paste("Perf optimal in CV: RMSEP=", vcpRMSEP[nbCPoptimun], " R2=", vcpR2[nbCPoptimun], sep = ""))
  
  return (nbCPoptimun)
  
}


# PLS - real
PLSTrainTest = function(dtrain, dtest, nbcp){
  
  modelpls = plsr(Aff~., data=dtrain, ncomp = nbcp)
  predplstest = predict(modelpls, ncomp = nbcp, newdata = dtest)
  predplstrain = predict(modelpls, ncomp = nbcp, newdata = dtrain)
  
  r2train = cor(dtrain[,"Aff"], predplstrain)
  r2test = cor(dtest[,"Aff"], predplstest)
  
  rmseptrain = vrmsep(dtrain[,"Aff"], predplstrain)
  rmseptest = vrmsep(dtest[,"Aff"], predplstest)
  
  
  print("****PLS model****")
  print(paste("NB components = ", nbcp, sep = ""))
  print(paste("Perf training (dim= ", dim(dtrain)[1], "*", dim(dtrain)[2], "):", sep = ""))
  print(paste("- R2=", r2train))
  print(paste("- RMSEP=", rmseptrain))
  
  print(paste("Perf test (dim = ", dim(dtest)[1], "*", dim(dtest)[2], "):", sep = ""))
  print(paste("- R2=", r2test))
  print(paste("- RMSEP=", rmseptest))
  
  return(modelpls)
  
}


