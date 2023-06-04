mcPeptides = function(text1, text2, repeats, potCsv, cEpitCsv, hEpitCsv, vsCsv, pref=""){
  Ssearch = read.csv(potCsv, na.strings = c("")) # read list of n-peptides that are absent in human peptidome
  cyto = read.csv(cEpitCsv,na.strings = c("")) # read list of Tc epitopes (Tarke)
  help = read.csv(hEpitCsv,na.strings = c("")) # read list of Th epitopes (Tarke)
  vsechny = read.csv(vsCsv, na.strings = c("")) # read list of all unique n-peptides present in the virus peptidome
  start = Sys.time()
  cytoFound = c()
  helpFound = c()
  notFound = c()
  matches = c()
  twoHitsFound = c()
  for(penta in Ssearch[!is.na(Ssearch$sequence),]$sequence) {
    c = grep(penta, cyto$sequence, value = T)
    h = grep(penta, help$sequence, value = T)
    if(length(c) != 0) {
      cytoFound[penta] = list(c)
    }
    if(length(h) != 0) {
      helpFound[penta] = list(h)
    }
    if(length(c) != 0 & length(h) != 0) {
      twoHitsFound[penta] = list(c)
    }
    if(length(c) != 0 || length(h) != 0) {
      matches = c(matches,c,h)
    }
    else {
      notFound = c(notFound, penta)
    }
  }
  ThE = length(helpFound)
  TcE = length(cytoFound)
  TTH = length(twoHitsFound)
  TE = ThE + TcE - TTH
  NE = length(notFound)
  end = Sys.time()
  cas1 = round(difftime(end, start, units="mins") * repeats * 2.16, digits = 2)
  print(paste("Odhad: ", cas1))
  varr = readline()
  for(penta in names(cytoFound)) {
    write(paste(penta,paste(cytoFound[[penta]],sep=",",collapse = ","),sep=","), "cytofound.csv", append = T)
  }
  for(penta in names(helpFound)) {
    write(paste(penta,paste(helpFound[[penta]],sep=",",collapse = ","),sep=","), "helpFound.csv", append = T)
  }
  
  for(penta in names(twoHitsFound)) {
    write(paste(penta,paste(twoHitsFound[[penta]],sep=",",collapse = ","),sep=","), "twoFound.csv", append = T)
  }
  all = union(cyto$sequence, help$sequence)
  unmatched = setdiff(all, matches)
  write(paste(unmatched,sep=",",collapse = ","), "notFound.csv")
  print(paste("ThE: ", ThE[1],", ", "TcE: ", TcE[1] ,", ",   "TE: ", TE, ",","Not found:", NE[1]))
  ThI = TcI = ThcI = TI = NI = c()
  for(i in  1:repeats){
    randPeptidy = vsechny$sequence[sample(1:length(vsechny$sequence), length(Ssearch$sequence[!is.na(Ssearch$sequence)]), replace = FALSE)]
    cat(" " ,i )
    cytoFound = c()
    helpFound = c()
    notFound = c()
    matches = c()
    twoHitsFound = c()
    for(penta in randPeptidy) {
      c = grep(penta, cyto$sequence, value = T)
      h = grep(penta, help$sequence, value = T)
      if(length(c) != 0) {
        cytoFound[penta] = list(c)
      }
      if(length(h) != 0) {
        helpFound[penta] = list(h)
      }
      if(length(c) != 0 & length(h) != 0) {
        twoHitsFound[penta] = list(c)
      }
      if(length(c) != 0 || length(h) != 0) {
        matches = c(matches,c,h)
      }
      else {
        notFound = c(notFound, penta)
      }
    }
    ThI = c(ThI, length(helpFound))
    TcI = c(TcI, length(cytoFound))
    TI = c(TI, length(helpFound) + length(cytoFound) - length(twoHitsFound))
    NI = c(NI,length(notFound))
  }
  ph = pc = phc = pt = pn = 0
  for(j in 1:repeats){
    ph = ifelse(ThI[j] > ThE[1], ph + 1, ph)
    pc = ifelse(TcI[j] > TcE[1], pc + 1,  pc)
    pt = ifelse(TI[j] > TE[1], pt + 1,  pt)
    pn = ifelse(NI[j] < NI[1], pn + 1, pn)
  }
  meanTh = mean(ThI)
  meanTc = mean(TcI)
  meanT = mean(TI)
  meanNI = mean(NI)
  sdTh = sd(ThI)
  sdTc = sd(TcI)
  sdT = sd(TI)
  sdNI = sd(NI)
  CohenTh =round((ThE[1])/sdTh,digits = 2)
  CohenT =round((TE[1])/sdT, digits = 2)
  CohenTc =round((TcE[1])/sdTc, digits = 2)
  CohenNI =round((NE[1])/sdNI,digits = 2)
  
  ph = ph/repeats
  pc = pc/repeats
  pt = pt/repeats
  pn = pn/repeats
  end = Sys.time()
  cas2 = round(difftime(end, start, units="mins"), digits = 2)
  sink(file="output3.csv", split=T)
  {
    print(paste("Text1: ", text1, ", ","Text2: ", text2, ", ", "Repeats: ", repeats,", ","Pool: ", vsCsv,", ","Potential targets: ", potCsv,
                ", ","Cytotoxic targets: ", cEpitCsv,", ", "Helper targets: ",hEpitCsv,", ", "Prefix: ",pref))
    print(paste("Th empiric: ", ThE[1],", ", "Tc empiric: ", TcE[1]      ,", ", "T empiric: ", TE, ",","Peptides without:", NE[1]))
    print(paste("Random Th", meanTh,", ","Random Tc",  meanTc, ", ", "Random T: ", meanT, ", ",  "random Pept. without: ", meanNI))
    print(paste("p Th: ", ph , ", ","p Tc: ", pc , ", ","p T: ", pt,  ",", "p without", pn))
    print(paste("sd d Th: ", sdTh , ", ","sd d Tc: ", sdTc , ", ","sd d T: ", sdT,  ",", "sd d without", sdNI))
    print(paste("Cohen d Th: ", CohenTh , ", ","Cohen d Tc: ", CohenTc , ", ","Cohen d T: ", CohenT,  ",", "Cohen d without", CohenNI, "Cas odhad: ", cas1, ",", "Cas: ", cas2))
  }
  unlink("output3.csv")
}
