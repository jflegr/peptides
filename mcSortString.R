mcSortStrings = function(text1, text2, repeats, s1_long, s2_short, s3_pool, nameFound = "stringsWithSubstrings.csv", nameNoFound = "stringsWithoutSubstrings.csv"){
  start = Sys.time()
  countUnfound = function(shorter, longer) {
    result=list()
    foundTotal = c()
    for(needle in shorter) {
      found = grep(needle, longer, value=T)
      foundTotal = union(found, foundTotal)
    }
    result$unfoundCount = length(longer) - length(foundTotal)
    result$found = foundTotal
    result$unfound = setdiff(longer, foundTotal)
    result$foundCount = length(foundTotal)
    return(result)
  }
  help = read.csv(s1_long,na.strings = c("")) # read list of strings, e.g. Th epitopes (Tarke)
  Ssearch = read.csv(s2_short, na.strings = c("")) # read list substring to be search for, e.g. n-peptides that are absent in human peptidome
  vsechny = read.csv(s3_pool, na.strings = c("")) # read list of all substrings for random selection, e.g. all unique n-peptides present in the virus peptidome
  result = countUnfound(Ssearch$sequence, help$sequence)
  nEmp = result$unfoundCount[1]
  fEmp = result$foundCount[1]
  write(paste(result$found,sep=",",collapse = ","), nameFound)
  write(paste(result$unfound,sep=",",collapse = ","), nameNoFound)
  ThI = c()
  NI = c()
  FI = c()
  for(i in  1:repeats){
    randPeptidy = vsechny$sequence[sample(1:length(vsechny$sequence), length(help$sequence), replace = FALSE)]
    result = countUnfound(Ssearch$sequence, randPeptidy)
    NI = c(NI,result$unfoundCount)
    FI = c(FI, result$foundCount)
  }
  ph  = pn = pf = 0
  for(j in 1:repeats){
    pn = ifelse(nEmp >= NI[j] , pn + 1, pn)
    pf = ifelse(fEmp <= FI[j], pf + 1, pf)
  }
  meanNI = mean(NI)
  sdNI = sd(NI)
  meanFI = mean(FI)
  sdFI = sd(FI)
  cohenDNI = (meanNI - nEmp)/sdNI
  cohenDFI = (meanFI - fEmp)/sdFI
  pn = pn/repeats
  pf = pf/repeats
  end = Sys.time()
  cas = round(difftime(end, start, units="mins"), digits = 2)
  {
    print(paste("What: ", text1, ", ","Length: ", text2, ", ", "repeats: ", repeats, ", ","Longer: ", s1_long,", ","Shoter: ", s2_short, ", ","Pool: ", s3_pool,   "Output with: ",  nameFound, ", ", "Output without: ", nameNoFound,", ", "Time: ", cas))
    print(paste("Total: ", length(help$sequence),", ","Strings without:", nEmp, ", ","Random without: ", meanNI, ", ","pn: ", pn,", ", "SD random not found: ", sdNI,", ", "Cohen d NF: ", cohenDNI))
    print(paste("Total: ", length(help$sequence),", ","Strings with:", fEmp, ", ","Random with: ", meanFI, ", ","pf: ", pf,", ", "SD random found: ", sdFI,", ", "Cohen d found: ", cohenDFI))
  }
  sink(file="outputMC.csv", split=T)
  {
    print(paste("What: ", text1, ", ","Length: ", text2, ", ", "Repeats: ", repeats, ", ","Longer: ", s1_long,", ","Shoter: ", s2_short, ", ","Pool: ", s3_pool,  "Output with: ",  nameFound, ", ", "Output without: ", nameNoFound,", ", "Time: ", cas))
    print(paste("Total: ", length(help$sequence),", ","Strings without:", nEmp, ", ","Random without: ", meanNI, ", ","pn: ", pn,", ", "SD random not found: ", sdNI,", ", "Cohen d NF: ", cohenDNI))
    print(paste("Total: ", length(help$sequence),", ","Strings with:", fEmp, ", ","Random with: ", meanFI, ", ","pf: ", pf,", ", "SD random found: ", sdFI,", ", "Cohen d found: ", cohenDFI))
  }
  unlink("outputMC.csv")
}
