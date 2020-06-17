### Output Latex

meaStdTabLat <- function(mea, std, rowNam, s){
  numCol <- ncol(mea)
  numRow <- nrow(mea)
  tab <- character(length = numRow)
  for(i in 1:numRow){
    texLin <- rowNam[i]
    for(j in 1:numCol){
      texLin <- paste0(texLin, " & $", signif(mea[i,j], s), "_{", signif(std[i,j], s), "}$")
    }
    texLin <- paste0(texLin, " \\\\")
    cat(texLin, "\n")
  }
}
