#Script en R para limpipar datos:
inicio_scrp <- as.numeric(Sys.time())
cnames = c("CHR","SNP","BP","A1","A2","FREQ","BETA","SE","P")
#args <- commandArgs(trailingOnly=TRUE)
setwd('C:/Users/luis_/Dropbox (Personal)/Documentos/Labortorio/Gabriel_Cuellar/ManhattanPlot/files')
gwas <- "../files/DataReal_filtrado.txt"#args[1]
out <- "../files/DataReal_filtrado.txt.fltdR"#args[2]
factor_division <- 1
xPx <-1148
yPx <- 649


library(data.table)

data <- fread(gwas, data.table=FALSE);

checkColumns <- all(cnames %in% colnames(data));
if(!checkColumns) {
  missingCols <- which(!(cnames %in% colnames(data)))
  message <- paste(cnames[missingCols], collapse=" and ")
  stop(paste("Error: Columns", message, "are missing from your data." ))
}

data <- data[complete.cases(data),]
if(nrow(data)<10) {
  stop("Many missing data or very few SNPs in data file")
}

checkFreq <- all(is.numeric(data$FREQ))
if(!checkFreq) {
  stop("FREQ column should contain only numbers")
}

checkFreqBounds <- all(data$FREQ >0 & data$FREQ <1)
if(!checkFreqBounds) {
  stop("FREQ column should contain only numbers between 0 and 1 exclusive")
}

checkBeta <- all(is.numeric(data$BETA))
if(!checkBeta) {
  stop("BETA column should contain only numbers")
}

checkSE <- all(is.numeric(data$SE))
if(!checkSE) {
  stop("SE column should contain only numbers")
}

checkP <- all(is.numeric(data$P))
if(!checkP) {
  data$P <- as.numeric(data$P)
  # stop("P-value column should contain only numbers")
}

checkPBounds <- all(data$P <1 & data$P >0)
if(!checkPBounds) {
  stop("P-value column should contain only numbers between 0 and 1 inclusive")
}
print(head(data))




inicio_algrt <- as.numeric(Sys.time())

chrm <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',
          '16','17','18','19','20','21','22','X','Y')
posc <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
  145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,
  90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
dataF_chrmPosic <- data.frame(chrm,posc,stringsAsFactors=FALSE)

chrm <- unique(data$CHR)
order_chrm <- order(as.numeric(gsub("[a-zA-Z]*(\\d+)", "\\1", chrm)))
posc<-as.vector(sapply(chrm,function(x)
{dataF_chrmPosic$posc[which(dataF_chrmPosic$chrm==gsub('[a-zA-Z]*(\\d+)','\\1',x))]}))
chrm <- chrm[order_chrm]
posc <- posc[order_chrm]
dataF_chrmPosic <- data.frame(chrm,posc,stringsAsFactors=FALSE)


xLen <- as.integer(xPx/factor_division)
yLen <- as.integer(yPx/factor_division)

vMax <- -log10(min(data$P))
vMin <- -log10(max(data$P))
pbPerPixel_X = sum(dataF_chrmPosic$posc)/xLen
valoresPerPixel_y  = (vMax-vMin)/yLen

posAnterior <- 0 
for (i in 1:nrow(dataF_chrmPosic)){
  tamChrm <- dataF_chrmPosic$posc[i]
  dataF_chrmPosic$posc[i] <- posAnterior
  posAnterior <- posAnterior + tamChrm
}


  
#vect_poscChrm<-as.vector(sapply(data$CHR,function(x)
 # {dataF_chrmPosic$posc[which(dataF_chrmPosic$chrm==x)]}))


matrix_indexIn <- matrix(nrow = yLen+1,ncol = xLen+1,byrow = T)
matrix_values <- matrix(data = 0,nrow = yLen+1,ncol = xLen+1,byrow = T)


inicio <- as.numeric(Sys.time())
vec_PosChr <- sapply(data$CHR,function(row)dataF_chrmPosic$posc[which(dataF_chrmPosic$chrm==row)])
print(as.numeric(Sys.time())-inicio)
#43.66129



vec_PosAbs <- data$BP + vec_PosChr
data[['AbslPost']] <- vec_PosAbs
vec_log <- -log10(data$P)
#data[['Log']] <- vec_log



inicio <- as.numeric(Sys.time())
for (i in 1:nrow(data)){
  row = data[i,]
  poscAbs <- row$AbslPost
  l <- -log10(row$P)
  #l <- row$Log
  x <-as.integer(poscAbs/pbPerPixel_X)+1
  y <- as.integer(l/valoresPerPixel_y) +1
  #print(c(y,x))
  if (matrix_values[y,x] < l){ 
    matrix_values[y,x] <- l
    matrix_indexIn[y,x] <- i
  }
}
print(as.numeric(Sys.time())-inicio)
#145.6163,125.4744


vec_indexSNPs <- as.vector(matrix_indexIn)[!is.na(as.vector(matrix_indexIn))]

cnames <- c(cnames,'AbslPost')
print(c(as.numeric(Sys.time())-inicio_algrt,as.numeric(Sys.time())-inicio_scrp,length(vec_indexSNPs)))
write.table(data[vec_indexSNPs,cnames], file=out, quote=F, row.names=F, col.names=T, sep="\t")



#################################################################################################









vec_xs <- as.integer(vec_PosAbs/pbPerPixel_X)+1
vec_ys <- as.integer(vec_log/valoresPerPixel_y) +1
mtx_datos <- matrix(c(1:nrow(data),vec_xs,vec_ys,vec_log),nrow = nrow(data),ncol = 4)

inicio <- as.numeric(Sys.time())
r <- apply(mtx_datos,1,function(x){
  if (matrix_values[x[3],x[2]] < x[4]){ 
    matrix_values[x[3],x[2]] <- x[4]
    matrix_indexIn[x[3],x[2]] <- x[1]
  }
})
as.numeric(Sys.time())-inicio
#334.4861


a <- by(data[1:10,], 1:10, function(row){
  lo <- -log10(row$P)
  c(as.integer(dataF_chrmPosic$posc[which(dataF_chrmPosic$chrm==row$CHR)]/pbPerPixel_X),
    as.integer(lo/valoresPerPixel_y),lo)
});a

lapply(a,FUN = function(x){
  xC <- as.integer(x[1])
  yC <- as.integer(x[2])
  v <- x[3]
  print(c(xC,yC,v,x))
  if(matrix_values[xC,yC] < v){
    matrix_values[xC,yC] <- v
    matrix_indexIn[xC,yC] <- as.integer(names(x))
  }
})


inicio <- as.numeric(Sys.time())
a <- by(data, 1:nrow(data), function(row){
  poscAbs <- row$BP+dataF_chrmPosic$posc[which(dataF_chrmPosic$chrm==row$CHR)]
  l <- -log10(row$P)
  x <-as.integer(poscAbs/pbPerPixel_X)+1
  y <- as.integer(l/valoresPerPixel_y) +1
  if (matrix_values[y,x] < l){ 
    matrix_values[y,x] <- l
    #matrix_indexIn[y,x] <- i
  }
})
as.numeric(Sys.time())-inicio
#149.606
