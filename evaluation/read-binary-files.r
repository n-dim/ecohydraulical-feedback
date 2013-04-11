nSteps <- 100
filepath <- file.choose()
bytes <- 4

readBinary <- function(filepath, nSteps, bytes){
  file <- file(filepath, open="rb")
  for(i in 1:nSteps){
    print(length <- readBin(con=file, what="integer", n=1)/bytes )#first 4 bytes indicate how many binarys follow
    print(readBin(con=file, what="integer", n=length))
    length <- readBin(con=file, what="integer", n=1)/bytes #there is another index number after each block (junk for us)
  }
  
  close(file)
}
readBinary(filepath, nSteps, bytes)