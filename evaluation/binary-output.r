file <- file("../model/bin.dat", open="rb")

length <- readBin(con=file, what="integer", n=1)/4
readBin(con=file, what="integer", n=length)
close(file)

file <- file("../model/bin2.dat", open="rb") 
length <- readBin(con=file, what="integer", n=1)/8
readBin(con=file, what="numeric", n=length)
close(file)
