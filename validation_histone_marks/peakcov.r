options(scipen = 10000)

# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Set arguments to default values
bedfile = NA # input motif
out = NA #output file
distance = 4000




#Parsing arguments and storing values
for (each.arg in args) {
	#input bed file
	if (grepl('^-intersect=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				bedfile <- arg.split[2]
		} else {
			stop('No input file')
		} 
	}

	if (grepl('^-out=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				out <- arg.split[2]
		} else {
			stop('no out')
		} 
	}

}

#=======================> DONE! 


inter = read.table(bedfile)
result = vector(mode = "numeric", length=distance)

genesuniq = length(unique(inter$V14))
peaksuniq = length(unique(inter$V4))

for (i in 1:length(inter$V1)) {
    
    
    
    if (inter$V16[i] == "+") {
		starting = inter$V2[i] - inter$V12[i]
		ending = inter$V3[i] - inter$V12[i]
    }
    
    if (inter$V16[i] == "-") {
		ending = distance - (inter$V2[i] - inter$V12[i])
		starting = distance - (inter$V3[i] - inter$V12[i])
    }
    
    starting[starting < 1] = 1
    ending[ending > distance] = distance
    
 
	pos1 = starting
	pos2 = ending
	result[pos1:pos2] = result[pos1:pos2] + 1
}

write((result / mean(result)), file = paste0(out, "/result.peakcov"), ncolumns = 1)
write(paste(genesuniq, peaksuniq, peaksuniq/genesuniq, sep = "	"), file = paste0(out, "/result.sens"), ncolumns = 3)
#write(paste(distan, sep = "	"), file = paste0(out, "/result.distance"), ncolumns = 1)
