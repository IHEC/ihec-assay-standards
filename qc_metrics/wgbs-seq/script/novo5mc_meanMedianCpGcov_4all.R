args=commandArgs(TRUE)
print(args)
for(j in 1:length(args)){
input<-args[j]
t<-read.table(input)
av<-mean(t$V2)
med<-median(t$V2)
print(input)
print(av)
print(med)
}
