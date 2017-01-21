rm(list=ls())

#par.vals - what you want to change it to. A matrix (use expand.grid)
#oat.nams - changed from. A matrix, same size
code.generate = function(par.vals, pat.nams, path="", path.opt=F,
                            templatename, out_name){ 
  count = 0
  for(i in 1:nrow(par.vals)){
    SourceCode=readLines(templatename) #read desired code

    for(j in 1:ncol(par.vals)){
      #replace with desired value
      SourceCode=gsub(pattern=pat.nams[i,j],replacement=par.vals[i,j],x=SourceCode) 
    }
    count = count + 1
    #output changed code
    path=ifelse(path==""|path.opt==F,"",paste(path,"/",sep=""))
    cat(SourceCode,file=paste(path, out_name,
                              paste(count,collapse=""),
                              '.R',sep=""),sep='\n')  
  }
}  

######################################################################

numvals = expand.grid(list(c(0.001,0.01,0.1), c(0.001,0.01,0.1)))
nams = cbind(rep("AA",times=nrow(numvals)),
    rep("BB",times=nrow(numvals)))

code.generate(numvals,nams,templatename='admm2.r', out_name='admm2_')


