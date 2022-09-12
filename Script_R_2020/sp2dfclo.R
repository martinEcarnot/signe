sp2dfclo <- function(sp,y1,y2,nam,lo)
{
  dat=data.frame(x=I(sp))
  class(dat$x)="matrix"
  if(!missing(nam)) rownames(dat$x)=nam
  if(!missing(lo)) colnames(dat)=lo
  dat=cbind(y1,y2,dat)
  colnames(dat)[1:2]=c(deparse(substitute(y1)),deparse(substitute(y2)))
  return(dat)
}
