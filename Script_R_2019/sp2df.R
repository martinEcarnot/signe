sp2df <- function(sp,y,nam,lo)
{

  dat=data.frame(x=I(sp))
  class(dat$x)="matrix"
  if(!missing(nam)) rownames(dat$x)=nam
  if(!missing(lo)) colnames(dat)=lo
  dat=cbind(y,dat)
  return(dat)
}
