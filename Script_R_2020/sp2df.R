sp2df <- function(sp,y)
{
  # browser()
  if(missing(y)) {
    dat=data.frame(x=I(sp))
  }
    else
    {
      dat=data.frame(x=I(sp),y)
    }
  class(dat$x)="matrix"
  # colnames(dat)=c("x","y")
  if(exists("nam")) rownames(dat)=nam   # rownames(dat$x)=nam
  if(exists("lo")) colnames(dat$x)=lo

  return(dat)
}
