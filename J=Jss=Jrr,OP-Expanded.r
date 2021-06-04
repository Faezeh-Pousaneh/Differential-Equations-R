## Mathematics Series, SIAM, Philadelphia, 1995].
## g'' = (g f' - f g')/eps
## f'''' = (-ff'''-gg')/eps
## g(0)=-1,f(0)=0,f'(0)=0, g(1)=1,f(1)=0,f'(1)=0
##
## 1st order system (y1=f, y3=h, y5=w)
## y1' = y2

## =============================================================================

# init <- c (NA, -1, NA,-1,NA,-0.02)
# yend   <- c(NA, 1, NA,1,NA,0.02)

# res_init  <- function (Y,yini,pars)  
# with (as.list(Y), {
	# f7<-Y[1]-1
	# f8<-Y[3]-1
	# f9<- -0.02
	# return(list(c(NA,f7,NA, f8,NA, f9)))
# })
 # res_end  <- function (Y,yini,pars)  with (as.list(Y), {
	# f10<--Y[1]+1
	# f11<--Y[3]+1
	# f12<- 0.02
# return(list(c(NA,f10,NA, f11,NA, f12)))
# })
require(bvpSolve)

T<- 1.002
J<-  (1./6.)

	p<-   0.0080
	alfa<-  (1/((9.25*p)^0.5)) * T
	b<-  T

	h0<- 0.001
	hL<- 0.001

	C0<-  0.001
	CL<-  0.001

	y10<-   -0.155
	y30<-   -0.0002
	y40<-   0.0002

A <- p*p/((1-2*p)^2)

L <-20
equ1<-function(t, Y,pars) {
 return(list(c  (

   
		  dy1  <-    Y[2],
		  dy2  <-   -6*Y[1]  +  T* (1/(J)) *Y[1]-  T* (1/(2*J))* (p/(b*b))  *(Y[3]/(1-Y[1]))^2 ,       
		 
		  
		  dy3  <-   Y[4],
		  dy4  <-   alfa*   (p/b)*(Y[3]/(1-Y[1])) 
		 
		)))
}


bound <- function(i, Y, pars) {
  with (as.list(Y), {
    if (i == 1) return ((Y[2])-Y[1]    +(h0/J))
    if (i == 2) return ((Y[4])         +alfa*C0)
    if (i == 3) return ((Y[2])+Y[1]    -(hL/J))
    if (i == 4) return (Y[4]           -alfa*CL)      
 })
}
L


x       <- seq(0, L,0.2)


xguess <- seq(0,L,len = 4)
yguess <- matrix(nc = 4,data = rep(c(y10,    (y10-(h0/J)),   y30 ,y40),4))
rownames(yguess) <- c("Solvent","y2","Psi", "y4")




print(system.time(
Sol <- bvptwp(yini = NULL, x = x, func = equ1, bound = bound,
              xguess = xguess, yguess = yguess, leftbc = 2,
              atol = 1e-10)			  
			 ))		  


Sol[,1]=Sol[,1]/L

for (i in 0:ncol(Sol))
write.table(Sol ,file="(-,-)compnumeric1.txt")

par(mfrow=c(1,1))
matplot((Sol[,1]), Sol[,2], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="solvent+solute")


#plot(Sol1, col = c("red", "blue", "green"),lwd=2,font=5, xlab = "z", mfrow = c(1,3),which = c("Solvent","Solute","Psi")
#)

legend("top", lty = 1:2, title = "T/Tc=1.02, h=0.02 (nm)^-3, charge density in units of (nm)^-2 is =" ,
legend = c(0.06))
