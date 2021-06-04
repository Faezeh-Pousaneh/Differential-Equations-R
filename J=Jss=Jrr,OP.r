
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
  dy2  <-   -6*Y[1]  +  T*( (1/(2*J)) * log(  (((1- 2*(((2*A-(-A*(tanh(Y[3]/b))^2+A)^0.5)*(1+Y[1]))/(4*A+(tanh(Y[3]/b))^2-1))+Y[1])/((1-Y[1])))) , base = exp(1)) 
			- (1/(2*J)) * log ((1-2*p) , base = exp(1))) ,         
 
  
  dy3  <-   Y[4],
  dy4  <-   alfa*   (((2*A-(-A*(tanh(Y[3]/b))^2+A)^0.5)*(1+Y[1]))/(4*A+(tanh(Y[3]/b))^2-1))   * tanh(Y[3]/b) 
 
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



Sol1<-matrix(nrow=nrow(Sol),ncol=ncol(Sol)+3)
for (i in 0:ncol(Sol))    Sol1[,i]<-Sol[,i]/L

Sol1[,ncol(Sol)+1]=   Sol[,2]- ( (((2*A-(-A*(tanh(Sol[,4]/b))^2+A)^0.5)*(1+Sol[,2]))/(4*A+(tanh(Sol[,4]/b))^2-1)) ) +p
Sol1[,ncol(Sol)+2]=  -  ( (((2*A-(-A*(tanh(Sol[,4]/b))^2+A)^0.5)*(1+Sol[,2]))/(4*A+(tanh(Sol[,4]/b))^2-1)) ) * tanh(Sol[,4]/b) 


for (i in 0:ncol(Sol1))
write.table(Sol1 ,file="(+,+),0.00002.txt")

par(mfrow=c(1,1))
matplot((Sol1[,1]), Sol1[,2], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="solvent+solute")
#matplot((Sol1[,1]), Sol1[,4], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="Phi.el")
#matplot((Sol1[,1]), Sol1[,6], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="Solvent")
#matplot((Sol1[,1]), Sol1[,7], type = "l", lty = 1,  col = c("blue"),lwd=2,font=5, xlab = "z", ylab = "Charge density")



legend("top", lty = 1:2, title = "T/Tc=1.02, h=0.02 (nm)^-3, charge density in units of (nm)^-2 is =" ,
legend = c(0.06))
