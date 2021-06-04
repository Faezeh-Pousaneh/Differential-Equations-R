

require(bvpSolve)
s=0.005
while_loop_num=1
Sol5<-matrix(nrow=((100-5)/5),ncol=5)
while ( s<0.05)  {
	while_loop_num=while_loop_num+1
	T<- 1 + s
	s=s*10
	J<-  (1./6.)

	p<-   0.001
	lambda<- 10
	b<- 1
	alfa<-  (1/((lambda^2)*p)) *b

	h0<- 0.001/6
	hL<- h0

	C0<-  -0.001
	CL<-  C0

	y10<-   -0.155
	y30<-   -0.0002
	y40<-   0.0002

	A <- p*p/((1-2*p)^2)

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


	
	loop_step=0
	z=0
	for(z in seq(5,100,5)) {
		L<-z
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
		for (i in 0:ncol(Sol))    Sol1[,i]<-Sol[,i]

		Sol1[,ncol(Sol)+1]=   Sol[,2]- ( (((2*A-(-A*(tanh(Sol[,4]/b))^2+A)^0.5)*(1+Sol[,2]))/(4*A+(tanh(Sol[,4]/b))^2-1)) ) +p
		Sol1[,ncol(Sol)+2]=  -  ( (((2*A-(-A*(tanh(Sol[,4]/b))^2+A)^0.5)*(1+Sol[,2]))/(4*A+(tanh(Sol[,4]/b))^2-1)) ) * tanh(Sol[,4]/b) 

		Sol1[,ncol(Sol)+3]=  0.5*J* ( -6*(Sol1[,2])*(Sol1[,2])+(Sol1[,3])*(Sol1[,3]))   -
							 0.5*(T/(alfa))*  (Sol1[,5]*Sol1[,5])+T*( Sol1[,7] * Sol1[,4])+
							 0+
							 T*(0.5*(1-2*(p-Sol1[,6]+Sol1[,2])+Sol1[,2])*log((0.5*(1-2*(p-Sol1[,6]+Sol1[,2])+Sol1[,2])) , base = exp(1))+
							 0.5*(1-Sol1[,2])*log((0.5*(1-Sol1[,2])) , base = exp(1))+
							  0.5*((p-Sol1[,6]+Sol1[,2]+Sol1[,7]))*log((0.5*((p-Sol1[,6]+Sol1[,2]+Sol1[,7]))) , base = exp(1))+
							 0.5*((p-Sol1[,6]+Sol1[,2]-Sol1[,7]))*log((0.5*((p-Sol1[,6]+Sol1[,2]-Sol1[,7]))) , base = exp(1))+
							 0
							 -(0.5*(1-2*p)*log((0.5*(1-2*p)) , base = exp(1))+
							 0.5*(1)*log((0.5*(1)) , base = exp(1))+
							 0.5*(p)*log((0.5*(p)) , base = exp(1))+
							 0.5*(p)*log((0.5*(p)) , base = exp(1))))-T*log((p/((1-2*p))), base = exp(1))*(-Sol1[,6]+Sol1[,2])-0.5*T*(log((1-2*p),base = exp(1)))*(Sol1[,2])
					
						
		Sol2<-matrix(nrow=1,ncol=1)
		Sol2[1,1]= 0		   
		for (j in 1:nrow(Sol1))  {
			Sol2[1,1]<- 0.5*sum((Sol1[j,1]-Sol1[j-1,1]) *  ( Sol1[j-1,8] +Sol1[j,8]))+
			Sol2[1,1] 
					 }
					 
					 
		Sol3<-matrix(nrow=1,ncol=1)
		Sol3[1,1]= 0.5*J*(Sol1[1,2])*(Sol1[1,2])-(h0)*(Sol1[1,2])+C0*Sol1[1,4]+
		0.5*J*(Sol1[nrow(Sol1),2])*(Sol1[nrow(Sol1),2])-(hL)*(Sol1[nrow(Sol1),2])+CL*Sol1[nrow(Sol1),4]
				   
		Sol4<-matrix(nrow=1,ncol=1)
		Sol4[1,1]= Sol3[1,1]+Sol2[1,1]
		Sol5[loop_step,1]= L
		Sol5[loop_step,while_loop_num]= Sol4[1,1]
		loop_step=loop_step+1
	}		 
	
	Sol5[,while_loop_num]=Sol5[,while_loop_num]-Sol5[(100-5)/5,while_loop_num]
}	

matplot((Sol5[,1]), (Sol5[,2]-Sol5[(100-5)/5,2]), type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z",ylab="Potential")
	#for (i in 0:ncol(Sol1))
		write.table(Sol5 ,file="full.txt")




#par(mfrow=c(1,3))
#matplot((Sol1[,1]), Sol1[,2], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="solvent+solute")
#matplot((Sol1[,1]), Sol1[,6], type = "l", lty = 1,  col = c("red"),lwd=2,font=5, xlab = "z", ylab="Solvent")
#matplot((Sol1[,1]), Sol1[,7], type = "l", lty = 1,  col = c("blue"),lwd=2,font=5, xlab = "z", ylab = "Charge density")


#plot(Sol1, col = c("red", "blue", "green"),lwd=2,font=5, xlab = "z", mfrow = c(1,3),which = c("Solvent","Solute","Psi")
#)

legend("top", lty = 1:2, title = "T/Tc=1.006,Rho=0.00108, h0=-0.01,hL= -0.01 (nm)^-3,Rn=1,Rs=1, charge density in units of (nm)^-2 is =" ,
legend = c(0.00035))
