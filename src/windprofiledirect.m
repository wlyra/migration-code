function [sigmaw,mdotR]=windprofiledirect(R,Rin)
%computes the direct wind profile for a given grid R with hole radius Rin


%Sigma_w is the dimensionless wind profile
%mdot is the value of mass-loss rate 

 x_in=R-Rin;
 
 
 	% coefficients
	a = -4.3822599312953403E-01;
	b = -1.0658387115013401E-01;
	c = 5.6994636363474982E-01;
	d = 1.0732277075017336E-02;
	e = -1.3180959703632333E-01;
	f = -1.3228570869396541E+00;

	temp = (a * exp(b*x_in) + c * exp(d*x_in) + e * exp(f*x_in));

	y = temp;
 
    y(x_in<0)=0;
    
    store=gradient(y,R);
    sigmaw=store./R.*exp(-((R-Rin)/57).^10);
    sigmaw=sigmaw/trapz(R,R.*sigmaw);
 
mdotR=cumtrapz(R,R.*sigmaw);

