function bisect,left,right,sigma,omega,mdot,alpha
;
    common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
;
    maxit=100
    it=0
;
    while (abs(right-left) gt epsi) do begin
;
        midpoint=.5*(right+left)
;
        func=get_phi(left    ,sigma,omega,mdot,alpha) & fleft=func.f
        func=get_phi(midpoint,sigma,omega,mdot,alpha) & fmid =func.f
;
; initialize the bisecting
;
        if (fleft*fmid lt 0) then begin
            right=midpoint
        endif else begin
            left=midpoint
        endelse
        it=it+1
        if (it gt maxit) then begin
            print,'exceed maximum of iterations'
            print,'left,right=',left,right
            stop
        end
    endwhile
;
    out=.5*(right+left)
;
    return,out
;
end
