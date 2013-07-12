function der,y
;
  common grid,rr,rrr,dr,dr1,dr2
  common dim,nx,mx,l1,l2
  common order,nghost,order
;
  if (order eq 6) then begin

      dy=dr1/60.*(   45*(y[l1+1:l2+1] - y[l1-1:l2-1])$
                     -9*(y[l1+2:l2+2] - y[l1-2:l2-2])$
                      + (y[l1+3:l2+3] - y[l1-3:l2-3]))
  endif else begin
      dy=.5*dr1*(y[2:mx-1]-y[0:mx-3])
  endelse
;
return,dy
;
end
