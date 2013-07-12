function der2,y
;
  common grid,rr,rrr,dr,dr1,dr2
  common dim,nx,mx,l1,l2
  common order,nghost,order
;
  if (order eq 6) then begin
      d2y =dr2/180.*(-490.* y[l1:l2]                    $
                     +270.*(y[l1+1:l2+1] + y[l1-1:l2-1])$
                     - 27.*(y[l1+2:l2+2] + y[l1-2:l2-2])$
                     +  2.*(y[l1+3:l2+3] + y[l1-3:l2-3]))
  endif else begin
      d2y=.5*dr2*(y[2:mx-1]-2*y[1:mx-2]+y[0:mx-3])  
  endelse
;
return,d2y
end
