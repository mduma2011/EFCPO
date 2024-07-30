function displaytsp(X, Y, Distance,TAU, ETA, nE, ant, tour)
 
  XList = zeros(nE,2);
  YList = zeros(nE,2);
  
  ant0 = [ant; ant(1)]; 
  figure;
  hold on; 
  title(sprintf('Best Tour = %d',tour));
  plot(X, Y,'-O','MarkerFaceColor','b','LineStyle','none');
   
  for i = 1:nE
      
      x1 = X(ant0(i)); x2 = X(ant0(i+1));
      y1 = Y(ant0(i)); y2 = Y(ant0(i+1));
      XList(i,1) = x1; XList(i,2) = x2;
      YList(i,1) = y1; YList(i,2) = y2;
      
      dx = round(( x1 + x2)/2, 2);
      dy = round(( y1 + y2)/2, 2);
      dist = Distance(ant0(i), ant0(i+1));
      tau = TAU(ant0(i), ant0(i+1));     
      strcatf = sprintf('%d:%0.3e (%d)',dist, tau, i);
      if x1 < x2
         text(dx+1,dy, strcatf); 
      else
         text(dx-1,dy, strcatf);  
      end      
  end 
  plot(XList,YList,'-','MarkerFaceColor','r');
  hold off;
  
end

