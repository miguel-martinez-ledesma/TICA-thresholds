% lm_plots ( t, y_dat, y_fit, sigma_y, cvg_hst, filename )
% Plot statistics of the results of a Levenberg-Marquardt least squares
% analysis with lm.m

function lm_plots ( t, y_dat, y_fit, sigma_y, cvg_hst, filename )

y_dat = y_dat(:);
y_fit = y_fit(:);

[max_it,n] = size(cvg_hst); n = n-3;

figure(101); % ---------- plot convergence history of parameters, chi^2, lambda
 clf
 subplot(211)
  plot( cvg_hst(:,1), cvg_hst(:,2:n+1), '-o','linewidth',4);
  for i=1:n
   text(1.02*cvg_hst(max_it,1),cvg_hst(max_it,1+i), sprintf('%d',i) );
  end
   ylabel('parameter values')
 subplot(212)
  semilogy( cvg_hst(:,1) , [ cvg_hst(:,n+2) cvg_hst(:,n+3)], '-o','linewidth',4)
   text(0.8*cvg_hst(max_it,1),0.2*cvg_hst(max_it,n+2), '\chi^2/(m-n+1)' );
   text(0.8*cvg_hst(max_it,1),0.5*cvg_hst(max_it,n+3), '\lambda' );
%  legend('\chi^2/(m-n+1)}', '\lambda', 3);
   ylabel('\chi^2/(m-n+1) and \lambda')
   xlabel('function calls')


figure(102); % ------------ plot data, fit, and confidence interval of fit
 clf
   plot(t,y_dat,'og', t,y_fit,'-b',...
        t,y_fit+1.96*sigma_y,'.k', t,y_fit-1.96*sigma_y,'.k');
    legend('y_{data}','y_{fit}','95% c.i.','',0);
    ylabel('y(t)')
    xlabel('t')
% subplot(212)
%  semilogy(t,sigma_y,'-r','linewidth',4);
%    ylabel('\sigma_y(t)')

figure(103); % ------------ plot histogram of residuals, are they Gaussean?
 clf
 hist(real(y_dat - y_fit),20)
  title('histogram of residuals')
  axis('tight')
  xlabel('y_{data} - y_{fit}')
  ylabel('count')
