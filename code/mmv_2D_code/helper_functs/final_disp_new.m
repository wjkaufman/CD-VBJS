function final_disp_new(out,opts)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

fprintf('Number of total iterations is %d. \n',out.total_iter);
fprintf('Termination data error: %5.3f\n',out.final_error);
fprintf('Final ||w||_1: %5.3g\n',out.final_wl1);
fprintf('Final ||Du-w||^2: %5.3g\n',out.final_Du_w);
if opts.phimin.joint
    fprintf('Final ||phi||^2: %5.3g\n',out.lam7(end));
    fprintf('Final ||Tphi||^2: %5.3g\n',out.lam6(end))
    if out.true_phase
        fprintf('Final ||phi-phi_true||_2: %5.3g\n',out.PE_err(end))
    end
end

if opts.disp

    if opts.phimin.joint
       for k = 1:opts.outer_iter  
           rel_phi((k-1)*opts.inner_iter+1:k*opts.inner_iter) = out.rel_chg_phi(k); 
       end
    end
    
    figure;
    if opts.phimin.joint
        S = suptitle('With Phase Estimation');
        set(S,'FontSize',16,'interpreter','latex');
    else
        S = suptitle('Without Phase Estimation');
        set(S,'FontSize',16,'interpreter','latex');
    end
    subplot(2,1,1); set(gca,'FontSize',14);
    plot(out.lam1,'LineWidth',1.5);  %plot lam1, ||W||_1
    hold on;
    plot(out.lam3.*out.mu,'LineWidth',1.5);  %plot lam3, mu||Au -f||^2
    plot(abs(out.f),'LineWidth',1.5);   %plot f, the objective function
    plot(2:opts.inner_iter:max(size(out.f)),...
        out.f(2:opts.inner_iter:end),'kx','Linewidth',2);
    L = legend('$||w||_1$','$\mu||\mathcal{F}f - \hat{f}||_2^2$',...
        'obj function','mlp update');
    set(L,'Location','NorthEast','interpreter','latex','FontSize',12);
%     h = xlabel('iteration');
%     set(h,'interpreter','latex','FontSize',18)
    hold off;
    
    subplot(2,1,2); set(gca,'FontSize',14);
    plot(out.rel_error,'LineWidth',1.5);
    hold on
    plot(out.rel_chg_inn,'LineWidth',1.5);
    plot(out.rel_lam2,'LineWidth',1.5);
    if opts.phimin.joint
        plot(rel_phi,'LineWidth',1.5);
        plot(2:opts.inner_iter:max(size(out.f)),...
            rel_phi(2:opts.inner_iter:end),'kx','Linewidth',2);
        L = legend('Rel error','Rel chg $f$','Rel chg $||\mathcal{L}f-w||_2$','Rel chg $\phi$','mlp update');
        set(L,'Location','NorthEast','interpreter','latex','FontSize',12);
    else
        plot(2:opts.inner_iter:max(size(out.f)),...
            0,'kx','Linewidth',2);
        L = legend('Rel error','Rel chg $f$','Rel chg $||\mathcal{L}f-w||_2$','mlp update');
        set(L,'Location','NorthEast','interpreter','latex','FontSize',12);
    end
    h = xlabel('iteration');
    set(h,'interpreter','latex','FontSize',16)
    hold off

end