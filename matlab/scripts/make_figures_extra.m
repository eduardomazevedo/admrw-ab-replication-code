%% This script generates additional figures for the appendix of the paper

%% Load and prepare data for the figure comparing weighted/unweighted ML

load('./matlab/mat/mle.mat')

mle_output = mle_output(strcmp(mle_output.algorithm, 'gradient knitro'), :);
           % Select the results corresponding to the gradient knitro
           % procedure
           
mle_output_weighted_unweighted ...
           = mle_output(mle_output.use_weights==1,[3,4,5]);
           % Select all the results corresponding to the weighted ML procedure 

beta_no_weights ...
           = mle_output(mle_output.use_weights==0,4); 
           % Select the betas corresponding to the unweighted ML
       
beta_no_weights.Properties.VariableNames ...
           = {'beta_no_weights'}; 
           
variance_matrix_no_weights ...
           = mle_output(mle_output.use_weights==0,5);
           % Select the variances corresponding to the unweighted ML
           
variance_matrix_no_weights.Properties.VariableNames ...
           = {'variance_matrix_no_weights'};        
       
mle_output_weighted_unweighted = [mle_output_weighted_unweighted,...
                                  beta_no_weights,...
                                  variance_matrix_no_weights]; 
                              
clear beta_no_weights variance_matrix_no_weights 

sortedNames = sort(mle_output_weighted_unweighted.Properties.VariableNames(2:end));

clear mle_output

mle_output  = [mle_output_weighted_unweighted(:,1),...
               mle_output_weighted_unweighted(:,sortedNames)];

clear sortedNames mle_output_weighted_unweighted          

%% Define the vector of indices for the short-run and long-run

%Short-run
aux(1).SRLR = [6, 1, 4, 8];

%Long-run
aux(2).SRLR = [3,7];

%% ML estimator of the degrees of freedom with standard errors

confidence     = .95;                        %Confidence level

for i_SRLR     = 1:2                         %loop over SR and LR metrics
   
    n_mets     = size(aux(i_SRLR).SRLR,2);   %number of SR or LR metrics
    
    for i_mets = 1:n_mets      %loop over total number of SR 
                                             %or LR metrics
       
        index  = aux(i_SRLR).SRLR(1,i_mets);
                    
        if i_SRLR ...
               == 1                             %Generates labels SR#
        
         
        ML(i_SRLR).lab(1,i_mets) ...
               =  {strcat('SR',num2str(i_mets-1))};
  
        else                                      %Generates labels LR#
      
        ML(i_SRLR).lab(1,i_mets) ...
               = {strcat('LR',num2str(i_mets))};
  
        end
        
    ML(i_SRLR).v(i_mets,1)  ...
               = mle_output.beta_no_weights(index,3);                          
                                            %Collect point estimators
                                            
    ML(i_SRLR).v(i_mets,2)  ...  
               = mle_output.variance_matrix_no_weights{index}(end,end)^.5;
                                            %Collect standard errors
                                            
    ML(i_SRLR).vweighted(i_mets,1) ...
               = mle_output.beta(index,3);
           
    ML(i_SRLR).vweighted(i_mets,2) ...
               = mle_output.variance_matrix{index}(end,end)^.5;
  
    clear index
        
    end
  
    clear n_mets   
    
end

%% 4) Plot ML results and standard errors

figure()

c            = norminv(1-(1-confidence)/2,0,1);

for i_metric = 1: size(ML(1).v,1)              %Plot lines for the CI's SR
    
    plot([i_metric i_metric], [ ML(1).v(i_metric,1)-...
                               c*ML(1).v(i_metric,2),...
                               ML(1).v(i_metric,1)-eps],'-b'); hold on
                           
    plot([i_metric+.1 i_metric+.1], [ ML(1).vweighted(i_metric,1)-...
                               c*ML(1).vweighted(i_metric,2),...
                               ML(1).vweighted(i_metric,1)-eps],'-r'); hold on
                           
    
    plot([i_metric i_metric], [ML(1).v(i_metric,1)+eps,...
                               ML(1).v(i_metric,1)+...
                               c*ML(1).v(i_metric,2)],'-b'); hold on 
                           
    plot([i_metric+.1 i_metric+.1], [ML(1).vweighted(i_metric,1)+eps,...
                               ML(1).vweighted(i_metric,1)+...
                               c*ML(1).vweighted(i_metric,2)],'-r'); hold on 
                           
end    

for i_metric2 = 1: size(ML(2).v,1)            %Plot lines for the CI's LR
    
    plot([i_metric+i_metric2 i_metric+i_metric2],...
                              [ ML(2).v(i_metric2,1)-...
                               c*ML(2).v(i_metric2,2),...
                               ML(2).v(i_metric2,1)-eps],'-b'); hold on
                           
    plot([i_metric+i_metric2+.1 i_metric+i_metric2+.1],...
                              [ ML(2).vweighted(i_metric2,1)-...
                               c*ML(2).vweighted(i_metric2,2),...
                               ML(2).vweighted(i_metric2,1)-eps],'-r'); hold on
                           
   plot([i_metric+i_metric2 i_metric+i_metric2],...
                              [ ML(2).v(i_metric2,1)+eps,...
                               ML(2).v(i_metric2,1)+...
                               c*ML(2).v(i_metric2,2)],'-b'); hold on 
                           
   plot([i_metric+i_metric2+.1 i_metric+i_metric2+.1],...
                              [ ML(2).vweighted(i_metric2,1)+eps,...
                               ML(2).vweighted(i_metric2,1)+...
                               c*ML(2).vweighted(i_metric2,2)],'-r'); hold on                         
                           
end  

plot(1:1:6,[ML(1).v(:,1);ML(2).v(:,1)],'o','MarkerFaceColor','b'); hold on
                                            %Plot point estimators
                                            
plot([1:1:6]+.1,[ML(1).vweighted(:,1);ML(2).vweighted(:,1)],'o','MarkerFaceColor','r'); hold on
                                            %Plot point estimators                                            
                                            
plot([0 size(ML(1).v,1)+size(ML(2).v,1)+1],[1 1],':b'); hold on
                                            % v=1

plot([0 size(ML(1).v,1)+size(ML(2).v,1)+1],[3 3],':b'); hold off
                                            % v=3

strValues1 = strtrim(cellstr(num2str(ML(1).v(:,1),3)));

text((1:size(ML(1).v,1))'-.4,ML(1).v(:,1),strValues1,'VerticalAlignment','bottom');

strValues2 = strtrim(cellstr(num2str(ML(2).v(:,1),3)));

text((size(ML(1).v,1)+1:size(ML(1).v,1)+size(ML(2).v,1))'-.4,...
      ML(2).v(:,1),strValues2,'VerticalAlignment','bottom');


xlim([0 size(ML(1).v,1)+size(ML(2).v,1)+1])

ylim([0 6])

set(gca,'XTick',1:6,'XTickLabel',horzcat(ML(1).lab,ML(2).lab))

xlabel('Metric')

ylabel('Tail Coefficient (\alpha)')

%% Save output

print(gcf,'-depsc2','./output/figures/mle-weighted-comparison.eps')

%% Preparation for the figure containing marginal rates of substitution

%% Positive part of the normal prior

f1     = @(s,m) (m.*(1-normcdf(-m./s))) + (s.*normpdf(-m./s));

[s,m]  = meshgrid(linspace(.001,.005,1000),linspace(-.005,0,1000));

f1val  = f1(s,m);

%% Positive part of the t-prior

f2    = @(s,m,alpha) (m.*tcdf(-m./s,alpha,'upper')) +...
                     (s.* (gamma((alpha-1)/2)/(2*gamma(alpha/2).*gamma(1/2))).*...
                     (alpha.^(alpha/2)).*((alpha + ((-m./s).^2)).^((1-alpha)/2))); 
                     %Formula derived in notes

f2val = f2(s,m,mle_output.beta_no_weights(6,3)); %alpha = ML for SSR in the paper

%% Isoquant curve for the normal

figure()

contour(s,m,f1val,20);

ax = gca;

ax.FontSize = 16; 

xlabel('Scale s')

ylabel('Mean M')

%% Save figure

print(gcf,'-depsc2','./output/figures/isoquant-normal.eps')

%% Isoquant curve for the t

figure()

contour(s,m,f2val,20);

ax = gca;

ax.FontSize = 16; 

xlabel('Scale s')

ylabel('Mean M')

%% Save figure

print(gcf,'-depsc2','./output/figures/isoquant-t.eps')
