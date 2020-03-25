%% This script generates the table showing the ML estimation results

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
aux(1).SRLR = [6, 1, 4, 8];    %Indices of the SR metrics of interest

%Long-run
aux(2).SRLR = [3,7];           %Indices of the LR metrics of interest

%% Create labels for the short-run and long-run metrics
% SR0, SR1, SR2, SR3, LR1, LR2

Labels  = cell(2*(size(aux(1).SRLR,2)+size(aux(2).SRLR,2)),1);

for i_m = 1:size(Labels,1)/2
    
    if i_m <= size(aux(1).SRLR,2)
    
    if i_m == 1
        
    Labels((2*i_m)-1) = {strcat('Success Rate (SR',num2str(i_m-1),')')};     
        
    else    
        
    Labels((2*i_m)-1) = {strcat('Short-Run Metric \#',num2str(i_m-1))};
                        
    %'(SR',num2str(i_m-1),')')};
        
    end
    
    else
        
    Labels((2*i_m)-1) = {strcat('Long-run metric \#',num2str(i_m-4))};
                         
    %'(LR',num2str(i_m-4),')')};
     
    end
    
end

%% Extract the unweighted estimator of degrees of freedom and standard errors

indexSR  = aux(1).SRLR;

indexLR  = aux(2).SRLR;

M  = cell(2*(size(aux(1).SRLR,2)+size(aux(2).SRLR,2)),1);   %"Mean parameter"

s  = cell(2*(size(aux(1).SRLR,2)+size(aux(2).SRLR,2)),1);   %"Standard deviation"

v  = cell(2*(size(aux(1).SRLR,2)+size(aux(2).SRLR,2)),1);   %"Degrees of freedom"

for i_m = 1:size(Labels,1)/2
   
    if i_m <= size(aux(1).SRLR,2)                           %Indicator for SR metrics
    
    Maux         =  mle_output.beta_no_weights(indexSR(1,i_m),1); 
    
    Maux_std     =  mle_output.variance_matrix_no_weights{indexSR(1,i_m)}(1,1)^.5;
        
    M((2*i_m))   = {strcat('(', num2str(Maux_std,3),...
                          ')')};
    
    M((2*i_m)-1) = {num2str(Maux,3)};      
                        
    clear Maux Maux_std 
    
    saux         = mle_output.beta_no_weights(indexSR(1,i_m),2);
    
    saux_std     = mle_output.variance_matrix_no_weights{indexSR(1,i_m)}(2,2)^.5;
    
    s((2*i_m))   = {strcat('(', num2str(saux_std,3),...
                          ')')};
    
     s((2*i_m)-1) = {num2str(saux,3)};                                                        
                                          
    clear saux saux_std 
    
    vaux         = mle_output.beta_no_weights(indexSR(1,i_m),3);
    
    vaux_std     = mle_output.variance_matrix_no_weights{indexSR(1,i_m)}(end,end)^.5;
    
    tv           = (3-vaux)/vaux_std;
                      
    v((2*i_m)) = {strcat('(', num2str(vaux_std,3),...
                          ')')}; 
    
                                                            %Third col is
                                                            %the estimator
                                                            %of v
        if tv > norminv(.999)
                                                            
        v((2*i_m)-1)   = {strcat(num2str(vaux,3),...
                          '**')}; 
                      
        elseif tv > norminv(.990)
            
        v((2*i_m)-1)   = {strcat(num2str(vaux,3),...
                          '*')}; 
        else
        
        v((2*i_m)-1)   = {num2str(vaux,3)};    
                      
        end
    
    clear vaux vaux_std tv
        
    else
        
    Maux         =  mle_output.beta_no_weights(indexLR(1,i_m-size(aux(1).SRLR,2)),1); 
    
    Maux_std     =  mle_output.variance_matrix_no_weights{i_m-size(aux(1).SRLR,2)}(1,1)^.5;
    
    M((2*i_m)) = {strcat('(',num2str(Maux_std,3),')')}; 
        
    M((2*i_m)-1)   = {num2str(Maux,3)};      
                      
    clear Maux Maux_std 
    
    saux         = mle_output.beta_no_weights(indexLR(1,i_m-size(aux(1).SRLR,2)),2);
    
    saux_std     = mle_output.variance_matrix_no_weights{i_m-size(aux(1).SRLR,2)}(2,2)^.5;
                      
    s((2*i_m)) = {strcat('(',num2str(saux_std,3),')')};
    
    s((2*i_m)-1)   = {num2str(saux,3)};
                      
    clear saux saux_std ts  
    
    vaux         = mle_output.beta_no_weights(indexLR(1,i_m-size(aux(1).SRLR,2)),3);
    
    vaux_std     = mle_output.variance_matrix_no_weights{i_m-size(aux(1).SRLR,2)}(end,end)^.5;
    
    tv           = (vaux-3)/vaux_std;
                      
    v((2*i_m)) = {strcat('(',num2str(vaux_std,3),')')};
    
    v((2*i_m)-1)   = {strcat('(', num2str(vaux,3),...
                          ')')};
                      
    if tv > norminv(.999)
                                                            
    v((2*i_m)-1)   = {strcat(num2str(vaux,3),...
                          '**')}; 
                      
    elseif tv > norminv(.990)
            
    v((2*i_m)-1)   = {strcat(num2str(vaux,3),...
                          '*')}; 
    else
        
    v((2*i_m)-1)   = {num2str(vaux,3)};    
                      
    end
    
    clear vaux vaux_std tv
                      
                      
    end
    
end

%% 

addpath './matlab/functions'

Table1 = table(Labels,v,M,s);

input.tableColumnAlignment      = 'c';

input.data                      = Table1;

input.dataFormat                = {'%i'};
    
input.tableBorders              = 0;

input.transposeTable            = 0;

input.makeCompleteLatexDocument = 1;

input.tableCaption              = '';

input.tableLabel                = '';

Table1tex                       = latexTable(input);

delete('./output/tables/mle-results.tex');

fid       =fopen('./output/tables/mle-results.tex','w');

[nrows,~] = size(Table1tex);

fprintf(fid,'%s\n','\begin{tabular}{lccc}');

fprintf(fid,'%s\n','\toprule');

fprintf(fid,'%s\n',strcat('&','$\alpha$', '&', '$M$','&','$s$','\\'));

fprintf(fid,'%s\n','\midrule');

for row = 7:nrows-5
    
    fprintf(fid,'%s\n',Table1tex{row,:});
    
end
fprintf(fid,'%s\n','\bottomrule');

fprintf(fid,'%s\n',Table1tex{nrows-4,:});

fclose(fid);
