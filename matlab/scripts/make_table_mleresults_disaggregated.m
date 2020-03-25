%% This script generates the table showing the ML estimation results, disaggregated by budget area

load('./matlab/mat/mle_disaggregated.mat')

%% Define the vector of indices for the budget areas

aux(1).IBA = [4, 5, 6];    %Indices of the budget areas

%% Create labels for the budget areas

Labels  = cell(2*(size(aux(1).IBA,2)),1);

Labels(1) = {'Ux'};

Labels(3) = {'Relevance'};

Labels(5) = {'Engagement'}; 
                        
%% Extract the estimator of degrees of freedom and standard errors

indexBA  = aux(1).IBA;

M  = cell(2*(size(aux(1).IBA,2)),1);   %"Mean parameter"

s  = cell(2*(size(aux(1).IBA,2)),1);   %"Standard deviation"

v  = cell(2*(size(aux(1).IBA,2)),1);   %"Degrees of freedom"

for i_m = 1:size(Labels,1)/2
   
    if i_m <= size(aux(1).IBA,2)                           %Indicator for SR metrics
    
    Maux         =  mle_output_disaggregated.beta(indexBA(1,i_m),1); 
    
    Maux_std     =  mle_output_disaggregated.variance_matrix{indexBA(1,i_m)}(1,1)^.5;
        
    M((2*i_m))   = {strcat('(', num2str(Maux_std,3),...
                          ')')};
    
    M((2*i_m)-1) = {num2str(Maux,3)};      
                        
    clear Maux Maux_std 
    
    saux         = mle_output_disaggregated.beta(indexBA(1,i_m),2);
    
    saux_std     = mle_output_disaggregated.variance_matrix{indexBA(1,i_m)}(2,2)^.5;
    
    s((2*i_m))   = {strcat('(', num2str(saux_std,3),...
                          ')')};
    
     s((2*i_m)-1) = {num2str(saux,3)};                                                        
                                          
    clear saux saux_std 
    
    vaux         = mle_output_disaggregated.beta(indexBA(1,i_m),3);
    
    vaux_std     = mle_output_disaggregated.variance_matrix{indexBA(1,i_m)}(end,end)^.5;
    
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

delete('./output/tables/mle-results_disaggregated.tex');

fid       =fopen('./output/tables/mle-results_disaggregated.tex','w');

[nrows,~] = size(Table1tex);

fprintf(fid,'%s\n','\begin{tabular}{lccc}');

fprintf(fid,'%s\n','\toprule');

fprintf(fid,'%s\n',strcat('Budget Area','&','$\alpha$', '&', '$M$','&','$s$','\\'));

fprintf(fid,'%s\n','\midrule');

for row = 7:nrows-5
    
    fprintf(fid,'%s\n',Table1tex{row,:});
    
end
fprintf(fid,'%s\n','\bottomrule');

fprintf(fid,'%s\n',Table1tex{nrows-4,:});

fclose(fid);
