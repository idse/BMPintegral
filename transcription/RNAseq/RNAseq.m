scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
cd(scriptPath);

% FKPM
% fragments per megabase per million reads
% normalized by gene length: not ideal

% normalized to what?
A = readtable('GSE77057_DESeq2_normalized_count_matrix.txt');

% sum of rows = fixed -> normalized across samples

% also check out single cell atlas
% https://www.ebi.ac.uk/gxa/sc/experiments
% and
% https://en.wikipedia.org/wiki/BioMart

% GTF = gene transfer format: count reads per gene (after mapping)
% Gencode V18 used 

% TPM = transcript per million
% count divided by gene length divided  by million read

% FPKM = 

% SAM = sequence alignment map
% sequence, where it mapped to, quality
% BAM = binary SAM

% paired end -> two fastq files

% trim_galore: detect adapter
% install cutadapt

% trim_galore --Paired on fastq files produces another fastq file trim_v1, trim_v2
% bioconda
 
% STAR make genome and transcript index
% overhang 124 for this specific case, set to 99 in general
% will take a few hours
% output .BWT ?

% STAR on raw data
% running mapping jobs
% - genomeLoad LoadAndKeep
% - outSAMtype BAM SortedByCoordinate
% -readFilesCommand zcat or gunzip
% Encode options: use

% HTSeq (Python alternative to featureCounts, which is R)

% from gencode FASTA file for genome sequence GRCh38.p13

%% write gene ID to file for external conversion (one time event)

ensemblID = table2cell(A(:,1));

% strip off .xx because that is the version of actual genomic sequence
for i = 1:numel(ensemblID)
    s = strsplit(ensemblID{i},'.');
    ensemblID{i} = s{1};
end
    
writecell(ensemblID,'geneID.txt')

%%
% convert Ensembl gene ID to gene name:
% https://biotools.fr/mouse/ensembl_symbol_converter
ID2name = readtable('ensembleID2genename.txt');

B = A;
B = cat(2, ID2name(:,2), B(:,[2 5 3 4]));

% for now remove rows without gene name
missingidx = ismissing(ID2name(:,2));
B(missingidx,:) = [];

writetable(B,'named_count_matrix.txt');

% sequence mapped to multiple parts
% half a read or throw out?

%% check which FGF related genes are expressed

FGF1fam = 'FGF1$|FGF2$';
FGF4fam = 'FGF4$|FGF5$|FGF6$';
FGF7fam = 'FGF3$|FGF7$|FGF10$|FGF22$';
FGF9fam = 'FGF9$|FGF16$|FGF20$';
FGF8fam = 'FGF8$|FGF17$|FGF18$';
canonicalFGF = [FGF1fam '|' FGF4fam '|' FGF7fam '|' FGF9fam '|' FGF8fam];

endocrineFGF = 'FGF23$|FGF21$|FGF15$|FGF19'; % 15 and 19 are the same
intracellFGF = 'FGF11$|FGF12$|FGF13$|FGF14$';

% AS = antisense
% IT = intronic transcript
% OP2P1 = oncogene partner 2 pseudogene 1
% BP = binding protein

disp('canonical FGF')
genenames = canonicalFGF;
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

disp('FGF receptors')
genenames = 'FGFR.$|FGFRL';
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

disp('perlecan (part of HSPG)')
genenames = '^SDC.$';
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

disp('fate markers')
genenames = 'TBXT$|TBX6$|EOMES|SOX17|^SOX2$|CDX2';
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

disp('EMT')
genenames = 'SNAI1$|^CDH1$|CDH2$';
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

disp('other signals')
genenames = 'WNT|NODAL|BMP4';
idx = ~cellfun(@isempty, regexp(table2cell(B(:,1)), genenames));
result = sortrows(B(idx,:));
nonzeroidx = any(table2array(result(:,3:end)) > 0, 2);
result = result(nonzeroidx,:);
disp(result)

% large numbers -> not normalized to gene length?

%% check if columns are normalized to the same total read count

justnumbers = table2array(A(:,2:end));
sum(justnumbers,1)

% 1.0e+08 *
%     1.2654    1.2354    1.2992    1.2882

