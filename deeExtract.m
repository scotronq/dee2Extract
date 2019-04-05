function [dataOut,srrOut] = deeExtract(species, accession, keyword)
%deeExtract Extracts RNA-seq data from DEE2
%   DEE2 is a database with uniformly processed RNA-seq data for nearly all
%   GEO submissions (http://dee2.io/about.html). This function queries,
%   extracts, and formats data from this database.
%   
%   Inputs
%   species:	e.g. "athaliana","celegans","dmelanogaster","drerio",
%               "ecoli","hsapiens","mmusculus","rnorvegicus","scerevisiae"
%   accession:	e.g. SRR, SRX, SRS, SRP, GSE, GSM
%   keyword:    keyword search in experiment title (e.g. 'muscle')
%   
%   Outputs
%   dataOut:	MATLAB table with gene names, info, and RNA-seq
%   srrOut:     output when keyword used. Table of experiments that match
%               keyword
%
%   Resources
%   DEE2:               http://dee2.io/about.html
%   Code modified from: https://github.com/markziemann/dee2/blob/master/getDEE2.R
%   
%   Example
%   [temp1,temp2] = deeExtract('hsapiens',{'SRR1818593';'SRR3192352';'SRR3646715'});
%   [temp1,temp2] = deeExtract('hsapiens',[], 'muscle');
%   
%   Version 1.0 (03/27/19)
%   Written by: Scott Ronquist
%   Contact: 	scotronq@umich.edu
%   Created: 	03/27/19
%   
%   Revision History:
%   v1.0 (03/27/19)
%   * deeExtract.m created

%% set default parameters
if ~exist('species','var')||isempty(species);species='hsapiens';end
if ~exist('accession','var')||isempty(accession);accession={''};end
if ~exist('keyword','var')||isempty(keyword);keyword=[];end

% put accession into cell array, if not already
if ~iscell(accession)
    accession = {accession};
end

% get temp file name
tempFN = tempdir;

%% load species metadata
fprintf('loading species metadata...\n')
metadataFN = websave(fullfile(tempFN,sprintf('%s_metadata.tsv.cut',species)),...
    sprintf('http://dee2.io/metadata/%s_metadata.tsv.cut',species));
metadata = readtable(metadataFN,'FileType','text');

% if keyword not empty, check and output
if ~isempty(keyword)
    % find SRR locs that match keyword
    switch keyword(1:3)
        case 'GSM'
            srrOut = metadata(contains(metadata.GSM_accession,keyword),:);
        case 'GSE'
            srrOut = metadata(contains(metadata.GSM_accession,keyword),:);
        otherwise
            srrOut = metadata(contains(metadata.experiment_title,keyword),:);
    end
    
    % if no results found, output error
    if isempty(srrOut)
        error('No results found for keyword "%s"',keyword)
    end
    
    % user input select samples
    [indx,~] = listdlg('ListString',join([srrOut.SRR_accession srrOut.experiment_title]),...
        'ListSize',[1000 500]);
    accession = srrOut.SRR_accession(indx);
end

%% determine accession type and loc
% get accession specific data
switch accession{1}(1:3)
    case 'SRR'
        metadataLoc = ismember(metadata.SRR_accession,accession);
    case 'SRX'
        metadataLoc = ismember(metadata.SRX_accession,accession);
    case 'SRS'
        metadataLoc = ismember(metadata.SRS_accession,accession);
    case 'SRP'
        metadataLoc = ismember(metadata.SRP_accession,accession);
    case 'GSE'
        metadataLoc = ismember(metadata.GSE_accession,accession);
    case 'GSM'
        metadataLoc = ismember(metadata.GSM_accession,accession);
    otherwise
        error('deeExtract: unrecognized accession type "%s"',accession{1})
end
srrList = metadata.SRR_accession(metadataLoc);
srrOut = metadata(metadataLoc,:);

% if no results found, output error
if isempty(srrOut)
    error('No results found for input accession')
end

%% download data from SRR accession
fprintf('loading SRR data...\n')
options = weboptions('Timeout',length(srrList)*5);
websave(fullfile(tempFN,'tempDEE.zip'),sprintf('http://dee2.io/cgi-bin/request.sh?org=%s%s',...
    species,sprintf('&x=%s',srrList{:})),options);
unzipFN = unzip(fullfile(tempFN,'tempDEE.zip'),tempFN);

% load and format output
% http://dee2.io/data/hsapiens/hsa_gene_info.tsv, hsa_tx_info.tsv
dataOut = struct;
dataOut.geneInfo = readtable(fullfile(tempFN,'GeneInfo.tsv'),'FileType','text');
dataOut.geneCount = readtable(fullfile(tempFN,'GeneCountMatrix.tsv'),'FileType','text');
dataOut.geneCount = [dataOut.geneInfo(:,ismember(dataOut.geneInfo.Properties.VariableNames,'GeneSymbol')),...
    dataOut.geneCount];

dataOut.txInfo = readtable(fullfile(tempFN,'TxInfo.tsv'),'FileType','text');
dataOut.txCount = readtable(fullfile(tempFN,'TxCountMatrix.tsv'),'FileType','text');
dataOut.txCount = [dataOut.txInfo(:,ismember(dataOut.txInfo.Properties.VariableNames,'GeneSymbol')),...
    dataOut.txCount];

dataOut.qcMatrix = readtable(fullfile(tempFN,'QC_Matrix.tsv'),'FileType','text');
dataOut.metadataSummary = readtable(fullfile(tempFN,'MetadataSummary.tsv'),...
    'FileType','text','Delimiter','\t');
dataOut.metadataFull = readtable(fullfile(tempFN,'MetadataFull.tsv'),...
    'FileType','text','Delimiter','\t');

srrOut.SRR_accession
% get kallisto TPM data
for iAcc = 1:length(srrOut.SRR_accession)
    fprintf('loading TPM data (%i/%i)...\n',iAcc,length(srrOut.SRR_accession))
    
    % download and read TPM data
    tpmFnZip = websave(fullfile(tempFN,sprintf('%s.ke.tsv.gz',srrOut.SRR_accession{iAcc})),...
        sprintf('http://dee2.io/data/%s/%s/%s.ke.tsv.gz',species,...
        srrOut.SRR_accession{iAcc},srrOut.SRR_accession{iAcc}));
    tpmFn = gunzip(tpmFnZip);
    tpmTable = readtable(tpmFn{1},'FileType','text');
    
    % delete temp files
    delete(tpmFnZip)
    delete(tpmFn{1})
    
    % match with dataOut
    if ~isequal(dataOut.txInfo.TxID,tpmTable{:,1})
        error('deeExtract: TPM list not equal "%s"',srrOut.SRR_accession{iAcc})
    end
    
    % initialize formatted TPM tables
    if iAcc == 1
        dataOut.txTpm = dataOut.txCount(:,1:2);
        dataOut.geneTpm = dataOut.geneCount(:,1:2);
    end
    dataOut.txTpm = [dataOut.txTpm,...
        tpmTable(:,contains(tpmTable.Properties.VariableNames,'tpm'))];
    
    % calculate gene-level
    [b, ~, n] = unique(dataOut.txInfo.GeneSymbol , 'first');
    tpmColumn  = accumarray(n ,...
        tpmTable{:,contains(tpmTable.Properties.VariableNames,'tpm')},...
        size(b) , @(x) sum(x));
    
    % match with geneInfo
    % fprintf('WARNING: need to double check this!!\n')
    [LIA,LOCB] = ismember(dataOut.geneTpm.GeneSymbol,b);
    LOCB(LOCB==0) = [];
    temp = nan(height(dataOut.geneTpm),1);
    temp(LIA) = tpmColumn(LOCB);
    dataOut.geneTpm = [dataOut.geneTpm table(temp)];
    dataOut.geneTpm.Properties.VariableNames{end} = srrOut.SRR_accession{iAcc};
end

% delete temp files
delete(metadataFN)
delete(fullfile(tempFN,'tempDEE.zip'))
for iFn = 1:length(unzipFN)
    delete(unzipFN{iFn})
end

end
