# dee2Extract
MATLAB tool to extract DEE2 data

# Resources
[Digital Expression Explorer 2 (DEE2)](http://dee2.io/)

# MATLAB function details
DEE2 is a database with uniformly processed RNA-seq data for nearly all
GEO submissions (http://dee2.io/about.html). This function queries,
extracts, and formats data from this database.

Inputs
species:	e.g. "athaliana","celegans","dmelanogaster","drerio",
"ecoli","hsapiens","mmusculus","rnorvegicus","scerevisiae"
accession:	e.g. SRR, SRX, SRS, SRP, GSE, GSM
keyword:    keyword search in experiment title (e.g. 'muscle')

Outputs
dataOut:	MATLAB table with gene names, info, and RNA-seq
srrOut:     output when keyword used. Table of experiments that match
keyword

Resources
DEE2:               http://dee2.io/about.html
Code modified from: https://github.com/markziemann/dee2/blob/master/getDEE2.R

Example
[temp1,temp2] = deeExtract('hsapiens',{'SRR1818593';'SRR3192352';'SRR3646715'});
[temp1,temp2] = deeExtract('hsapiens',[], 'muscle');
