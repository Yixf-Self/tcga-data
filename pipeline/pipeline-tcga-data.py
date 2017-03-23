#################################################################
#################################################################
############### TCGA Data Pipeline ##############################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, os, urllib, json, glob, xmltodict, collections, itertools
import pandas as pd
import rpy2.robjects as robjects
import pandas.rpy.common as com
import numpy as np
import scipy.stats as ss

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineTcgaData as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
ensemblAnnotationFile = 'f2-rawdata.dir/annotations/ensembl_hgnc_annotation.txt'
proteomicsFile = 'f2-rawdata.dir/proteomics/TCGA-PANCAN19-L4.csv'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-tcga-data.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Download Data
#######################################################
#######################################################

#############################################
########## 1. Create payload
#############################################

@follows(mkdir('f1-files.dir'))

@originate('f1-files.dir/tcga-payload.txt')

def getTcgaPayload(outfile):

	# Create payload
	payload = '''

	{
		"op":"and",
		"content":[
			{
				"op":"=",
				"content":{
					"field":"access",
					"value":"open"
				}
			},
			{
				"op":"in",
				"content":{
					"field":"files.data_type",
					"value":["Gene Expression Quantification", "Biospecimen Supplement", "Clinical Supplement"]
				}
			}
		]
	}

	'''

	# Encode
	payloadEncoded = urllib.quote(payload.replace('\n', '').replace('\t', ''))

	# Write file
	with open(outfile, 'w') as openfile:
		openfile.write(payloadEncoded)

#############################################
########## 2. Download JSON
#############################################

@transform(getTcgaPayload,
		   suffix('payload.txt'),
		   'files.json')

def getTcgaFiles(infile, outfile):

	# Read filter
	with open(infile, 'r') as openfile:
		filterString = openfile.read().strip()

	# Create command string
	commandString = ''' curl 'https://gdc-api.nci.nih.gov/files?
							format=JSON&
							pretty=true&
							expand=analysis,cases.project,associated_entities&
							size=100000&
							filters=%(filterString)s' > %(outfile)s ''' % locals()

	# Run command
	os.system(commandString.replace('\n','').replace('\t', ''))

#############################################
########## 3. Make Table
#############################################

@transform(getTcgaFiles,
		   suffix('.json'),
		   '.txt')

def makeTcgaFileTable(infile, outfile):

	# Read infile
	with open(infile, 'r') as openfile:
		tcgaData = json.loads(openfile.read())

	# Create filtered dict
	resultDict = {}

	# Loop through results
	for result in tcgaData['data']['hits']:
		
		# Add result
		resultDict[result['file_id']] = {'project_id': result['cases'][0]['project']['project_id'],
										 'data_type': result['data_type'],
										 'workflow_type': result['analysis']['workflow_type'] if 'analysis' in result.keys() else None,
										 'file_name': result['file_name'],
										 'entity_submitter_id': result['associated_entities'][0]['entity_submitter_id']}

	# Convert to dataframe
	resultDataframe = pd.DataFrame(resultDict).T

	# Save
	resultDataframe.to_csv(outfile, sep='\t', index_label='file_id')

#############################################
########## 4. Make Payloads
#############################################

@follows(mkdir('f1-files.dir/payloads'))

@subdivide(makeTcgaFileTable,
		   formatter(),
		   'f1-files.dir/payloads/payload-*.txt',
		   'f1-files.dir/payloads/payload')

def splitPayloads(infile, outfiles, outfileRoot):

	# Read file
	fileDataframe = pd.read_table(infile)

	# Filter
	fileDataframe = fileDataframe[fileDataframe['workflow_type'] != 'HTSeq - FPKM']
	fileDataframe = fileDataframe[fileDataframe['data_type'] != 'Biospecimen Supplement']

	# Set index
	fileDataframe.index = [x.lower().replace(' ', '_') for x in fileDataframe['data_type']]

	# Loop through index
	for dataType in set(fileDataframe.index):

		# Get file Ids
		fileIds = fileDataframe.loc[dataType, 'file_id'].tolist()

		# Set counter
		i = 0
		
		# Split
		for idList in np.array_split(fileIds, max([len(fileIds)/100, 1])):
			
			# Increase counter
			i += 1

			# Get ID string
			idString = '",\n"'.join(idList)

			# Create payload
			payloadString = '{"ids":["%(idString)s"]}\n' % locals()

			# Write payload
			outfile = '{outfileRoot}-{dataType}-{i}.txt'.format(**locals())
			with open(outfile, 'w') as openfile:
				openfile.write(payloadString)

#############################################
########## 5. Download
#############################################

@follows(mkdir('f2-rawdata.dir/download'))

@transform(splitPayloads,
		   regex(r'.*/payload-(.*).txt'),
		   r'f2-rawdata.dir/download/download-\1.tar.gz')

def downloadTcgaData(infile, outfile):

	commandString = 'curl -o %(outfile)s --remote-name --remote-header-name --request POST --header "Content-Type: application/json" --data @%(infile)s "https://gdc-api.nci.nih.gov/data"' % locals()

	print 'Doing '+outfile+'...'

	os.system(commandString)

#############################################
########## 6. Decompress
#############################################

@follows(mkdir('f2-rawdata.dir/data'))

@transform(downloadTcgaData,
		   regex(r'.*/(.*).tar.gz'),
		   r'f2-rawdata.dir/data/\1')

def decompressTcgaData(infile, outfile):

	# Create dir
	os.makedirs(outfile)

	# Run command
	os.system('tar -zxvf %(infile)s -C %(outfile)s' % locals())

#######################################################
#######################################################
########## S2. Prepare Data
#######################################################
#######################################################

#############################################
########## 1. Make directories
#############################################

@subdivide(makeTcgaFileTable,
		   formatter(),
		   'T*-*')

def makeDirectories(infile, outfiles):

	# Read data
	fileDataframe = pd.read_table(infile)

	# Loop through projects
	for projectId in set(fileDataframe['project_id']):
		outDir = '{projectId}'.format(**locals())
		if not os.path.exists(outDir):
			os.makedirs(outDir)

#############################################
########## 1. Transcriptomics
#############################################

def transcriptionJobs():

	# Get infile
	infile = 'f1-files.dir/tcga-files.txt'

	# Read files
	fileDataframe = pd.read_table(infile)

	# Add workflow
	fileDataframe['workflow'] = [x.split(' - ')[1].lower() if type(x) == str and ' - ' in x else x for x in fileDataframe['workflow_type']]

	# Filter
	fileDataframe = fileDataframe.replace('fpkm', np.nan).dropna()

	# Loop
	for projectId in set(fileDataframe['project_id']):
		for workflow in set(fileDataframe['workflow']):
			outfile = '%(projectId)s/%(projectId)s-%(workflow)s.txt' % locals()
			yield [[infile, ensemblAnnotationFile], outfile]

@follows(makeDirectories)
@files(transcriptionJobs)

def prepareTranscriptionData(infiles, outfile):

	# Report
	print 'Doing '+outfile+'...'

	# Split infiles
	tableFile, ensemblAnnotationFile = infiles

	# Read files
	fileDataframe = pd.read_table(tableFile)
	ensemblAnnotationDataframe = pd.read_table(ensemblAnnotationFile).dropna().rename(columns={'Gene ID': 'ensembl_gene_id', 'HGNC symbol': 'gene_symbol'})

	# Fix workflow
	fileDataframe['workflow'] = [x.split(' - ')[1].lower() if type(x) == str and ' - ' in x else x for x in fileDataframe['workflow_type']]

	# Get data
	projectId = os.path.dirname(outfile)
	workflow = '-'.join(outfile.replace('.txt', '').split('-')[3:])

	# Filter
	fileDataframe = fileDataframe.replace('fpkm', np.nan).dropna().set_index(['project_id', 'workflow']).loc[(projectId, workflow)]

	# Get file paths
	fileDataframe['file_path'] = [glob.glob('f2-rawdata.dir/data/*/%(x)s/*' % locals())[0] for x in fileDataframe['file_id']]

	# Get sample id
	fileDataframe['sample_id'] = ['-'.join(x.split('-')[:4])[:-1] for x in fileDataframe['entity_submitter_id']]

	# Get file dict
	fileDict = fileDataframe[['file_path', 'sample_id']].set_index('file_path').to_dict()['sample_id']

	# Initialize list
	dataList = []

	# Loop through data
	for filePath, entitySubmitterId in fileDict.iteritems():
		expressionDataframe = pd.read_table(filePath, names=['ensembl_gene_id', 'expression'])
		expressionDataframe['sample_id'] = entitySubmitterId
		dataList.append(expressionDataframe)

	# Prepare matrix
	concatenatedDataframe = pd.concat(dataList)

	# Fix rownames
	concatenatedDataframe['ensembl_gene_id'] = [x.split('.')[0] for x in concatenatedDataframe['ensembl_gene_id']]

	# Annotate
	annotatedDataframe = concatenatedDataframe.merge(ensemblAnnotationDataframe, on='ensembl_gene_id', how='inner')

	# Aggregate
	aggregatedDataframe = annotatedDataframe.groupby(['gene_symbol', 'sample_id']).mean().astype(annotatedDataframe['expression'].dtype).reset_index()

	# Pivot
	resultDataframe = aggregatedDataframe.pivot(index='gene_symbol', columns='sample_id', values='expression')

	# Fix gene IDs
	resultDataframe.index = [x.split('.')[0] for x in resultDataframe.index]

	# Write
	resultDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#############################################
########## 2. Survival
#############################################

@subdivide(makeTcgaFileTable,
		   formatter(),
		   '*/*-survival.txt',
		   '-survival.txt')

def makeSurvivalMatrices(infile, outfiles, outfileRoot):

	# Read files
	fileDataframe = pd.read_table(infile, index_col='project_id')

	# Filter
	fileDataframe = fileDataframe[fileDataframe['data_type'] == 'Clinical Supplement']

	# Get files
	fileDataframe['file_path'] = [glob.glob('f2-rawdata.dir/data/*/%(x)s/*' % locals())[0] for x in fileDataframe['file_id']]

	# Get tumor types
	tumorTypes = set(fileDataframe.index)

	# Loop through projects
	for projectId in tumorTypes:

		# Report
		print 'Doing '+projectId+'...'
		
		# Get file paths
		filePaths = [x for x in fileDataframe.loc[projectId, 'file_path'] if x.split('.')[-1] == 'xml']

		# Filter
		if len(filePaths) > 1:

			# Get tumor type
			tumorType = projectId.split('-')[-1].lower()
			
			# Initialize list
			dataList = []
			
			# Loop through files
			for filePath in filePaths:

				# Get patient dict
				with open(filePath, 'r') as openfile:
					patientDict = xmltodict.parse(openfile.read())[tumorType+':tcga_bcr'][tumorType+':patient']

				# Get data
				patientBarcode = patientDict['shared:bcr_patient_barcode']['#text']
				vitalStatus = patientDict['clin_shared:vital_status']['#text'] if '#text' in patientDict['clin_shared:vital_status'].keys() else None
				daysToDeath = patientDict['clin_shared:days_to_death']['#text'] if '#text' in patientDict['clin_shared:days_to_death'].keys() else None
				daysToLastFollowup = patientDict['clin_shared:days_to_last_followup']['#text'] if '#text' in patientDict['clin_shared:days_to_last_followup'].keys() else None	
				
				# Add data
				dataList.append([patientBarcode, vitalStatus, daysToDeath, daysToLastFollowup])

			# Convert to dataframe
			survivalDataframe = pd.DataFrame(dataList, columns=['patient_barcode', 'vital_status', 'days_to_death', 'days_to_last_followup'])

			# Add last checked
			# survivalDataframe['last_checked'] = [days_to_death if not np.isnan(float(days_to_death)) else days_to_last_followup for days_to_death, days_to_last_followup in survivalDataframe[['days_to_death', 'days_to_last_followup']].as_matrix()]

			# Save
			outfile = '{projectId}/{projectId}{outfileRoot}'.format(**locals())
			survivalDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 3. Proteomics
#############################################

@subdivide(proteomicsFile,
		   formatter(),
		   '*/*-rppa.txt',
		   '-rppa.txt')

def makeRppaMatrices(infile, outfiles, outfileRoot):

	# Read dataframe
	rppaDataframe = pd.read_csv(infile, index_col='Cancer_Type').drop(['Sample_Type', 'SetID'], axis=1)

	# Fix samples
	rppaDataframe['Sample_ID'] = ['-'.join(x.split('-')[:4])[:-1] for x in rppaDataframe['Sample_ID']]

	# Get unique tumor types
	tumorTypes = set(rppaDataframe.index)

	# Loop
	for tumorType in tumorTypes:

		# Get subset
		rppaDataframeSubset = rppaDataframe.loc[tumorType].reset_index(drop=True)

		# Melt and aggregate
		rppaDataframeMelt = pd.melt(rppaDataframeSubset, id_vars='Sample_ID').dropna().groupby(['Sample_ID', 'variable']).mean().reset_index()

		# Cast
		rppaDataframeCast = rppaDataframeMelt.pivot(index='variable', columns='Sample_ID', values='value').fillna(0)

		# Save
		tumorType = 'TCGA-'+tumorType
		outfile = '{tumorType}/{tumorType}{outfileRoot}'.format(**locals())
		rppaDataframeCast.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S2. Analyze Data
#######################################################
#######################################################

#############################################
########## 1. DESeq2
#############################################

@subdivide(glob.glob('*/*-counts.txt'),
		   regex(r'.*/(.*)-counts.txt'),
		   r'\1/deseq/\1-deseq-*.txt',
		   r'\1/deseq/\1-deseq-')

def runDESeq(infile, outfiles, outfileRoot):

	# Report
	print 'Doing '+infile+'...'

	# Read dataframe
	countDataframe = pd.read_table(infile, index_col='gene_symbol')

	# Sample counts
	sampleCounts = collections.Counter([x.split('-')[-1] for x in countDataframe.columns])

	# Make annotation dataframe
	annotationDataframe = pd.DataFrame.from_dict([{'sample_id': x, 'sample_type': x.split('-')[-1]} for x in countDataframe.columns]).set_index('sample_id')

	# Sample counts
	sampleCounts = collections.Counter([x.split('-')[-1] for x in countDataframe.columns])

	# Get comparisons
	comparisons = [list(x[::-1]) for x in itertools.combinations([key for key, value in sampleCounts.iteritems() if value >= 5], 2)]

	# Loop through comparisons
	for comparison in comparisons:
		
		# Filter
		annotationDataframeSubset = annotationDataframe[annotationDataframe['sample_type'].isin(comparison)]
		countDataframeSubset = countDataframe[annotationDataframeSubset.index]

		# Run function
		deseqDataframe = r.runDESeq2(com.convert_to_r_dataframe(countDataframeSubset), com.convert_to_r_dataframe(annotationDataframeSubset), '~ sample_type')

		# Convert to dataframe
		deseqDataframe = com.convert_robj(deseqDataframe)

		# Get comparison string
		comparisonString = 'v'.join(comparison)

		# Get outfile
		outfile = '{outfileRoot}{comparisonString}.txt'.format(**locals())

		# Create outdir
		outDir = os.path.dirname(outfile)
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		# Write
		deseqDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#############################################
########## 2. Differential Expression
#############################################

@subdivide(glob.glob('*/*-fpkm-uq.txt'),
		   regex(r'.*/(.*)-fpkm-uq.txt'),
		   r'\1/de/\1-de-*.txt',
		   r'\1/de/\1-de-')

def runDifferentialExpression(infile, outfiles, outfileRoot):

	# Report
	print 'Doing '+infile+'...'

	# Read dataframe
	expressionDataframe = pd.read_table(infile, index_col='gene_symbol')

	# Z-score
	expressionDataframe = ((expressionDataframe.T - expressionDataframe.T.mean())/expressionDataframe.T.std()).T

	# Create dicts
	sampleDict = {x: x.split('-')[-1] for x in expressionDataframe.columns}
	sampleTypeDict = {sampleType:[key for key, value in sampleDict.iteritems() if value == sampleType] for sampleType in set(sampleDict.values())}
	sampleCountDict = {x:len(sampleTypeDict[x]) for x in sampleTypeDict.keys()}

	# Get comparisons
	comparisons = [list(x[::-1]) for x in itertools.combinations([key for key, value in sampleCountDict.iteritems() if value >= 5], 2)]

	# Loop through comparisons
	for comparison in comparisons:

		# Initialize result dict
		resultDict = {x:{} for x in expressionDataframe.index}

		# Get samples
		treatedSamples = sampleTypeDict[comparison[0]]
		controlSamples = sampleTypeDict[comparison[1]]

		# Loop through genes
		for index, rowData in expressionDataframe.iterrows():
			# run ttest
			tResults = ss.ttest_ind(rowData[treatedSamples], rowData[controlSamples])
			# add data
			resultDict[index]['statistic'] = tResults.statistic
			resultDict[index]['pvalue'] = tResults.pvalue
			
		# Convert to dataframe
		resultDataframe = pd.DataFrame(resultDict).T

		# Get comparison string
		comparisonString = 'v'.join(comparison)

		# Get outfile
		outfile = '{outfileRoot}{comparisonString}.txt'.format(**locals())

		# Create outdir
		outDir = os.path.dirname(outfile)
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		# Write
		resultDataframe.sort_values('pvalue').to_csv(outfile, sep='\t', index_label='gene_symbol')
		print 'Done '+outfile+'!'

#############################################
########## 3. Survival Association
#############################################

@transform(glob.glob('*/*-survival.txt'),
		   regex(r'.*/(.*)-survival.txt'),
		   add_inputs(r'\1/\1-fpkm-uq.txt'),
		   r'\1/\1-survival_association.txt')

def getSurvivalAssociation(infiles, outfile):

	# Split infiles
	survivalFile, fpkmFile = infiles

	# Report
	print 'Doing '+survivalFile+'...'

	# Read dataframes
	expressionDataframe = pd.read_table(fpkmFile, index_col='gene_symbol')
	survivalDataframe = pd.read_table(survivalFile, index_col='patient_barcode')

	# Fix survival dataframe
	survivalDataframe['last_checked'] = [days_to_death if not np.isnan(days_to_death) else days_to_last_followup for days_to_death, days_to_last_followup in survivalDataframe[['days_to_death', 'days_to_last_followup']].as_matrix()]
	survivalDataframe['event'] = [x=='Dead' for x in survivalDataframe['vital_status']]
	
	# Get column types
	cols = [x for x in expressionDataframe.columns if x.split('-')[-1]=='01']
	expressionDataframeSubset = expressionDataframe[cols]

	# Fix patient columns
	expressionDataframeSubset.columns = ['-'.join(x.split('-')[:3]) for x in expressionDataframeSubset.columns]

	# Aggregate
	expressionDataframeSubset = expressionDataframeSubset.T.reset_index().groupby('index').mean().T

	# Get common patients
	commonPatients = list(set(expressionDataframeSubset.columns).intersection(set(survivalDataframe.index)))
	expressionDataframeSubset = expressionDataframeSubset[commonPatients]

	# Result
	resultDict = {}

	# Loop through gene symbols
	for geneSymbol in expressionDataframeSubset.index:

		# Split samples
		split_cols = np.array_split(expressionDataframeSubset.loc[geneSymbol].to_frame().sort_values(geneSymbol).index.tolist(), 3)
		groupDataframe = pd.DataFrame([[x, i] for i in range(len(split_cols)) for x in split_cols[i]], columns=['patient_barcode','group']).set_index('patient_barcode')

		# Check order
		survivalDataframeSubset = survivalDataframe.loc[groupDataframe.index, ['event', 'last_checked']]

		# Get pvalue
		try:
			survivalP = r.getSurvivalAssociation(com.convert_to_r_dataframe(survivalDataframeSubset), com.convert_to_r_dataframe(groupDataframe))
		except:
			survivalP = [1]

		# Add to dict
		resultDict[geneSymbol] = {'pvalue': survivalP[0]}

	# Convert to dataframe
	resultDataframe = pd.DataFrame(resultDict).T

	# Save
	resultDataframe.sort_values('pvalue').to_csv(outfile, sep='\t', index_label='gene_symbol')


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
