
import argparse, platform
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
import concurrent.futures, numpy, sys
import os.path
import re

# Validates user-provided command-line arguments
class ParameterValidator():
	def __init__(self, args):
		self.args = args
		self.check_python()
		self.check_biopython()
		self.check_args()

	# Check python 3.2 is available
	def check_python(self):
		if platform.python_version_tuple()[0:2] >= ('3', '0'): # >= 3.0
			out('Python v. '+ str(platform.python_version()) + ' found [OK]')
		else:
			raise RuntimeError('Python 3.2 required')

	# BioPython must be installed
	def check_biopython(self):
		try:
			import Bio
			out('BioPython v.' + Bio.__version__ + ' found [OK]')
		except ImportError:
			raise ImportError('Make sure BioPython is installed')

	# Checks user-provided arguments are valid
	def check_args(self):
		return all([self.test_num_workers(), self.test_mutual_matrices(),
				self.test_valid_matrix()])

	# Test either a custom matrix or in-built matrix is selected
	def test_mutual_matrices(self):
		if self.args['custom'] and self.args['matrix']:
			raise IOError('Either a custom or in-built matrix must be selected')
		else:
			return True

	# Test a valid substitution matrix is selected
	def test_valid_matrix(self):
		all_matrices = set(MatrixInfo.available_matrices) # al sub. matrices
		if self.args['matrix'] in all_matrices or self.args['custom']:
			return True
		else:
			err = 'An in-built matrix (below) or custom matrix must be provided\n'+\
			'http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html'
			raise IOError(err)

	# Test a valid number of workers are provided
	def test_num_workers(self):
		if self.args['n'] >= 1:
			return True
		else:
			raise IOError('>= 1 worker processes must be provided')

# Helper-function to write a string
def out(s):
	sys.stdout.write(s+'\n')

# Trivial function to write parameter arguments to a file 
def write_args(args):
	outhandle = open('param_args.tab', 'w')
	for i in args:
		outhandle.write(i+'\t'+str(args[i])+'\n') # write args
		outhandle.flush()
	outhandle.close()
	out('') # write new line

# Trivial function to parse a fasta file
def parse_fasta(fname):
	queries = list(SeqIO.parse(fname, 'fasta')) # easy indexing
	for query in queries:
		query.name = clean_header(query.name)
	out(str(len(queries)) + ' queries parsed [OK]')
	return queries # return set of fasta entries

def clean_header(header):
	m = re.search('^([^ \\|]+) \\|.*',header)
	if m:
		return m.group(1)
	return header
				
# Maps the aligned two bases against a user-selected substitution matrix
def get_score(cA, cB, submatrix):
	if (cA, cB) in submatrix:
		return submatrix[(cA, cB)]
	else:
		return submatrix[(cB, cA)] # returns score

# Performs Needleman-Wunsch alignment given two sequences, s1 and s2	
def needle(seq1, seq2, gap, submat):
	l1, l2 = len(seq1.seq), len(seq2.seq)	
	m = numpy.zeros((l1+1, l2+1)) # create matrix for storing counts
	for i in range(0, l1 + 1): # set each row by the desired gap
		m[i][0] = gap * i
	for j in range(0, l2 + 1): # set each column by the desired gap
		m[0][j] = gap * j
	for i in range(1, l1 + 1): # per base-pair in sequence 1 ...
		for j in range(1, l2 + 1): # per base-pair in sequence 2, align them
			score = m[i - 1][j - 1] + get_score(seq1.seq[i-1], seq2.seq[j-1], submat)
			left = m[i - 1][j] + gap # upwards
			up = m[i][j - 1] + gap # left score
			m[i][j] = max(score, left, up) # get max of all three scores

	a1, a2 = '', '' # when complete, identify both aligned sequences
	i,j = l1,l2
	while i > 0 and j > 0: # walk-back to the index [0][0] of the m
		score_current = m[i][j]
		score_diag = m[i-1][j-1]
		score_up = m[i][j-1]
		score_left = m[i-1][j]

		# if the score is a match, walk-back one index in both i and j
		if score_current == score_diag + get_score(seq1.seq[i-1], seq2.seq[j-1], submat):
			a1 += seq1.seq[i-1]
			a2 += seq2.seq[j-1]
			i -= 1
			j -= 1
		# if score is a gap in sequence 2, only walk back on i
		elif score_current == score_left + gap:
			a1 += seq1.seq[i-1]
			a2 += '-'
			i -= 1
		# if score is a gap in sequence 1, only walk back on j
		elif score_current == score_up + gap:
			a1 += '-'
			a2 += seq2.seq[j-1]
			j -= 1
	# walk-back to index 0 for both i and j; either could be reached first
	while i > 0:
		a1 += seq1.seq[i-1]
		a2 += '-'
		i -= 1
	while j > 0:
		a1 += '-'
		a2 += seq2.seq[j-1]
		j -= 1
	return [ m[-1][-1], a1+'\l2'+a2, seq2.description ]

# Function to parse custom scoring matrix.
def parse_custom_matrix(fname):
	""" 
	A schema is organized such that you have 3 columns: A, B, C.
	Columns A and B represents the query and target base, respectively.
	Column C represents the real-number assigned to the mismatch given A and B.
	This file must be tab-delimited.
	Example:
	A	T	-5
	A	C	2
	A	G	-2
	A	A	8
	T	A	2 
	... 
	etc. 
	"""
	submat = {} # the substitution matrix
	for line in open(fname):
		line = line.strip()
		if len(line) == 0: # if an empty line, terminate analysis
			break
		else: # pertains to parsing of the custom substition matrix
			line = line.split('\t')
			if len(line) < 3:
				raise IndexError('Custom matrix must have 3x columns')
			else: # set the Query (A) and Target (B) scores
				a, b, score = line[0], line[1], float(line[2])
				submat[(a, b)] = score
	return submat

# Function to parse a preexisting version of the file that will be written to, returning a list of sequences already aligned
def parse_output(fname):
	alreadyDone = []
	for line in open(fname):
		line = line.strip()
		if len(line) == 0:
			break
		else:
			sequenceName = line.split('\t',1)[0]
			alreadyDone.append(sequenceName)
	return alreadyDone
	
# Concurrent alignment given a sequence, query, and a set of sequences, baseline
def run_factory(target, queries, params):
	results = [] # K => target, V => aligned queries 
	for query in queries:
		out = needle(target, query, params['gap'], params['submat'])
		results.append(out)
	return target, results

# Performs the high-level functions which drive concurrent execution
def initializer(queries, args):
	if not args['matrix']: # if no in-built matrix, parse custom matrix
		submat = parse_custom_matrix(fname=args['custom']) # custom matrix
	else:
		submat = getattr(MatrixInfo, args['matrix']) # get substitution matrix
		
	# get all completed jobs
	outputFile = args['o']
	if outputFile == 'scores.tab' or not os.path.exists(outputFile):
		alreadyRun = []
	else:
		alreadyRun = parse_output(outputFile)

	executor = concurrent.futures.ProcessPoolExecutor(max_workers=args['n'])
	futures = [] # create collection to store all concurrent jobs in
	params = {'gap': args['gap'], 'submat': submat} # alignment parameters
	tmpCount = 0
	for target in queries: # per fasta entry, create a concurrent job for it
		if not target.name in alreadyRun:
			futures.append(executor.submit(run_factory, target, queries, params))
			
	# get all completed jobs
	if outputFile == 'scores.tab' or not os.path.exists(outputFile):
		outhandle = open(outputFile, 'w')
	else:
		outhandle = open(outputFile, 'a')
		
	for counter, future in enumerate(concurrent.futures.as_completed(futures)):
		# print percentage complete
		perc = round((float(counter+1+len(alreadyRun)) / len(queries)) * 100, 4)
		sys.stdout.write('\r[%d%% complete] ' % (perc))

		target, results = future.result()
		results = sorted(results, key=lambda x: x[-1]) # sort by query (last item)
		if len(alreadyRun) == 0 and counter == 0:
			headers = [h[-1] for h in results]
			outhandle.write('Input\t' + '\t'.join(headers) + '\n')
			outhandle.flush()
		# get scores for each sequence and output to file
		scores = '\t'.join([str(s[0]) for s in results])
		outhandle.write(target.description + '\t' + scores + '\n')
		outhandle.flush()
	outhandle.close()

if __name__ == '__main__':
	desc = 'Script to execute exhaustive brute-force pairwise alignment'
	u='%(prog)s [options]' # command-line usage
	p = argparse.ArgumentParser(description=desc, add_help=False, usage=u)
	param_reqd = p.add_argument_group('Required Parameters')
	param_opts = p.add_argument_group('Optional Parameters')

	# Specify required arguments
	param_reqd.add_argument('-f', metavar='FILE', required=True,
				help='Input fasta file [na]')

	# Specify optional arguments
	param_opts.add_argument('--gap', metavar='INT', default=-8, type=int,
				help='Gap open and extension penalty [-8]')

	param_opts.add_argument('-custom', metavar='FILE', default=None,
				help='Custom substitution matrix [na]')
	param_opts.add_argument('-matrix', metavar='STR', default=None,
				help='Matrix name; see Biopython MatrixInfo for all matrices [na]')
	param_opts.add_argument('-n', metavar='INT', default=2, type=int,
				help='Number of worker processes [2]')
	param_opts.add_argument('-o', metavar='FILE', default='scores.tab', 
				help='File to which output should be writen/appended')
	param_opts.add_argument('-h','--help', action='help',
				help='Show this help screen and exit')
	args = vars(p.parse_args()) # parse user-provided Arguments
	try:
		ParameterValidator(args)
		queries = parse_fasta(fname=args['f']) # parse fasta file
		initializer(queries, args)
		write_args(args) # write arguments to a file
	except (IOError, KeyboardInterrupt, IndexError) as e:
		out(e+'\n')
		
