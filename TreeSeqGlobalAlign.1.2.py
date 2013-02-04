
import numpy, sys, math

# Creates a dictionary for a given tree sequence linking each T node to an associated A node
def create_ta_dictionary(seq,nodeTypes,submatrix,gap_cost):
	# Dictionary with (key,value) pairs of (T-node index, A-node index) as well as (str(T-node index), total gap cost)
	taDict = {}
	AStack = []
	AStack.append(-1) # Begin the AStack with a sentinal value
	CostStack = []
	CostStack.append(0)
	# Visit each character in the sequence
	for index in range(0,len(seq)):
		# Add the cost of the current node to the current gap cost register
		CostStack[-1] += submatrix[seq[index],'-']

		# If the current node is an A-node
		if seq[index] in nodeTypes['A']:
			# push it onto the top of the stack
			AStack.append(index)
			# push on a new gap cost register
			CostStack.append(0)
			# The cost of an A-node at the front of a gap is added on in
			# NeedlemanWunsch.calculate_gap when determining whether the A is gapped or is
			# matched to a C-node

		# If the current node is a T-node
		elif seq[index] in nodeTypes['T']:
			# Set the dictionary to contain the (T-node index, associated A-node pair)
			taDict[index] = AStack.pop()
			# Get the cost of the gap between the T and A nodes
			currCost = CostStack.pop()
			# Add that cost to the dictionary
			taDict[str(index)] = currCost
			# Add the cost to the enclosing AT pair
			if (len(CostStack) > 0):
				CostStack[-1] += currCost
				
	return taDict

# Implementation of global alignment - Needleman-Wunsch
class NeedlemanWunsch():
	def __init__(self, s1, s2, costs, submat, nodeTypes):
		self.seq1 = s1 # sequence 1
		self.seq2 = s2 # sequence 2
		self.costs = costs # dictionary of all costs (i.e. penalties)
		self.submat = submat # substitution matrix
		self.create_node_types(nodeTypes)
		self.create_residue_specific_gapcost()
		self.TADict1 = create_ta_dictionary(s1.seq, nodeTypes, submat, costs['gap'])
		self.TADict2 = create_ta_dictionary(s2.seq, nodeTypes, submat, costs['gap'])
		self.scoreMat = None # references score matrix
		self.directionMat = None # references diag(0),left(1),up(2) matrix
		self.leftMat = None # references diag(0),left(1),up(2) matrix
		self.upMat = None # references diag(0),left(1),up(2) matrix
		self.backPos = {} # the backtrace position from one position to its prior (contains integer pairs)
		self.align1 = '' # alignment string for sequence 1
		self.align2 = '' # alignment string for sequence 2
		self._aligner()
	
	def create_node_types(self,nodeTypes):
		self.nodeTypes = {}
		for nodeType in nodeTypes.keys():
			for residue in nodeTypes[nodeType]:
				self.nodeTypes[residue] = nodeType
	
	# Fills in all gap-residue pairs either with flipped order entry if it exists, else with the default gap cost
	# Also fills in flipped residue-residue scores
	def create_residue_specific_gapcost(self):
		newToSubmat = {}
		for pair in self.submat.keys():
			if pair[0] == '-' and pair[1] == '-':
				# do nothing, this is useless and shouldn't happen
				pass
			elif pair[0] == '-' or pair[1] == '-':
				# Note that the given residue has an associated gap cost
				if pair[0] == '-' and (pair[1],'-') not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.submat[pair]
				elif pair[1] == '-' and (pair[1],'-') not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.submat[pair]
			else:
				# Fill the residueDict so none are missed
				if (pair[0],'-') not in self.submat.keys() and ('-',pair[0]) not in self.submat.keys():
					newToSubmat[pair[0],'-'] = self.costs['gap']
					newToSubmat['-',pair[0]] = self.costs['gap']
				if (pair[1],'-') not in self.submat.keys() and ('-',pair[1]) not in self.submat.keys():
					newToSubmat[pair[1],'-'] = self.costs['gap']
					newToSubmat['-',pair[1]] = self.costs['gap']
				if (pair[1],pair[0]) not in self.submat.keys():
					newToSubmat[pair[1],pair[0]] = self.submat[pair]
		for pair in newToSubmat.keys():
			self.submat[pair] = newToSubmat[pair]
	
	# return top (highest) alignment score given sequence 1 and 2
	def get_top_score(self):
		return self.scoreMat[-1][-1]
	
	# alignment string resultant from sequence 1 and 2
	def get_alignment(self):
		return (self.align1, self.align2)
	
	# Create parser-friendly output given a NW alignment 
	def prettify(self):
		return [ self.get_top_score(), self.get_alignment(), self.seq2.name ]

	def determine_open_extend(self,i,j,m,directionM,dirScoreM,currentGapCost,gapDirection):
		gapScore = m[i][j] + currentGapCost
		# Add on the gap open cost if the prior position is not gapped in the same direction (gapping the same sequence)
		if math.isnan(dirScoreM[i,j][0]): # Previous position can't gap
			gapScore += self.costs['gapopen']
			scoreExtendPair = gapScore,False
		elif directionM[i,j] is not gapDirection: # Previous position didn't choose gap
			extendGapScore = dirScoreM[i][j][0] + currentGapCost
			openGapScore = gapScore + self.costs['gapopen']
			if openGapScore >= extendGapScore: # new gap is best, go with previous position's choice
				scoreExtendPair = openGapScore,False
			else: 
				# continuation is best, override previous position choice; if current position is on path, previous position will gap
				scoreExtendPair = extendGapScore,True
		else: # Previous position did choose gap, so this is automatically a continuation
			scoreExtendPair = gapScore,True
		return scoreExtendPair
		
	# Calculates the gap cost in a given direction from a given position, which depends on the node type
	def calculate_gap(self,i,j,seq1,seq2,m,directionM,dirScoreM,TADict,gapDirection):
		isExtend = False
		# CType: gap one
		if self.nodeTypes[seq1[i-1]] == 'C':
			# The prior position assuming a gap (index based on m)
			gapPosi = i - 1
			gapPosj = j
			# Determine the appropriate gap score depending on whether this opens or extends a gap
			gapScore,isExtend = self.determine_open_extend(gapPosi,gapPosj,m,directionM,dirScoreM,get_gapcost(seq1[i-1],self.submat),gapDirection)

		# TType: gap until paired A
		elif self.nodeTypes[seq1[i-1]] == 'T': 
			# The prior position assuming a gap (index based on m)
			gapPosi = TADict[i-1]
			# Case where this is the last T; handle sentinal and get cost of front-gap
			if TADict[i-1] is -1:
				gapPosi = 0
			gapCostMajor = TADict[str(i-1)] # Cost of the gap from the T-node up to the A-node
			gapCostStart = get_gapcost(seq1[gapPosi],self.submat) # Cost of the A-node that starts the gap

			# Calculate the total gap cost assuming the associated A is also gapped
			gapScoreGapFinish,isExtend = self.determine_open_extend(gapPosi,j,m,directionM,dirScoreM,gapCostMajor+gapCostStart,gapDirection)

			# If seq2 character is C-type, determine whether to match T-paired A and C, or to just gap the A
			# This will not happen if this is the last T (TADict[i-1] is -1), as the whole sequence must be gapped
			if self.nodeTypes[seq2[j-1]] == 'C' and TADict[i-1] is not -1:
				# Calculate the total gap cost assuming the associated A-node matches a C-node
				ACScore = get_score(seq1[gapPosi],seq2[j-1],self.submat)
				gapScoreACFinish = m[gapPosi][j-1] + gapCostMajor + ACScore + self.costs['gapopen']
				# Determine which gap produces a higher overall score, and use that for this position's gap score
				if gapScoreACFinish >= gapScoreGapFinish:
					gapScore = gapScoreACFinish
					gapPosj = j-1
					isExtend = False
				else:
					gapScore = gapScoreGapFinish
					gapPosj = j
			# If seq2 character is not a C, then use the total gap score assuming the A is also gapped
			else:
				gapScore = gapScoreGapFinish
				gapPosj = j
		else: # AType, no gapping allowed
			gapScore = None
			gapPosi = None
			gapPosj = None
		dirScoreM[i,j] = gapScore,isExtend
		return [gapScore, gapPosi, gapPosj]

	# Execute alignment
	def _aligner(self):
		l1, l2 = len(self.seq1.seq), len(self.seq2.seq)	
		self.scoreMat = numpy.zeros((l1+1, l2+1)) # create matrix for storing counts
		self.directionMat = numpy.zeros((l1+1, l2+1))
		# Each position contains a 2-tuple of the respective gap score and True if the gap is a continuation, False if it is new
		self.leftMat = numpy.zeros((l1+1, l2+1), dtype=('f16,b1')) 
		self.upMat = numpy.zeros((l1+1, l2+1), dtype=('f16,b1'))
		self.scoreMat[0][0] = 0
		for i in range(1, l1 + 1): # set each row by the desired gap
			self.scoreMat[i][0] = self.costs['gap'] * i + self.costs['gapopen']
			self.leftMat[i][0][0] = None
			self.upMat[i][0][0] = None
		for j in range(1, l2 + 1): # set each column by the desired gap
			self.scoreMat[0][j] = self.costs['gap'] * j + self.costs['gapopen']
			self.leftMat[0][j][0] = None
			self.upMat[0][j][0] = None
		for i in range(1, l1+1): # per base-pair in sequence 1 ...
			for j in range(1, l2+1): # per base-pair in sequence 2, align them
			
				if (self.nodeTypes[self.seq1.seq[i-1]] == 'C' and self.nodeTypes[self.seq2.seq[j-1]] == 'A') or (self.nodeTypes[self.seq1.seq[i-1]] == 'A' and self.nodeTypes[self.seq2.seq[j-1]] == 'C'):
					score = None # no match if one is a C type and the other is an A type
				elif (self.nodeTypes[self.seq1.seq[i-1]] == 'T') ^ (self.nodeTypes[self.seq2.seq[j-1]] == 'T'):
					score = None # no match if one is a T type and the other is not
				else:
					score = self.scoreMat[i - 1][j - 1] + get_score(self.seq1.seq[i-1], self.seq2.seq[j-1], self.submat)
				
				# Cost for gapping left (over sequence 1)
				left, lefti, leftj = self.calculate_gap(i,j,self.seq1.seq,self.seq2.seq,self.scoreMat,self.directionMat,self.leftMat,self.TADict1,1)
				# Cost for gapping up (over sequence 2)
				up, upj, upi = self.calculate_gap(j,i,self.seq2.seq,self.seq1.seq,self.scoreMat.T,self.directionMat.T,self.upMat.T,self.TADict2,2)
				
				#out(str(i)+' '+str(j)+' match='+str(score)+' left='+str(left)+' up='+str(up))
				if score is not None and (left is None or score >= left) and (up is None or score >= up):
					# Node match is allowed and produces the best score
					self.scoreMat[i][j] = score
					self.directionMat[i][j] = 0
					self.backPos[i,j] = i-1,j-1
					#out('MATCH')
				elif left is not None and (up is None or left >= up):
					# Gapping left is allows and produces the best score
					self.scoreMat[i][j] = left
					self.directionMat[i][j] = 1
					self.backPos[i,j] = lefti,leftj
					#out('LEFT')
				else:
					# Gapping up is allows and produces the best score
					self.scoreMat[i][j] = up
					self.directionMat[i][j] = 2
					self.backPos[i,j] = upi,upj
					#out('UP')
		
		i, j = l1, l2 # for trace-back process 
		keepGapping = 0
		while i > 0 and j > 0: # walk-back to the index [0][0] of the m
			# if score is a gap in sequence 2 (direction is 1), only walk back on i
			if keepGapping == 1 or keepGapping == 0 and self.directionMat[i][j] == 1:
				# Check whether the choice of gapping left requires the leftward position to also gap left
				if self.leftMat[i,j][1]:
					keepGapping = 1
				else:
					keepGapping = 0
					
				# If the node being gapped is a T-node, appropriately gap the entire subtree
				if self.nodeTypes[self.seq1.seq[i-1]] == 'T':
					previ = self.backPos[i,j][0]
					prevj = self.backPos[i,j][1]
					while i > previ+1: # Walk back with gaps until the one position greater than the final position
						self.align1 += self.seq1.seq[i-1]
						self.align2 += '-'
						i -= 1
					if prevj < j: # If prevj < j, then the gap is preceeded by an A-C match, so add that to the alignment
						self.align1 += self.seq1.seq[i-1]
						self.align2 += self.seq2.seq[j-1]
						i -= 1
						j -= 1
					else: # otherwise the A is also gapped
						self.align1 += self.seq1.seq[i-1]
						self.align2 += '-'
						i -= 1
				# otherwise just gap the current node (it will be a C-node)
				else:
					self.align1 += self.seq1.seq[i-1]
					self.align2 += '-'
					i -= 1
			# if score is a gap in sequence 1 (direction is 2), only walk back on j
			elif keepGapping == 2 or keepGapping == 0 and self.directionMat[i][j] == 2:
				# Check whether the choice of gapping up requires the leftward position to also gap up
				if self.upMat[i,j][1]:
					keepGapping = 2
				else:
					keepGapping = 0
					
				if self.nodeTypes[self.seq2.seq[j-1]] == 'T':
					previ = self.backPos[i,j][0]
					prevj = self.backPos[i,j][1]
					while j > prevj+1:
						self.align1 += '-'
						self.align2 += self.seq2.seq[j-1]
						j -= 1
					if previ < i:
						self.align1 += self.seq1.seq[i-1]
						self.align2 += self.seq2.seq[j-1]
						i -= 1
						j -= 1
					else:
						self.align1 += '-'
						self.align2 += self.seq2.seq[j-1]
						j -= 1
				else:
					self.align1 += '-'
					self.align2 += self.seq2.seq[j-1]
					j -= 1
			# if the score is a match, walk-back one index in both i and j
			elif self.directionMat[i][j] == 0:
				keepGapping = 0
				self.align1 += self.seq1.seq[i-1]
				self.align2 += self.seq2.seq[j-1]
				i -= 1
				j -= 1
		
		# walk-back to index 0 for both i and j; either could be reached first
		while i > 0:
			self.align1 += self.seq1.seq[i-1]
			self.align2 += '-'
			i -= 1
		while j > 0:
			self.align1 += '-'
			self.align2 += self.seq2.seq[j-1]
			j -= 1
		# Reverse the alignment strings as they are assembled backwards
		self.align1 = self.align1[::-1]
		self.align2 = self.align2[::-1]

# Helper-function to write a string
def out(s):
	sys.stdout.write(s+'\n')

# Maps the aligned two bases against a user-selected substitution matrix
def get_score(cA, cB, submatrix):
	return submatrix[(cA, cB)]
#	if (cA, cB) in submatrix:
#		return submatrix[(cA, cB)]
#	else:
#		return submatrix[(cB, cA)] # returns score

# Determines whether the given character pair is in the substitution matrix
def has_score(cA, cB, submatrix):
	return (cA, cB) in submatrix or (cB, cA) in submatrix

# Returns the gap cost of the given character
def get_gapcost(char,submatrix):
	return get_score(char,'-',submatrix)

#def get_gapcost(char,submatrix,gapdefault):
#	if has_score('-',char,submatrix):
#		# Returns a character/residue specific gap cost
#		return get_score('-',char,submatrix)
#	else:
#		# Returns a non-specific gap cost
#		return gapdefault

# Reads in a file containing node types (A,C,T) and a string of specific characters of that type
# Format for a given line should be: "<node-type>:<character string>", for example: "C:BRPD", or the simple case: "C:C"
def parse_nodetypes(fname):
	nodeTypes = {}
	for line in open(fname):
		line = line.strip()
		if len(line) == 0:
			break
		else:
			vals = line.split(':')
			nodeTypes[vals[0]] = vals[1]
	return nodeTypes
	
def default_nodetypes():
	return {'A':'A','C':'C','T':'T'}

