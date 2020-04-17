import cobra

# Returns a list the size of genes in cobra model (model) that contains the cobra reactions that would be 
# constrained to zero if gene i is knocked out alone.

def singleKOList(model):
	# initialize
	geneRxnList = [[] for gene in range(len(model.genes))]	
	i = 0

	# iterate through genes
	for gene in model.genes:
		geneName = gene._id
		# iterate through reactions coupled to gene
		geneRxns = gene._reaction
		for reaction in geneRxns:
			rxnGeneRule = reaction._gene_reaction_rule
			ruleLength = len(rxnGeneRule)

			# find indices of gene in reaction-gene rule
			start = rxnGeneRule.find(geneName)
			end = start + len(geneName) - 1  # end gives the last index of the name
			
			# check if gene alone can knockout reaction
			if start == 0:
				if ruleLength == len(geneName) or rxnGeneRule[end+2] == 'o':
					geneRxnList[i].append(reaction)
					continue
			elif end == ( ruleLength - 1):
				if rxnGeneRule[start-2] == 'r':
					geneRxnList[i].append(reaction)
					continue
			else:
				if rxnGeneRule[start-2] == 'r' or rxnGeneRule[end+2] == 'o':
					geneRxnList[i].append(reaction)
					continue
			
		# update counter
		i += 1			
	
	# we are done
	return geneRxnList	
