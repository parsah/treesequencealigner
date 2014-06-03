''' 
Runs some finer statistical analysis of MSAs. The initial version is meant to give a value
for the weighted about of sequence material between characters in the consensus to represent
the variability of those regions. Consensus characters with no gaps between represent a 
level of variability no greater than that defined by the consensus threshold.
'''

from parameter import ConsensusStatsCommandParser
from buffer import XMLBuildReader, status_message
import pairwise
from msa import MultipleSequenceDriver
from sequence import NeuriteSequence

version = 0.1

def extract_and_run_stats(args):
    msa = XMLBuildReader(args['build']).parse()
    consensus = msa.consensuses[0]
    ungapped = consensus.ungapped_consensus
    alignments = msa.alignments

    contributions = []
    for alignment in alignments:
        contributions.append(1/len(alignment.replace('-','')))

    current_var_sum = 0
    variability_sums = []
    conservation = []
    # Go through each position in the alignment (based on the composite or gapped consensus)
    for position in range(len(msa.composite)):
        # If the current position is a gap, add contributions to the current variability sum
        if consensus.seq[position] == '-':
            # Go through each alignment to check whether it has a character at the current position
            for alignment_num in range(len(alignments)):
                # If it does have a character (not a gap '-'), add its variability contribution
                if not alignments[alignment_num][position] == '-':
                    current_var_sum += contributions[alignment_num]
        # The current position is a consensus character, so store the variability sum preceeding it
        else:
            # if the position within the ungapped consensus is needed, just use the current length of variability_sums
            variability_sums.append(current_var_sum)
            current_var_sum = 0

            # determine the current position's conservation level
            conservation_sum = 0
            for alignment_num in range(len(alignments)):
                # If it does have a character (not a gap '-'), add its conservation contribution
                if not alignments[alignment_num][position] == '-':
                    conservation_sum += 1
            conservation.append(conservation_sum/len(alignments))

    if 'newick' in args.keys():
        status_message('Generating newick string','OK')
        newick_string = generate_consensus_newick(consensus,conservation,variability_sums)
        handle = open(args['newick'],'w')
        handle.write(newick_string)
        handle.close()
    

#TODO - A stack is going to be required for this. Look at Java code to recreate (shouldn't be too hard)
def generate_consensus_newick(consensus,conservation,variability_sums):
    newick_string = ""
    cons_str = consensus.ungapped_consensus
    a_stack = [(1,0)]
    print(cons_str)
    for position in range(len(cons_str)):
        #print(str(position)+' '+cons_str[position])
        #print(str(a_stack))
        if cons_str[position] == 'A':
            a_stack.append((0,0))
            node_string = ')'
        elif cons_str[position] == 'C':
            node_string = ',:1)'
            a_val = a_stack.pop()
            a_stack.append((a_val[0],a_val[1]+1))
        elif cons_str[position] == 'T':
            node_string = '(:1,:1)'
            not_done = 1
            while not_done and len(a_stack) > 0:
                a_val = a_stack.pop()
                if a_val[0] == 0:
                    node_string = ','+a_val[1]*'('+node_string
                    a_stack.append((1,1))
                    not_done = 0
                else:
                    node_string = a_val[1]*'('+node_string
            
        newick_string = node_string + cons_str[position] + '-' + str(round(conservation[position],3)) + ':' + str(round(variability_sums[position],3)) + newick_string

    newick_string += ';'
    return newick_string

if __name__ == '__main__':
    try:
        args = ConsensusStatsCommandParser().parse_args()
#        ConsensusStatsArgumentValidator(args) # test all arguments are correct

        print('capellini - v.' + str(version) + '\n=============')

        extract_and_run_stats(args)

        status_message('Consensus statistics computation complete ', 'OK')

    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
