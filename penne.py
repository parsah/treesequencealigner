''' 
Enables the ability to analyze multiple sequence alignment (MSA) builds and
identify over-represented domains within such analyses.
'''

from parameter import DomainArgumentValidator, DomainCommandParser, InputWrapperState
from domain import DomainSetBuilder, DomainAbundanceBuilder, DomainPrettyPrinter, DomainExtractionDriver, merge_counts
from buffer import XMLBuildReader, status_message
import pairwise
from msa import MultipleSequenceDriver, ConsensusFilterFactory
from pairwise import PairwiseDriver
#from scipy import cluster
from sequence import NeuriteSequence

version = 0.2 
# Updates: Run MSA (if necessary); Wrapper to run MSA multiple times (shuffled or
# with multiple subsets of the given fasta files), and clustering of domains
max_cluster_dist = 2

def extract_and_analyze_domains(targets,baselines,input_state):
    domain_set = set() # domains from targets and baselines
    target_domains = {}
    baseline_domains = {}
    args = input_state.get_args()
    node_types = input_state.node_types
    subsmat = input_state.subsmat
    num_runs = args['runs']

    output_file = args['o']
#    output_base = args['o']
#    if output_base is not None:
#        consensus_file = output_base + ".consensuses.txt"
#        domain_file = output_base + ".domains.txt"

    ded = DomainExtractionDriver(targets,baselines,node_types,subsmat,num_runs,input_state)
    ded.start()
    domain_set = ded.domain_set
    target_domains = ded.target_domains
    baseline_domains = ded.baseline_domains
    target_consensuses = ded.target_consensuses
    baseline_consensuses = ded.baseline_consensuses

    # Cluster domains
    if args['cluster']:
        combined_occurrences = merge_counts(target_domains,baseline_domains)
        #domain_list = list(domain_set)
        #domain_seeds, domain_clusters = cluster_domains(combined_occurrences,node_types,subsmat,args['max_cluster_dist'],args)
        domain_clusters = cluster_domains(combined_occurrences,node_types,subsmat,args)

        # Compile counts for each cluster from the individual domains
        target_clust_counts = {}
        baseline_clust_counts = {}
        #for clust_num in range(len(domain_clusters)):
        for seed in domain_clusters:
            target_count = 0
            baseline_count = 0
            for domain in domain_clusters[seed]:
                if domain in target_domains.keys():
                    target_count += target_domains[domain]
                if domain in baseline_domains.keys():
                    baseline_count += baseline_domains[domain]
            target_clust_counts[seed] = target_count
            baseline_clust_counts[seed] = baseline_count
            #target_clust_counts["DC"+str(clust_num+1)] = target_count
            #baseline_clust_counts["DC"+str(clust_num+1)] = baseline_count
            #target_clust_counts.append(target_count)
            #baseline_clust_counts.append(baseline_count)
        target_counts = target_clust_counts
        baseline_counts = baseline_clust_counts
    else:
        target_counts = target_domains
        baseline_counts = baseline_domains

    # Calculate probability of target being different from baseline for each domain (cluster)
    dab = DomainAbundanceBuilder(target_counts,baseline_counts)
    domain_matrices = dab.build();

    # Print out domain cluster seeds, their domains, and their p-value
    domain_tuples = []
    for domain_mat in domain_matrices:
        cluster_name = domain_mat.name
        #i = int(cluster_name.replace("DC",""))-1
        #cluster_seed = domain_seeds[i]
        cluster_seed = cluster_name
        target_counts = target_clust_counts[cluster_name]
        baseline_counts = baseline_clust_counts[cluster_name]
        pval = domain_mat.get_hypergeometric_prob()
#        domain_tuples.append((i, cluster_name, cluster_seed, pval, target_counts, baseline_counts))
        domain_tuples.append((cluster_name, pval, target_counts, baseline_counts))

    sorted_tuples = sorted(domain_tuples, key=lambda domain: domain[1])
    if output_file is not None:
        handle = open(output_file,'w')
        handle.write("Cluster Seed,pVal,Target Counts,Baseline Counts,Domains\n")
        for tup in sorted_tuples:
            handle.write(tup[0]+','+str(round(tup[1], 6))+','+str(tup[2])+','+str(tup[3])+', "'+str(domain_clusters[tup[0]])+'"\n')
            handle.flush()
        handle.close()

        handle = open(output_file+".consensuses.txt",'w')
        handle.write('TARGETS:\n')
        for consensus in target_consensuses:
            handle.write(consensus.seq.replace('-','')+'\n')
        handle.write('\nBASELINES:\n')
        for consensus in baseline_consensuses:
            handle.write(consensus.seq.replace('-','')+'\n')
        handle.close()
    else:
        for tup in sorted_tuples:
            tuple_str = tup[0]+','+str(round(tup[1], 6))+','+str(tup[2])+','+str(tup[3])
            print(tuple_str)
            print(domain_clusters[tup[0]])
            print()
            print('TARGET CONSENSUSES:\n')
            for consensus in target_consensuses:
                print(consensus.seq.replace('-',''))
            print('\nBASELINES:')
            for consensus in baseline_consensuses:
                print(consensus.seq.replace('-',''))

    
    #dpp = DomainPrettyPrinter(domains=domain_matrices,pval=args['p'],out=args['o'])
    #dpp.display()

    # Run statistics on groups and produce output
    #for clust_num in range(len(domain_clusters)):
    #    print("Domain Cluster: "+str(domain_clusters[clust_num]))
    #    dab = DomainAbundanceBuilder(target_clust_counts[clust_num],baseline_clust_counts[clust_num])
    #    domain_matrices = dab.build();
    #    dpp = DomainPrettyPrinter(domains=domain_matrices,pval=args['p'],out=args['o'])
    #    dpp.display()

def cluster_domains(domain_occurrences,node_types,subsmat,args):
    '''
    A class that finds pairwise edit distances between domains and clusters them, does so for each domain length.
    Returns a domain-seed keyed dictionary of domain lists (clusters).
    '''
    # Set the minimum size at which domains should be clustered
    min_domain_size = 7 # maybe should be 8
    print("Clustering domains")
    # Run pairwise alignment
    domain_args = {'f':None,'f2':None,'a':None,'subsmat':subsmat,'gap':-1,'gapopen':0,'matrix':None,'custom':None,'o':None,'n':args['n'],'node_types':None}
    input_state = InputWrapperState(domain_args)
    input_state.subsmat = subsmat
    domains = list([domain for domain in domain_occurrences.keys()])

    # Cluster results via average distance criterion
    #linkages = cluster.hierarchy.average(distance_mat)

    ### Split domains up by size first, then cluster those that are large enough
    # Determine max length
    min_length = min(len(domain) for domain in domains)
    max_length = max(len(domain) for domain in domains)

    all_clusters = {}

    # Move short domains into their own clusters
    for length in range(min_length,min_domain_size):
        domains_sub = list([domain for domain in domains if len(domain) == length])
        for domain in domains_sub:
            all_clusters[domain] = [domain]

    # Cluster domains that are long enough for each size
    for length in range(min_domain_size,max_length+1):
        domains_sub = list([domain for domain in domains if len(domain) == length])

        #domains_ns = list([NeuriteSequence("D"+str(i),domains[i]) for i in range(len(domains))])
        domains_ns = list([NeuriteSequence(domains_sub[i],domains_sub[i]) for i in range(len(domains_sub))])
        domain_id_map = {domains_sub[i]:i for i in range(len(domains_sub))}
        driver = PairwiseDriver(domains_ns, domains_ns, input_state, store_pairwise=True, score_type='num_gaps')
        driver.start()
        distance_mat = driver.get_score_matrix()

        '''
        Determine clusters from hierarchy ???
        Use abundance sort and then UCLUST 
        - First and subsequent sequences are seeds unless they are close enough to an existing seed, 
          in which case they are merged with the seed's cluster.
        - Domains shorter than min_domain_size go in their own separate clusters.
        - Domains that are exact sub- or super-sequences of a domain already in a cluster are excluded
        '''
        # First sort domains by occurrence frequency, requires creating tuples from domain and occurrence count
        domain_tuples = []
        for domain in domains_sub:
            domain_tuples.append((domain,domain_occurrences[domain]))
        sorted_tuples = sorted(domain_tuples, key=lambda domain_tup: (domain_tup[1],len(domain_tup[0])), reverse=True)
        sorted_domains = list([tup[0] for tup in sorted_tuples])
    
        if len(sorted_domains) > 0:
            # UCLUST implementation
            separates = {}
            clusters = {sorted_domains[0]:list([sorted_domains[0]])}
            for domain in sorted_domains[1:]:
                domain_pos = domain_id_map[domain]
                closest_cluster = -1,1000
                if len(domain) < min_domain_size:
                    separates[domain] = list([domain])
                else:
                    for seed in clusters:
                        seed_pos = domain_id_map[seed]
                        dist = max(distance_mat[domain_pos][seed_pos],distance_mat[seed_pos][domain_pos])
                        if dist <= max_cluster_dist and dist < closest_cluster[1]:
                            closest_cluster = seed,dist

                    if closest_cluster[0] == -1:
                        # If not close enough to any existing seed, add new cluster
                        clusters[domain] = list([domain])
                    else:
                        # Otherwise add it to the cluster it's closest to
                        clusters[closest_cluster[0]].append(domain)
           
            # Copy new clusters to all_clusters
            for seed in clusters:
                all_clusters[seed] = clusters[seed]

    return all_clusters


if __name__ == '__main__':
    try:
        args = DomainCommandParser().parse_args()
        DomainArgumentValidator(args) # test all arguments are correct

        if args['mode'] == 'single':
            print('penne - v.' + str(version) + '\n=============')
            cons_query = XMLBuildReader(args['query']).parse()
            cons_baseline = XMLBuildReader(args['baseline']).parse()
        
            # next, yield domains for both query and baseline datasets. 
            dsb_query = DomainSetBuilder(win=args['win'], max_gap=args['max_g'], 
                         is_enum=args['enumerate'], consensus=cons_query,
                         is_strip=args['strip'])
            dsb_baseline = DomainSetBuilder(win=args['win'], max_gap=args['max_g'], 
                         is_enum=args['enumerate'], consensus=cons_baseline,
                         is_strip=args['strip'])
            domains_query = dsb_query.build() # build abundance counts
            domains_baseline = dsb_baseline.build()
            status_message('Identifying domains', 'OK')
            db = DomainAbundanceBuilder(query=domains_query, baseline=domains_baseline)
            domains = db.build() # build contingency matrices 
            dpp = DomainPrettyPrinter(domains = domains, pval = args['p'],
                                  out=args['o'])
            dpp.display() # pretty-print domains
            status_message('Domain over-representation computation complete ', 'OK')
        else:
            args.update({'f':args['query'],'f2':args['baseline'],'a':None})
            input_state = InputWrapperState(args)
            #input_state.assign_matrix() # parse in-built or custom matrix
            targets = input_state.parse_fasta(input_state.fname)
            baselines = input_state.parse_fasta(input_state.fname2)
            extract_and_analyze_domains(targets,baselines,input_state)
            status_message('Domain analysis via multiple runs complete ', 'OK')

    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
