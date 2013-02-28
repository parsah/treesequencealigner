package factory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import model.AlignmentResult;
import model.FASTASequence;
import model.NeedlemanWunsch;
import runtime.Presenter;

/**
 * An AlignmentFactory object performs all-to-all pairwise sequence
 * alignment using a thread-pool. Such computation is made possible given
 * a user-provided FASTA file and a positive number of desired worker
 * threads. The end result of such work is an AlignmentResult object.
 * This object encapsulates a number of behaviors and states such as
 * the query and target sequence, alignment score, and traceback.
 * */
public class AlignmentFactory implements Callable< List<AlignmentResult>> {

	private FASTASequence query;
	private List<FASTASequence> targets;

	/**
	 * Each instance of AlignmentFactory is such that a specific query object
	 * is iteratively mapped against all other targets. This all-to-all
	 * approach is the basis behind having a concurrent implementation, 
	 * considering the fact that dataset sizes can only increase.
	 * @param query Query FASTA object
	 * @param targets List of FASTA objects the query will be mapped against.
	 * */
	public AlignmentFactory(FASTASequence query, List<FASTASequence> targets) {
		this.setQuery(query);
		this.setTargets(targets);
	}

	@Override
	public List<AlignmentResult> call() throws Exception {
		// Iterate over each target and align the query to this target.
		List<AlignmentResult> results = new ArrayList<AlignmentResult>();
		for (FASTASequence target: this.getTargets()) {
			NeedlemanWunsch nw = new NeedlemanWunsch(this.getQuery(), target);
			nw.align();
			// wrap results from the alignment
			results.add(new AlignmentResult(this.getQuery(), target, 
					nw.getScore(), ""));
		}
		return results; // exhaustive results for the respective query.
	}

	/**
	 * Begin all-to-all sequence alignment whereby a thread-worker pool is
	 * defined apriori.
	 * */
	public static void exhaustiveAlignment(List<FASTASequence> seqs) {
		int nThreads = Runtime.getRuntime().availableProcessors();
		Presenter.out(nThreads + "x worker threads created [OK]", true);
		ExecutorService service = Executors.newFixedThreadPool(2);
		try {
			// create a list of row-specific alignments.
			int counter = 1;
			List<Future<List<AlignmentResult>>>  futures = 
					new ArrayList< Future<List<AlignmentResult>>>();
			for (FASTASequence seq: seqs) {
				Callable<List<AlignmentResult>> out = new AlignmentFactory(seq, seqs);
				futures.add(service.submit(out)); // add output once processed.
			}

			for (Future<List<AlignmentResult>> future: futures) {
				
				List<AlignmentResult> result = future.get();
				Collections.sort(result);
				//				System.out.println((double)counter % seqs.size());
				if (counter % 3 == 0) {
					System.out.printf("%d / %d alignments complete\n", 
							counter, seqs.size());
					//				}
				}
				++counter;
//				System.out.println((int)((double)counter % seqs.size()));
			}

			Presenter.out("Analysis complete. Thread-pool closed. [OK]", true);
		} catch (ExecutionException e) {
			System.out.println("Cannot retrieve result of concurrent task.");
		} catch (InterruptedException e) {
			System.out.println("Thread-pool interrupted. Analysis aborted.");
		} finally {
			service.shutdown();
		}

	}

	/**
	 * @return the query
	 */
	public FASTASequence getQuery() {
		return query;
	}

	/**
	 * @param query the query to set
	 */
	private void setQuery(FASTASequence query) {
		this.query = query;
	}

	/**
	 * @return the targets
	 */
	public List<FASTASequence> getTargets() {
		return targets;
	}

	/**
	 * @param targets the targets to set
	 */
	private void setTargets(List<FASTASequence> targets) {
		this.targets = targets;
	}

}
