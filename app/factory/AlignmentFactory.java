package factory;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import runtime.Presenter;

import model.AlignmentResult;
import model.FASTASequence;

/**
 * An AlignmentFactory object performs all-to-all pairwise sequence
 * alignment using a thread-pool. Such computation is made possible given
 * a user-provided FASTA file and a positive number of desired worker
 * threads. The end result of such work is an AlignmentResult object.
 * This object encapsulates a number of behaviors and states such as
 * the query and target sequence, alignment score, and traceback.
 * */
public class AlignmentFactory implements Callable<AlignmentResult> {
	
	@Override
	public AlignmentResult call() throws Exception {
		return null;
	}
	
	/**
	 * Begin all-to-all sequence alignment whereby a thread-worker pool is
	 * defined apriori.
	 * */
	public static void exhaustiveAlignment(List<FASTASequence> sequences) {
		int numThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService service = Executors.newFixedThreadPool(numThreads);
		Presenter.out(numThreads + "x worker threads created [OK]", true);
		for (FASTASequence seq: sequences) {
			System.out.println(seq);
		}
		
		
		service.shutdown();
		Presenter.out("Analysis complete. Thread-pool closed. [OK]", true);
		
	}

}
