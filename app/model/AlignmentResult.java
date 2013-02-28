package model;


/**
 * An AlignmentResult wraps the query and target FASTA sequence which
 * were alignment together, producing a corresponding alignment score.
 * */
public class AlignmentResult implements Comparable<AlignmentResult>{

	private FASTASequence query;
	private FASTASequence target;
	private double score;
	private String alignmentString;

	public AlignmentResult(FASTASequence query, FASTASequence target,
			double score, String alignment) {
		this.setAlignmentString(alignment);
		this.setQuery(query);
		this.setScore(score);
		this.setTarget(target);
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
	 * @return the target
	 */
	public FASTASequence getTarget() {
		return target;
	}

	/**
	 * @param target the target to set
	 */
	private void setTarget(FASTASequence target) {
		this.target = target;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	private void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the alignmentString
	 */
	public String getAlignmentString() {
		return alignmentString;
	}

	/**
	 * @param alignmentString the alignmentString to set
	 */
	private void setAlignmentString(String alignmentString) {
		this.alignmentString = alignmentString;
	}

	@Override
	public int compareTo(AlignmentResult o) {
		if (this.getQuery().getSeqID() < o.getQuery().getSeqID()) {
			return -1;
		}
		else if (this.getQuery().getSeqID() > o.getQuery().getSeqID()) {
			return 1;
		}
		else {
			return 0;			
		}
	}
	
	public String toString() {
		return this.getQuery() + "_" + this.getTarget();
	}


}
