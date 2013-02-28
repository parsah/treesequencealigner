package model;

public class NeedlemanWunsch {

	private Matrix matrix;
	private FASTASequence query;
	private FASTASequence target;

	public NeedlemanWunsch(FASTASequence query, FASTASequence target) {
		this.setMatrix(null);
		this.setQuery(query);
		this.setTarget(target);
		this.build();
	}

	/**
	 * An initial matrix must have initial gap extension costs preset. These
	 * costs are provided at runtime and facilitate the crux of this algorithm.d
	 * */
	private void build() {
		double[][] matrix = new 
				double[this.getQuery().length() + 1]
						[this.getTarget().length() + 1];
		this.setMatrix(new Matrix(matrix)); // initial matrix

		for (int i = 0; i <= this.getQuery().length(); i++) {
			for (int j = 0; j <= this.getTarget().length(); j++) {
				if (i == 0) {
					this.getMatrix().getData()[i][j] = j; // extension costs
				} else if (j == 0) {
					this.getMatrix().getData()[i][j] = i; // extension costs
				} else {
					this.getMatrix().getData()[i][j] = 0;
				}
			}
		}
	}

	public void align() {
		for (int i = 1; i <= this.getQuery().length(); i++) {
			for (int j = 1; j <= this.getTarget().length(); j++) {
				double scoreDiag = this.getMatrix().getData()[i-1][j-1] + 
						getWeight(i, j);
				double scoreLeft = this.getMatrix().getData()[i][j-1] - 1;
				double scoreUp = this.getMatrix().getData()[i-1][j] - 1;
				double score = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
				this.getMatrix().getData()[i][j] = score;
			}
		}
	}

	public void traceBack() {
		String mAlignmentSeqA = "";
		String mAlignmentSeqB = "";

		int i = this.getQuery().length();
		int j = this.getTarget().length();
		while (i > 0 && j > 0) {                        
			if (this.getMatrix().getData()[i][j] == this.getMatrix().getData()[i-1][j-1] + this.getWeight(i, j)) {                          
				mAlignmentSeqA += this.getQuery().getSequence().charAt(i-1);
				mAlignmentSeqB += this.getTarget().getSequence().charAt(j-1);
				i--;
				j--;                            
				continue;
			} else if (this.getMatrix().getData()[i][j] == this.getMatrix().getData()[i][j-1] - 1) {
				mAlignmentSeqA += "-";
				mAlignmentSeqB += this.getTarget().getSequence().charAt(j-1);
				j--;
				continue;
			} else {
				mAlignmentSeqA += this.getQuery().getSequence().charAt(i-1);
				mAlignmentSeqB += "-";
				i--;
				continue;
			}
		}
		mAlignmentSeqA = new StringBuffer(mAlignmentSeqA).reverse().toString();
		mAlignmentSeqB = new StringBuffer(mAlignmentSeqB).reverse().toString();
		System.out.println(mAlignmentSeqA);
		System.out.println(mAlignmentSeqB);
	}


	private int getWeight(int i, int j) {
		char charQuery = this.getQuery().getSequence().charAt(i-1);
		char charTarget = this.getTarget().getSequence().charAt(j-1);
		if (charQuery == charTarget) {
			return 1;
		} else {
			return -1;
		}
	}



	/**
	 * Given a populated score-matrix, the score of the alignment is found
	 * at the bottom right-most index. Based on initial gap open and extension
	 * costs, this total score can be either positive or negative.
	 * @return alignment score.
	 * */
	public double getScore() {
		int matrixHeight = this.getMatrix().getHeight();
		int matrixWidth = this.getMatrix().getWidth();
		return this.getMatrix().getData()[matrixHeight-1][matrixWidth-1];
	}


	/**
	 * @return the matrix
	 */
	public Matrix getMatrix() {
		return matrix;
	}

	/**
	 * @param matrix the matrix to set
	 */
	private void setMatrix(Matrix matrix) {
		this.matrix = matrix;
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

}
