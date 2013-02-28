package model;

public class FASTASequence {
	private String sequence; // FASTA sequence
	private String header; // FASTA header
	private int seqID; // integer to help ID the sequence

	// Creates a bare-bones FASTA object
	public FASTASequence() {
		this.setHeader(""); // initially, no header
		this.setSequence(""); // initially, no sequence
		this.setSeqID(0);
	}

	/**
	 * Pertains to the actual FASTA sequence object
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * Set an actual sequence string
	 * @param sequence the sequence to set
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	/**
	 * Return the FASTA header which identifies this entry
	 * @return the header
	 */
	public String getHeader() {
		return header;
	}

	/**
	 * Set the FASTA header
	 * @param header the header to set
	 */
	public void setHeader(String header) {
		this.header = header;
	}

	/**
	 * Returns the length of the FASTA sequence entry
	 * @return sequence length
	 * */
	public int length() {
		return this.getSequence().length();
	}
	
	@Override
	public String toString() {
		return "\"" + this.getHeader() + "\"";
	}
	
	/**
	 * Represents a behavior to model object hashes
	 * @return object hash solely based on the header's hashcode
	 * */
	@Override
	public int hashCode() {
		return this.getHeader().hashCode();
	}
	
	/**
	 * Equality function to determine if two FASTA objects are
	 * the same. Such equality is based solely on their respective
	 * headers only.
	 * @return boolean to represent equality
	 * */
	@Override
	public boolean equals( Object obj ) {
		FASTASequence fasta = (FASTASequence)obj;
		if( fasta.getHeader() == this.getHeader() )
			return true;
		else {
			return false;
		}
	}

	/**
	 * @return the sequence ID
	 */
	public int getSeqID() {
		return seqID;
	}

	/**
	 * @param seqID the sequenceNum to set
	 */
	public void setSeqID(int seqID) {
		this.seqID = seqID;
	}

}
