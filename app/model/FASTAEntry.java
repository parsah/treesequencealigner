package model;

public class FASTAEntry {
	private String sequence; // FASTA sequence
	private String header; // FASTA header

	// Creates a bare-bones FASTA object
	public FASTAEntry() {
		this.setHeader(""); // initially, no header
		this.setSequence(""); // initially, no sequence
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
	public int getSequenceLength() {
		return this.getSequence().length();
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
		FASTAEntry fasta = (FASTAEntry)obj;
		if( fasta.getHeader() == this.getHeader() )
			return true;
		else {
			return false;
		}
	}

}
