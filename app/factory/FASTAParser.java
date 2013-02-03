package factory;

import java.io.File;
import java.util.ArrayList;

import model.FASTAEntry;

/**
 * A FASTAParser class parses a user-provided FASTA file and produces
 * a list of FASTAEntry objects
 * */
public class FASTAParser {
	private boolean successful; // whether parsing was successful.
	private ArrayList<FASTAEntry> fastaEntries; // Parsed FASTA objects
	private File fastaFile; // input fasta file
	
	/**
	 * Create a bare-bone FASTA-parser object which ultimately produces
	 * a list of parsed FASTA objects
	 * */
	public FASTAParser(File input) {
		this.setSuccessful(false); // no parsing has been performed
		this.setFastaEntries(new ArrayList<FASTAEntry>()); // empty list
		this.setFastaFile(input); // set the user-provided input
	}
	
	/**
	 * Parse a user-provided FASTA file and save resultant FASTA objects
	 * in a list for further analyses.
	 * */
	public void parse() {
		
	}
	
	/**
	 * Returns the number of parsed FASTA objects
	 * @return count as-to the number of parsed FASTAEntry objects
	 * */
	public int getNumSequences() {
		return this.getFastaEntries().size();
	}
	
	/**
	 * Returns whether FASTA parsing was a success
	 * @return the successful
	 */
	public boolean isSuccessful() {
		return successful;
	}
	
	/**
	 * Specify whether FASTA parsing was successful.
	 * @param successful the successful to set
	 */
	public void setSuccessful(boolean successful) {
		this.successful = successful;
	}
	
	/**
	 * Return list of FASTA objects
	 * @return the fastaEntries
	 */
	public ArrayList<FASTAEntry> getFastaEntries() {
		return fastaEntries;
	}
	
	/**
	 * Sets a list of FASTA objects
	 * @param fastaEntries the fastaEntries to set
	 */
	public void setFastaEntries(ArrayList<FASTAEntry> fastaEntries) {
		this.fastaEntries = fastaEntries;
	}

	/**
	 * Returns a reference to the user-provided FASTA file
	 * @return the fastaFile
	 */
	public File getFastaFile() {
		return fastaFile;
	}

	/**
	 * Specify a user-provided input FASTA file to parse
	 * @param fastaFile the fastaFile to set
	 */
	public void setFastaFile(File fastaFile) {
		this.fastaFile = fastaFile;
	}
}
