package factory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import model.FASTAEntry;
import runtime.Presenter;

/**
 * A FASTAParser class parses a user-provided FASTA file and produces
 * a list of FASTAEntry objects
 * */
public class FASTAParser {
	private ArrayList<FASTAEntry> fastaEntries; // Parsed FASTA objects
	private File fastaFile; // input fasta file

	/**
	 * Create a bare-bone FASTA-parser object which ultimately produces
	 * a list of parsed FASTA objects
	 * */
	public FASTAParser(File input) {
		this.setFastaEntries(new ArrayList<FASTAEntry>()); // empty list
		this.setFastaFile(input); // set the user-provided input
	}

	/**
	 * Parse a user-provided FASTA file and save resultant FASTA objects
	 * in a list for further analyses.
	 * @throws IOException 
	 * */
	public void parse() throws IOException {
		if (this.getFastaFile().isFile()) {
			Presenter.out("Parsing " + this.getFastaFile(), true);
			BufferedReader reader = new BufferedReader( // read FASTA file
					new FileReader(this.getFastaFile()));
			String eachLine = ""; // represents each line in the FASTA
			FASTAEntry entry = new FASTAEntry(); // current entry
			while ((eachLine = reader.readLine()) != null) {
				if (eachLine.startsWith(">")) {
					entry = new FASTAEntry();
					entry.setHeader(eachLine.substring(1, eachLine.length()));
					this.getFastaEntries().add(entry);
				}
				else {
					entry.setSequence(entry.getSequence() + eachLine);						
				}
			}
			reader.close(); // close input file reader
			this.verifyParser(); // verify parsing was a success
		}
		else {
			throw new FileNotFoundException("Invalid FASTA [ERROR]");
		}
	}
	
	/**
	 * Verify the parser processed a valid FASTA file.
	 * */
	private void verifyParser() {
		if (this.getNumSequences() == 0) { // if no seqs parsed, exit
			throw new NullPointerException("Input file is not FASTA [ERROR]");
		}
		else { // display the number of sequences parsed
			Presenter.out(this.getNumSequences() +" sequences parsed [OK]", true);
		}
	}

	/**
	 * Returns the number of parsed FASTA objects
	 * @return count as-to the number of parsed FASTAEntry objects
	 * */
	public int getNumSequences() {
		return this.getFastaEntries().size();
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
