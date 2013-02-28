package factory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import model.FASTASequence;
import runtime.Presenter;

/**
 * A FASTAParser class parses a user-provided FASTA file and produces
 * a list of FASTAEntry objects
 * */
public class FASTAParser {
	private ArrayList<FASTASequence> fastaEntries; // Parsed FASTA objects
	private File fastaFile; // input fasta file

	/**
	 * Create a bare-bone FASTA-parser object which ultimately produces
	 * a list of parsed FASTA objects
	 * @throws IOException 
	 * */
	public FASTAParser(File input) throws IOException {
		this.setFastaEntries(new ArrayList<FASTASequence>()); // empty list
		this.setFastaFile(input); // set the user-provided input
	}

	/**
	 * Parse a user-provided FASTA file and save resultant FASTA objects
	 * in a list for further analyses.
	 * @throws IOException 
	 * */
	public void parse() throws IOException {
		if (this.getFastaFile().exists()) {
			Presenter.out("Parsing " + this.getFastaFile(), true);
			BufferedReader reader = new BufferedReader( // read FASTA file
					new FileReader(this.getFastaFile()));
			String eachLine = ""; // represents each line in the FASTA
			FASTASequence entry = new FASTASequence(); // current entry
			while ((eachLine = reader.readLine()) != null) {
				if (eachLine.startsWith(">")) {
					entry = new FASTASequence();
					entry.setHeader(eachLine.substring(1, eachLine.length()));
					this.getFastaEntries().add(entry);
				}
				else {
					entry.setSequence(entry.getSequence() + eachLine);						
				}
			}
			reader.close(); // close input file reader
			this.validate(); // ensure objects have sequences and are sound
		}
		else {
			throw new FileNotFoundException("Invalid FASTA [ERROR]");
		}
	}

	/**
	 * Verify the parser processed a sound FASTA file.
	 * */
	private void validate() {
		if (this.getNumSequences() == 0) { // if no seqs parsed, exit
			throw new NullPointerException("Input is not of FASTA format [ERROR]");
		}
		else { // iterate over each FASTA entry and make sure it has sequence
			boolean is_parser_valid = true;
			for (FASTASequence e: this.getFastaEntries()) {
				if (e.length() == 0) {
					is_parser_valid = false; // parser invalid if no sequence
				}
			}
			if (is_parser_valid == false) { // if a FASTA entry is invalid
				throw new NullPointerException("Sequence missing for FASTA entry [ERROR]");
			}
			else { // if all FASTA entries are good, proceed
				Presenter.out(this.getNumSequences() +" sequences parsed [OK]", true);
			}
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
	public ArrayList<FASTASequence> getFastaEntries() {
		return fastaEntries;
	}

	/**
	 * Sets a list of FASTA objects
	 * @param fastaEntries the fastaEntries to set
	 */
	public void setFastaEntries(ArrayList<FASTASequence> fastaEntries) {
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
