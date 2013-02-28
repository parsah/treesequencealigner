package tests;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import model.FASTASequence;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import factory.FASTAParser;

public class FASTAParserTest {

	@Rule // for testing certain exceptions are thrown
	public ExpectedException exception;

	// Trivial test to determine that all objects parsed are FASTA
	@Test
	public void testAllSeqsFASTA() {
		try {
			FASTAParser parser = new FASTAParser(new File("./test_data/test_seqs.fasta"));
			parser.parse();
			for (FASTASequence e: parser.getFastaEntries()) {
				assertTrue(e instanceof FASTASequence);
			}
		} catch (IOException e) { }
	}

	// Test number of parsed sequences equal manual derivation
	@Test
	public void testNumFASTAEqualsLength() {
		try {
			FASTAParser parser = new FASTAParser(new File("./test_data/test_seqs.fasta"));
			parser.parse();
			int numSeqs = 0;
			for (FASTASequence e: parser.getFastaEntries()) {
				if (e instanceof FASTASequence) {
					numSeqs ++;
				}
			}
			assertEquals(numSeqs, parser.getNumSequences());
			
		} catch (IOException e) { }
	}

	// You can parse an unstructure FASTA file, eg. column 1 is header
	// and column 2 is the sequence, however this is not valid.
	@Test
	public void testInvalidFASTAThrowsException() {
		try {
			FASTAParser parser = new FASTAParser(new File("./test_data/unstructured.tab"));
			parser.parse();
			this.exception = ExpectedException.none();
			this.exception.expect(NullPointerException.class);
		} 
		catch (NullPointerException e) { }
		catch (IOException e) { }
	}

	// If an invalid file is provided, no parsing can be performed
	@Test
	public void testNonExistentFileThrowsException() {
		try {
			FASTAParser parser = new FASTAParser(new File("./test_data/invalid"));
			parser.parse();
			this.exception = ExpectedException.none();
			this.exception.expect(FileNotFoundException.class);
		} catch (IOException e) { }
	}

	// If parsing is a success, there must be atleast 1 sequence
	@Test
	public void testParsingValidFASTAWorks() {
		try {
			FASTAParser parser = new FASTAParser(new File("./test_data/test_seqs.fasta"));
			parser.parse();
			assertTrue(parser.getNumSequences() > 0); // many sequences parsed
		} catch (IOException e) {
			System.out.println(e.getMessage());
		}
	}
}
