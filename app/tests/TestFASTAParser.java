package tests;

import static org.junit.Assert.*;

import org.junit.Test;

public class TestFASTAParser {
	
	// If parsing was not possible, there must be no sequences
	@Test
	public void testNoSeqsIfFailed() {
	}
	
	// If an invalid file is provided, no parsing can be performed
	@Test
	public void testNoParsingIfInvalidFile() {
	}
	
	// If parsing is a success, there must be atleast 1 sequence
	@Test
	public void testManySeqsIfSuccess() {
	}
	
	// When parsing is complete and successful, both success state must
	// be true and so too the presence of many sequences
	@Test
	public void testSuccessImpliesManySeqs() {
	}

}
