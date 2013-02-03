package tests;

import static org.junit.Assert.*;

import model.FASTAEntry;

import org.junit.Test;

public class TestFASTAEntry {
	
	@Test
	public void testFASTALength() {
		// Test a FASTA entry is of the right size
		FASTAEntry entry = new FASTAEntry();
		entry.setSequence("ATCTTACTTCATCATCAT");
		assertEquals(entry.getSequenceLength(), 18);
	}
	
	@Test
	public void testNoFASTASequence() {
		// No supplied sequence yields null
		FASTAEntry entry = new FASTAEntry();
		assertNull(entry.getSequence());
	}
	
	@Test
	public void testTwoSameHeadersEquality() {
		// Test whether two sequences with same headers are equal
		FASTAEntry fastaA = new FASTAEntry();
		fastaA.setHeader("A"); // create a random fasta header
		FASTAEntry fastaB = new FASTAEntry();
		fastaB.setHeader("A"); // create another fasta header
		assertEquals(fastaA, fastaB); // equality is based on header ID
	}
	
	@Test
	public void testTwoDifferentHeadersEquality() {
		// Test whether two sequences with different headers are not equal
		FASTAEntry fastaA = new FASTAEntry();
		fastaA.setHeader("A"); // create a random fasta header
		FASTAEntry fastaB = new FASTAEntry();
		fastaB.setHeader("AA"); // create another fasta header
		assertNotSame(fastaA, fastaB); // equality is based on header ID
	}

}
