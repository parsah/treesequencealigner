package tests;

import static org.junit.Assert.*;

import model.FASTASequence;

import org.junit.Test;

public class FASTAEntryTest {
	
	@Test
	public void testFASTALength() {
		// Test a FASTA entry is of the right size
		FASTASequence entry = new FASTASequence();
		entry.setSequence("ATCTTACTTCATCATCAT");
		assertEquals(entry.length(), 18);
	}
	
	@Test
	public void testNoFASTASequence() {
		// No supplied sequence yields null
		FASTASequence entry = new FASTASequence();
		assertEquals(entry.getSequence(), "");
	}
	
	@Test
	public void testTwoSameHeadersEquality() {
		// Test whether two sequences with same headers are equal
		FASTASequence fastaA = new FASTASequence();
		fastaA.setHeader("A"); // create a random fasta header
		FASTASequence fastaB = new FASTASequence();
		fastaB.setHeader("A"); // create another fasta header
		assertEquals(fastaA, fastaB); // equality is based on header ID
	}
	
	@Test
	public void testTwoDifferentHeadersEquality() {
		// Test whether two sequences with different headers are not equal
		FASTASequence fastaA = new FASTASequence();
		fastaA.setHeader("A"); // create a random fasta header
		FASTASequence fastaB = new FASTASequence();
		fastaB.setHeader("AA"); // create another fasta header
		assertNotSame(fastaA, fastaB); // equality is based on header ID
	}

}
