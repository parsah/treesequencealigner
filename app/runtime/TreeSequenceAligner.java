package runtime;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import factory.FASTAParser;

public class TreeSequenceAligner {

	public static void main(String[] args) throws IOException {
		try {
			System.out.println(Presenter.getJVMStats());
            FASTAParser parser = new FASTAParser(new File(new File("").
                    getAbsolutePath()+"/test_data/test_seqs.fasta"));
			parser.parse();
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
		} catch (IOException e) {
			System.out.println(e.getMessage());
		} catch (NullPointerException e) {
			System.out.println(e.getMessage());
		}
	}
}
