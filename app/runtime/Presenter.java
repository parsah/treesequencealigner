package runtime;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * A helper-class solely to deal with presentation, i.e. formatting,
 * displaying messages and text.
 * */
public abstract class Presenter {
	
	/**
	 * Show text besides some formatted text, i.e. date or time
	 * */
	public static void out(String text, boolean newline) {
		Date d = new Date(); // get current date
		String time = new SimpleDateFormat("E MMM dd HH:mm:ss aaa").format(d);
		if (newline == true) { // move to next line
			System.out.println("[" + time + "]  " + text);	
		}
		else { // do not buffer to next line
			System.out.print("[" + time + "]  " + text + " ... ");	
		}
	}
	
	/**
	 * Return the JVM version and corresponding architecture
	 * @return JVM version and architecture string
	 * */
	public static String getJVMStats() {
		return "Java (" + System.getProperty("sun.arch.data.model") + 
				"bit) " + System.getProperty("java.version");
	}

}
