package model;

/**
 * A Matrix is a data-structure with an arbitrary number of rows and
 * columns. Such an object also has row-names and respective column 
 * names. These two attributes are optional but are built-upon especially
 * with more specific sub-class functionalities.
 */

public class Matrix {
	
	private double[][] matrix; // encapsulates actual dataset.
	private String columnNames; // column names
	private String rowNames; // row names
	
	/**
	 * Create a bare-bones Matrix object.
	 * */
	public Matrix() {
		this.setColumnNames(null);
		this.setMatrix(null);
		this.setRowNames(null);
	}
	
	/**
	 * Create a Matrix object given only a 2D-array of values.
	 * */
	public Matrix(double[][] data) {
		this.setMatrix(data);
	}
	
	/**
	 * Helper-method which interates through a given Matrix and prints-out
	 * all the values in each cell.
	 * */
	public void debug() {
		System.out.println("--- MATRIX DEBUGGER ---");
		for (int i = 0; i < this.getHeight(); i++) {
			for (int j = 0; j < this.getWidth(); j++) {
				System.out.print(this.getMatrix()[i][j] + "\t");
			}
			System.out.println();
		}
		System.out.println("--- END MATRIX DEBUGGER ---");
	}
	
	/**
	 * Returns the width of a 2D matrix
	 * @return width of Matrix (number of columns)
	 * */
	public int getWidth() {
		return this.getMatrix()[0].length;
	}
	
	/**
	 * Returns the height of a 2D matrix
	 * @return height of Matrix (number of rows)
	 * */
	public int getHeight() {
		return this.getMatrix().length;
	}
	
	/**
	 * Given a specific row index, return an array of values from that
	 * respective Matrix.
	 * @return list of values referencing a given row.
	 * */
	public double[] getRow(int row) {
		return this.getMatrix()[row];
	}
	
	/**
	 * Get values for a specific column, given the desired column number.
	 * @return list of values referencing a specific column.
	 * */
	public double[] getColumn(int column) {
		double[] col = new double[this.getHeight()];
		for (int i = 0; i < this.getHeight(); i++) {
			double value = this.getRow(i)[column]; // get value per column
			col[i] = value; // set the column value
		}
		return col;
	}

	/**
	 * @return the matrix
	 */
	public double[][] getMatrix() {
		return matrix;
	}

	/**
	 * @param matrix the matrix to set
	 */
	public void setMatrix(double[][] matrix) {
		this.matrix = matrix;
	}

	/**
	 * @return the columnNames
	 */
	public String getColumnNames() {
		return columnNames;
	}

	/**
	 * @param columnNames the columnNames to set
	 */
	public void setColumnNames(String columnNames) {
		this.columnNames = columnNames;
	}

	/**
	 * @return the rowNames
	 */
	public String getRowNames() {
		return rowNames;
	}

	/**
	 * @param rowNames the rowNames to set
	 */
	public void setRowNames(String rowNames) {
		this.rowNames = rowNames;
	}

}
