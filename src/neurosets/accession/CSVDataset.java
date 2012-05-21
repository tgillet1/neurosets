package neurosets.accession;

public class CSVDataset {
	private int numCols;
	private int numRows;
	
	/**
	 * CSV File objects which are parsed are encapsulated in a CSVDataset
	 * object for easy manipulation. Such an object contains behaviors and
	 * states to delve into its row and column values.
	 * */
	public CSVDataset() {
		
	}
	
	/**
	 * @return the numCols
	 */
	public int getNumCols() {
		return numCols;
	}
	/**
	 * @param numCols the numCols to set
	 */
	public void setNumCols(int numCols) {
		this.numCols = numCols;
	}
	/**
	 * @return the numRows
	 */
	public int getNumRows() {
		return numRows;
	}
	/**
	 * @param numRows the numRows to set
	 */
	public void setNumRows(int numRows) {
		this.numRows = numRows;
	}

}
