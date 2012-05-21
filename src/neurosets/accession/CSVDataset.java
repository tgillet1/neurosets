package neurosets.accession;

import java.util.ArrayList;
import java.util.List;

public class CSVDataset {
	private int numCols;
	private int numRows;
	private List<String> headers; // first line of the CSV file
	
	/**
	 * CSV File objects which are parsed are encapsulated in a CSVDataset
	 * object for easy manipulation. Such an object contains behaviors and
	 * states to delve into its row and column values.
	 * */
	public CSVDataset() {
		this.setHeaders(new ArrayList<String>());
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

	/**
	 * @return the headers
	 */
	public List<String> getHeaders() {
		return headers;
	}

	/**
	 * @param headers the headers to set
	 */
	private void setHeaders(List<String> headers) {
		this.headers = headers;
	}

}
