package neurosets.parse;

import java.io.FileFilter;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * Contains functions to parse CSV files containing neurite information.
 * To fully execute dataset parsing, two CSV files are required: one must have
 * metadata while the other must have neurite branch-lengths given a motif
 * 
 * @author Parsa Hosseini
 * 
 * */
public class CSVParser {
	
	private FileNameExtensionFilter createFilter() {
		// apply an implicit file-filter
		FileNameExtensionFilter filter = new FileNameExtensionFilter(
				"CSV file (*.csv)", "csv");
		return filter;
	}
	
	// Loads a CSV file using a file-filter for screen anything not .csv
	public void loadFile(String title) {
		JFileChooser chooser = new JFileChooser("./demo");
		chooser.setDialogTitle(title);
		chooser.setFileFilter(this.createFilter());
		
		// clause if the user selects a file
		if (chooser.showOpenDialog(null) == 0) {
			
		}
	}
	
	/**
	 * Parse a CSV metadata file such that the following conditions are met:
	 * 1) The first row are column names such that column number 2 ... N 
	 * represent motifs. 2) Each row on column 1 must represent a unique ID
	 * for a sequence. Please see demo motifs .csv file for more details.
	 * */
	public void parseMetadata() {
		this.loadFile("Select CSV file containing motifs");
	}

}
