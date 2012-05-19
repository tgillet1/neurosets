package neurosets.parse;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;
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
	private File fileMotifCounts; // CSV file containing motif counts
	private File fileMetadata; // Meta-data for each accession
	
	// Apply an implicit file-filter so that only CSV files can be selected
	private FileFilter createFilter() {
		FileFilter filter = new FileNameExtensionFilter("CSV file", "csv");
		return filter;
	}
	
	// Loads a CSV file using a file-filter for screen anything not .csv
	private File selectFile(String title) {
		File csvFile = null; // the selected CSV file
		JFileChooser chooser = new JFileChooser("./demo");
		chooser.setDialogTitle(title);
		chooser.setFileFilter(this.createFilter());
		chooser.setAcceptAllFileFilterUsed(false); // choose only CSV files
		
		// clause if the user selects a file
		if (chooser.showOpenDialog(null) == 0) {
			csvFile = chooser.getSelectedFile();
		}
		
		// return the selected CSV file
		return csvFile;
	}
	
	/**
	 * Two CSV files are required: A meta-data file as well as a motif-count
	 * file. Both files must be provided and successfully parsed.
	 * */
	public void loadMetadata() {
		// Simple string to indicate to the user what file is being sought
		String dialogTitle = "Select CSV file containing meta-data";
		File fileMetadata = this.selectFile(dialogTitle);
		this.setFileMetadata(fileMetadata);
	}
	
	/**
	 * Load the CSV which contains motif-counts. These counts represent
	 * branch-lengths across a tree given a specific motif. Such counts pertain
	 * to specific accessions. See motif CSV for detailed data organization.
	 * */
	public void loadMotifCounts() {
		// Select only the CSV containing motif counts
		String dialogTitle = "Select motif-counts file";
		File fileMotifCounts = this.selectFile(dialogTitle);
		this.setFileMotifCounts(fileMotifCounts);
	}
	
	/**
	 * Determines whether the selected CSV files are valid or not; useful
	 * since both motif counts and meta-data must be provided for any project.
	 * 
	 * @return boolean representing validity of the selected CSV files
	 * */
	public boolean isSanitized() {
		if ((this.getFileMetadata() != null) && 
				(this.getFileMotifCounts() != null)) {
			return true;
		}
		else {
			return false;
		}
	}

	public File getFileMetadata() {
		return fileMetadata;
	}

	private void setFileMetadata(File fileMetadata) {
		this.fileMetadata = fileMetadata;
	}

	public File getFileMotifCounts() {
		return fileMotifCounts;
	}

	private void setFileMotifCounts(File fileMotifCounts) {
		this.fileMotifCounts = fileMotifCounts;
	}

}
