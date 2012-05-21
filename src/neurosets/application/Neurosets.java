package neurosets.application;

import java.io.FileNotFoundException;
import java.io.IOException;

import neurosets.parse.CSVParser;

// Represents the main application class; drives analysis
class Neurosets {

	// Initiate analysis
	public static void main(String args[]) {
		try {
			CSVParser csvParser = new CSVParser();
			csvParser.loadData();
			if (csvParser.isSanitized()) {
				csvParser.parseMotifCounts(csvParser.getFileMotifCounts());
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
		} catch (IOException e) {
			System.out.println(e.getMessage());
		}
	}

}
