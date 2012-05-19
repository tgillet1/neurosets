package neurosets.application;

import neurosets.parse.CSVParser;

// Represents the main application class; drives analysis
class Neurosets {
	
	// Initiate analysis
	public static void main(String args[]) {
		CSVParser csvParser = new CSVParser();
		csvParser.loadMetadata();
		csvParser.loadMotifCounts();
		System.out.println(csvParser.isSanitized());
	}

}
