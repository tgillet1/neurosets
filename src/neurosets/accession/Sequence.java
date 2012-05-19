package neurosets.accession;

/**
 * Represents sequence which is used to model motif counts and meta-data.
 * Each sequence has an accession name, a set of motif counts, and a list of
 * meta-data entries mapping to this specific object.
 * 
 * @author Parsa Hosseini
 * */
public class Sequence {
	private String name; // accession name
	
	/**
	 * Creates a Sequence object given the objects accession.
	 * 
	 * @param name The Sequence accession name
	 * */
	public Sequence(String name) {
		this.setName(name);
	}

	public String getName() {
		return name;
	}

	private void setName(String name) {
		this.name = name;
	}

}
