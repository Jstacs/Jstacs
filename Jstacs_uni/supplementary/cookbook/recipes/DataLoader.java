package supplementary.cookbook.recipes;

import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.bioJava.BioJavaAdapter;
import de.jstacs.data.bioJava.SimpleSequenceIterator;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.StringExtractor;

/**
 * This class exemplarily shows how to load data.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DataLoader {

	/**
	 * @param args args[0] contains the home path, i.e., the path to myfile.fa
	 */
	public static void main(String[] args) throws Exception {
		String home = args[0];
		
		//load DNA sequences in FastA-format
		DataSet data = new DNADataSet( home+"myfile.fa" ); 
		
		//create a DNA-alphabet
		AlphabetContainer container = DNAAlphabetContainer.SINGLETON;
		
		//create a DataSet using the alphabet from above in FastA-format
		data = new DataSet( container, new SparseStringExtractor( home+"myfile.fa", StringExtractor.FASTA ));
		
		//create a DataSet using the alphabet from above
		data = new DataSet( container, new SparseStringExtractor( home+"myfile.txt" ));
		
		//defining the ids, we want to obtain from NCBI Genbank:
		GenbankRichSequenceDB db = new GenbankRichSequenceDB();
		
		SimpleSequenceIterator it = new SimpleSequenceIterator(
				db.getRichSequence( "NC_001284.2" ),
				db.getRichSequence( "NC_000932.1" )
				);
		 
		//conversion to Jstacs DataSet
		data = BioJavaAdapter.sequenceIteratorToDataSet( it, null );
	}
}