/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */
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