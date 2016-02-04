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

package projects.slim;

import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;


public class SlimDimontPredictorWebParameterSet extends SlimDimontPredictorParameterSet {

	public SlimDimontPredictorWebParameterSet() throws Exception {
		
		super();
		this.parameters.remove( 0 );
		this.parameters.remove( 0 );
		this.parameters.remove( 0 );
		this.parameters.remove( 0 );
		
		this.parameters.add( 0, new FileParameter( "SlimDimont", "The trained SlimDimont classifier, i.e. the &quot;SlimDimont&quot; output of a previous SlimDimont run.", "xml", true ) );
		this.parameters.add( 1, new FileParameter( "Input sequences", "The input sequences for de-novo motif discovery (can be uploaded using &quot;GetData&quot; -&gt; &quot;Upload File&quot;), annotated FastA format. The required format is described in the help section.", "fasta", true ) );
		
	}

	public DataSet getInputSequences() throws Exception {
		String filename = (String)this.parameters.get( "Input sequences" ).getValue();
		return SparseSequence.getDataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( filename, '>', new SplitSequenceAnnotationParser(":",";") ) );
		
	}
	
}
