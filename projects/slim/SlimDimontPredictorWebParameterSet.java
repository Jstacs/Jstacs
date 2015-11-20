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
