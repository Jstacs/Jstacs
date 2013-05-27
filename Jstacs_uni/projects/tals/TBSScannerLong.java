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
package projects.tals;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import projects.tals.TBSScanner.ResultList;
import projects.tals.TBSScanner.TBSScannerParameterSet;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.galaxy.MultilineSimpleParameter;

/**
 * Class that scans input sequences for putative TAL effector target sites given a {@link TALgetterDiffSM},
 * adapted for large input data.
 * 
 * @author Jan Grau
 * @see TBSScanner
 *
 */
public class TBSScannerLong {

	
	/**
	 * Class for the parameters of the {@link TBSScannerLong}.
	 * @author Jan Grau
	 *
	 */
	public static class TBSScannerLongParameterSet extends ParameterSet{

		/**
		 * Creates a new {@link TBSScannerLongParameterSet}.
		 * @throws Exception if creation fails
		 */
		public TBSScannerLongParameterSet() throws Exception {
			super();
			this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{"Use a previously uploaded file", "Paste sequences in FastA format"}, new ParameterSet[]{
			                                                                                                                                         new SimpleParameterSet( new FileParameter( "FastA file", "The sequences to scan for TAL effector target sites, FastA format", "fasta", true ) ),
			                                                                                                                                         new SimpleParameterSet( new MultilineSimpleParameter( "FastA sequences", "The sequences to scan for TAL effector target sites, FastA format", true ) )
			}, "Input sequences", "You can either use a previously uploaded file (see task &quot;GetData&quot; -&gt; &quot;Upload File&quot;) or paste in sequences in FastA format", true ) );
			
			this.parameters.add( new MultilineSimpleParameter( "RVD sequence", "Sequence of RVDs, seperated by '-'", true, "NI-HD-HD-NG-NN-NK-NK" ) );
					
			this.parameters.add( new SimpleParameter( DataType.INT, "Maximum number of target sites", "Limits the total number of reported target sites in all input sequences, ranked by their score.", true, new NumberValidator<Comparable<? extends Number>>( 1, 10000 ), 100 ) );
			
		}
		
		/**
		 * Returns all parameters in the set.
		 * @return the parameters
		 */
		public Parameter[] getAllParameters(){
			return this.parameters.toArray( new Parameter[0] );
		}
		
		/**
		 * Returns all parameters in the set.
		 * @return the parameters
		 */
		public void addParameter(int i, Parameter par){
			this.parameters.add( i, par );
		}

		/**
		 * Returns the sequence of RVDs of the TAL effector.
		 * @return the RVD sequence
		 */
		public String getTALSequence(){
			return (String)this.parameters.get( "RVD sequence" ).getValue();
		}

		/**
		 * Sets the path to the input files.
		 * @param path the path
		 * @throws Exception if the path could not be set
		 */
		public void setInputPath(String path) throws Exception {
			((SelectionParameter)this.parameters.get( "Input sequences" )).setValue( "Use a previously uploaded file" );
			((FileParameter)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 )).setValue( new FileParameter.FileRepresentation( path ) );
		}
		
		/**
		 * Returns the input sequences.
		 * @return the input sequences
		 * @throws Exception if the input sequences could not be read
		 */
		public BufferedReader getInputSequences() throws Exception {
			if( ((SelectionParameter)this.parameters.get( "Input sequences" )).getSelected() == 0){
				String filename = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				return new BufferedReader( new FileReader( filename ) );
			}else{
				String content = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				return new BufferedReader(new StringReader( content ));
			}
		}

		/**
		 * Returns the maximum number of target sites.
		 * @return the maximum number
		 */
		public int getN(){
			return (Integer) this.parameters.get( "Maximum number of target sites" ).getValue();
		}
		
	}
	
	/**
	 * Scans the input sequences given in <code>dsRead</code> for TAL effector target sites given parameters <code>params</code> and a {@link TALgetterDiffSM} <code>model</code>.
	 * @param model the model
	 * @param params the parameters
	 * @param dsRead the input data
	 * @return the predicted target sites
	 * @throws Exception if something went wrong
	 */
	public static ResultList[] scan(TALgetterDiffSM model, TBSScannerLongParameterSet params, BufferedReader dsRead) throws Exception {
		if( !params.hasDefaultOrIsSet() ) {
			System.err.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}



		//RVD-Alph
		String[] alph={"NI","NG","NN","NS","N*","ND","NK","NC","NV","NA","NH","HD","HG","HA","H*","HH","HI","HN","S*","SN","SS","IG","YG","NP","NT","IS"};
		//System.out.println("alph: "+alph.length);
		AlphabetContainer alphabetsRVD= new AlphabetContainer( new DiscreteAlphabet(true,alph) );

		Sequence tal = Sequence.create( alphabetsRVD, params.getTALSequence(), "-" );

		

		ResultList rl = new ResultList( params.getN() );
		ResultList rl2 = new ResultList( params.getN() );
	
		StringBuffer lastHeader = new StringBuffer();
		
		Pair<int[][],Sequence[]> curr = null;
		
		while( (curr = readNextSequences( dsRead, lastHeader )) != null){

			int[][] offsets = curr.getFirstElement();
			Sequence[] seqs = curr.getSecondElement();

			for(int s=0;s<seqs.length;s++){
				Sequence seq = seqs[s];

				String id = (String) seq.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue();

				id = id.trim();

				for(int j=0;j<seq.getLength()-tal.getLength();j++){
					Sequence sub = seq.getSubSequence( j, tal.getLength()+1 );
					sub = sub.annotate( true, seq.getAnnotation() );
					sub = sub.annotate( true, new ReferenceSequenceAnnotation( "seq", tal ) );
					double score = model.getLogScoreFor( sub );

					if(rl.better( score )){




						String str = model.getMatchString( sub );

						Result[] tr = new Result[]{
						                           new CategoricalResult( "ID", "", id ),
						                           new NumericalResult( "Position", "", j+offsets[0][s] ),
						                           new NumericalResult( "Distance to end", "", (seq.getLength()-(j+tal.getLength()+1))+offsets[1][s] ),
						                           new CategoricalResult( "Sequence", "", sub.toString() ),
						                           new CategoricalResult( "Matches", "", str ),//TODO
						                           new NumericalResult( "Score", "", score )

						};					

						rl.add( new ResultSet( tr ),score );
						tr[0] = new TBSScanner.GeneLinkResult( "ID", "", id );
						rl2.add( new ResultSet( tr ), score );

					}
				}
			}
		}
		
		return new ResultList[]{rl,rl2};
	}

	
	private static Pair<int[][],Sequence[]> readNextSequences(BufferedReader read, StringBuffer lastHeader) throws Exception {
		
		String str = null;
		
		StringBuffer line = new StringBuffer();
		
		IntList starts = new IntList();
		IntList ends = new IntList();
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		
		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		while( (str = read.readLine()) != null || line.length() > 0 ){
			if(str != null){
				str = str.trim();
			}
			if(str == null || str.startsWith( ">" )){//next sequence
				String header = lastHeader.toString();
				if(str != null){
					lastHeader.delete( 0, lastHeader.length() );
					lastHeader.append( str.substring( 1 ) );
				}
				if(line.length() > 0){//we have a sequence
					String seqStr = line.toString();
					line.delete( 0, line.length() );
					Matcher match = acgt.matcher( seqStr );
					while(match.find()){
						int start = match.start();
						int end = match.end();
						SequenceAnnotation annotation = new SequenceAnnotation( "unparsed comment line", "unparsed comment line", new CategoricalResult( "unparsed comment", "unparsed comment", header ) );
						Sequence seq = new SparseSequence( DNAAlphabetContainer.SINGLETON, seqStr.substring( start, end ) );
						seq = seq.annotate( false, annotation );
						seqs.add( seq );
						
						starts.add( start );
						ends.add( seqStr.length()-end );
					}
					return new Pair<int[][],Sequence[]>(new int[][]{starts.toArray(),ends.toArray()},seqs.toArray( new Sequence[0] ));
				}
			}else{
				line.append( str );
			}	
		}
		return null;
	}

	
	private static HashSet<Sequence> makeHash(DataSet d){
		HashSet<Sequence> set = new HashSet<Sequence>();
		for(int i=0;i<d.getNumberOfElements();i++){
			set.add( d.getElementAt( i ) );
		}
		return set;
	}
	
}
