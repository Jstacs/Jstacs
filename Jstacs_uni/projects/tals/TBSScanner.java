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
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.HomogeneousMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.galaxy.MultilineSimpleParameter;

/**
 * Class that scans input sequences for putative TAL effector target sites given a {@link TALgetterDiffSM}.
 * 
 * @author Jan Grau
 *
 */
public class TBSScanner {

	/**
	 * Granularity of p-values.
	 * @author Jan Grau
	 *
	 */
	public enum PVals{
		/**
		 * Do not compute p-values
		 */
		NONE,
		/**
		 * Determine coarse p-values
		 */
		COARSE,
		/**
		 * Determine fine p-values
		 */
		FINE
	}
	
	/**
	 * {@link CategoricalResult} that links gene names to other resources.
	 * @author Jan Grau
	 *
	 */
	public static class GeneLinkResult extends CategoricalResult{
		
		private static Object[][] map = {{Pattern.compile( "Os[0-9]+[a-zA-Z]+[0-9]+(\\.[0-9]+)?"),"http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_$$$"},
		                                 {Pattern.compile( "AT[0-9]+[a-zA-Z]+[0-9]+(\\.[0-9]+)"),"http://www.arabidopsis.org/servlets/TairObject?name=$$$&type=gene"},
		                                 {Pattern.compile( "AT[0-9]+[a-zA-Z]+[0-9]+"),"http://www.arabidopsis.org/servlets/TairObject?name=$$$&type=locus"}
		                                 };

		
		
		/**
		 * Create a new {@link GeneLinkResult}
		 * @param name the name
		 * @param comment the comment
		 * @param result the result
		 */
		public GeneLinkResult( String name, String comment, String result ) {
			super( name, comment, result );
		}

		/**
		 * Creates a new {@link GeneLinkResult} from XML.
		 * @param representation the XML representation
		 * @throws NonParsableException if <code>representation</code> could not be parsed
		 */
		public GeneLinkResult( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		@Override
		public String getValue() {
			
			String str = super.getValue().toString();
			for(int i=0;i<map.length;i++){
				Matcher m = ((Pattern)map[i][0]).matcher( str );
				if(m.find()){
					int start = m.start();
					int end = m.end();
					String url = ((String)map[i][1]).replaceAll( "\\$\\$\\$", str.substring( start, end ) );
					return str.substring( 0, start )+"<a href=\""+url+"\" target=\"_blank\">"+str.substring( start, end )+"</a>"+str.substring( end );
				}
			}
			return str;
		}
		
		
		
	}
	
	private static class BackgroundDistribution{
		
		private static final int r = 10;
		
		double[] bestScores;
		double[] restScores;
		double pivot;
		double pivotP;
		
		private BackgroundDistribution(DataSet template, AbstractDifferentiableStatisticalModel model, Sequence tal, int n, PVals pvals) throws Exception{
			int m = Math.min( 1000, n );
			
			int myR = r;
			if(pvals == PVals.NONE){
				throw new Exception();
			}else if(pvals == PVals.COARSE){
				myR = 1;
			}
			
			this.pivot = 0;
			this.restScores = new double[n];
			
			int off = 0;
			bestScores = new double[m*myR];
			for(int i=0;i<myR;i++){
				double[] scores = getSortedScores( template, model, tal, n );
				System.arraycopy( scores, scores.length-m, bestScores, off, m );
				off += m;
				
				pivot += scores[scores.length-m];
				for(int j=0;j<scores.length;j++){
					restScores[j] += scores[j];
				}
			}
			Arrays.sort( bestScores );
			pivot /= myR;
			pivotP = (double)m/(double)n;
			for(int i=0;i<restScores.length;i++){
				restScores[i] /= myR;
			}
		}
		
		private double[] getSortedScores(DataSet template, AbstractDifferentiableStatisticalModel model, Sequence tal, int n) throws Exception{

			HomogeneousMM mm = new HomogeneousMM( new HomMMParameterSet( template.getAlphabetContainer(), 4, "", (byte)2 ) );

			mm.train( template );

			DataSet gen = mm.emitDataSet( 1, n+tal.getLength() );

			Sequence gs = gen.getElementAt( 0 );
			double[] scores = new double[n];
			for(int j=0;j<n;j++){
				Sequence sub = gs.getSubSequence( j, tal.getLength()+1 );
				sub = sub.annotate( true, new ReferenceSequenceAnnotation( "seq", tal ) );
				scores[j] = model.getLogScoreFor( sub );
			}
			Arrays.sort( scores );
			return scores;
		}
		
		private double getPValue(double score){
			if(score <= pivot){
				int idx = Arrays.binarySearch( restScores, score );
				if(idx < 0){
					idx = -idx -1;
				}
				idx = restScores.length-idx;
				return (double)idx/(double)restScores.length;
			}else{
				int idx = Arrays.binarySearch( bestScores, score );
				if(idx < 0){
					idx = -idx -1;
				}
				idx = bestScores.length-idx;
				double locP = (double)idx/(double)bestScores.length;
				return pivotP*locP;
			}
		}
		
	}
	
	/**
	 * Class for the parameters of the {@link TBSScanner}.
	 * @author Jan Grau
	 *
	 */
	public static class TBSScannerParameterSet extends ParameterSet{

		/**
		 * Creates a new {@link TBSScannerParameterSet}.
		 * @throws Exception if creation fails
		 */
		public TBSScannerParameterSet() throws Exception {
			super();
			this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{"Use a previously uploaded file", "Paste sequences in FastA format"}, new ParameterSet[]{
			                                                                                                                                         new SimpleParameterSet( new FileParameter( "FastA file", "The sequences to scan for TAL effector target sites, FastA format", "fasta", true ) ),
			                                                                                                                                         new SimpleParameterSet( new MultilineSimpleParameter( "FastA sequences", "The sequences to scan for TAL effector target sites, FastA format", true ) )
			}, "Input sequences", "You can either use a previously uploaded file (see task &quot;GetData&quot; -&gt; &quot;Upload File&quot;) or paste in sequences in FastA format", true ) );
			
			this.parameters.add( new MultilineSimpleParameter( "RVD sequence", "Sequence of RVDs, seperated by '-'", true, "NI-HD-HD-NG-NN-NK-NK" ) );
			
			this.parameters.add( new SimpleParameter( DataType.INT, "Upstream offset", "Number of positions ignored at 5' end of each sequence", true, 0 ) );
			this.parameters.add( new SimpleParameter( DataType.INT, "Downstream offset", "Number of positions ignored at 3' end of each sequence", true, 0) );
			
			this.parameters.add( new SimpleParameter( DataType.INT, "Maximum number of target sites", "Limits the total number of reported target sites in all input sequences, ranked by their score.", true, new NumberValidator<Comparable<? extends Number>>( 1, 10000 ), 100 ) );
			
			SimpleParameter pValThresh = new SimpleParameter( DataType.DOUBLE, "p-Value", "Filter the reported hits by a maximum p-Value. A value of 0 or 1 switches off the filter.", true, new NumberValidator<Comparable<? extends Number>>( 0.0, 1.0 ), 1.0E-6 );
			SelectionParameter sp = new SelectionParameter( DataType.PARAMETERSET, new String[]{"No p-Values (fastest)","Coarse p-Values (faster but less accurate)","Fine-grained p-Values (slower but more accurate)"}, new ParameterSet[]{new SimpleParameterSet(  ),new SimpleParameterSet( pValThresh ), new SimpleParameterSet( pValThresh )}, "Computation of p-Values",
					"Mode to compute p-values for predicted target sites. If no p-values are computed, filtering for p-values is not available.", true);
			sp.setDefault( "Fine-grained p-Values (slower but more accurate)" );
			this.parameters.add( sp );
			
			
		}
		
		/**
		 * Returns all parameters in the set.
		 * @return the parameters
		 */
		public Parameter[] getAllParameters(){
			return this.parameters.toArray( new Parameter[0] );
		}
		
		/**
		 * Adds a parameter to the set.
		 * @param i the index where the parameter is added
		 * @param par the parameter
		 */
		public void addParameter(int i, Parameter par){
			this.parameters.add( i, par );
		}

		/**
		 * Returns the selected value of {@link PVals}.
		 * @return the selected value
		 */
		public PVals computePValues(){
			int sel = ((SelectionParameter)parameters.get( "Computation of p-Values" )).getSelected();
			if(sel == 0){
				return PVals.NONE;
			}else if(sel == 1){
				return PVals.COARSE;
			}else if(sel == 2){
				return PVals.FINE;
			}else{
				throw new RuntimeException("Computation of p-Values is required parameter.");
			}
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
		public DataSet getInputSequences() throws Exception {
			if( ((SelectionParameter)this.parameters.get( "Input sequences" )).getSelected() == 0){
				String filename = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				return new DataSet( new AlphabetContainer(new DiscreteAlphabet( true, "A", "C", "G", "T", "N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V" )), new SparseStringExtractor( filename, '>', new SimpleSequenceAnnotationParser() ) );
			}else{
				String content = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				//System.out.println(content);
				return new DataSet( new AlphabetContainer(new DiscreteAlphabet( true, "A", "C", "G", "T", "N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V" )), new SparseStringExtractor( new BufferedReader(new StringReader( content )), '>', null, new SimpleSequenceAnnotationParser() ) );
			}
		}

		/**
		 * Returns the maximum number of target sites.
		 * @return the maximum number
		 */
		public int getN(){
			return (Integer) this.parameters.get( "Maximum number of target sites" ).getValue();
		}
		
		/**
		 * Returns the threshold on the p-values.
		 * @return the threshold
		 */
		public double getPValue(){
			ParameterSet ps = (ParameterSet)parameters.get( "Computation of p-Values" ).getValue();
			if(ps.getNumberOfParameters() == 0){
				return 0;
			}else{
				return (Double)ps.getParameterAt( 0 ).getValue();
			}
		}

		/**
		 * Returns the upstream offset.
		 * @return the offset
		 */
		private int getFirstPosition(){
			return (Integer) this.parameters.get( "Upstream offset" ).getValue();
		}

		/**
		 * Returns the downstream offset.
		 * @return the offset
		 */
		private int getDownstreamOffset(){
			return (Integer) this.parameters.get( "Downstream offset" ).getValue();
		}

	}

	/**
	 * Class for a sorted list of results.
	 * @author Jan Grau
	 *
	 */
	public static class ResultList{

		private ComparableElement<ResultSet, Double>[] list;
		private int curr;

		/**
		 * Creates a new {@link ResultList} with at most <code>n</code> elements.
		 * @param n the maximum number of elements
		 */
		public ResultList(int n){
			this.list = new ComparableElement[n];
			this.curr = 0;
		}

		/**
		 * Returns if <code>values</code> is greater than the currently worst score
		 * @param value the value to test
		 * @return the result
		 */
		public boolean better(double value){
			if(list[curr] == null){
				return true;
			}else{
				if(-value < list[curr].getWeight()){
					return true;
				}
			}
			return false;
		}

		/**
		 * Adds the result <code>res</code> with score <code>value</code> to this set.
		 * @param res the result
		 * @param value the value
		 */
		public void add(ResultSet res, double value){
			if(list[curr] == null){
				list[curr] = new ComparableElement<ResultSet, Double>( res, -value );
				if(curr < list.length-1){
					curr++;
					Arrays.sort( list, 0, curr );
				}else{
					Arrays.sort( list );
				}
				
			}else{
				if(-value < list[curr].getWeight()){
					list[curr] = new ComparableElement<ResultSet, Double>( res, -value );
					Arrays.sort( list );
				}
			}
		}

		/**
		 * Returns the predicted binding sites of this {@link ResultList}.
		 * @return the sites
		 * @throws Exception if the data set of sites could not be created
		 */
		public DataSet getBindingSites() throws Exception {
			if(getNumberOfResults() == 0){
				return null;
			}
			Sequence[] seqs = new Sequence[getNumberOfResults()];
			for(int i=0;i<seqs.length;i++){
				seqs[i] = Sequence.create( DNAAlphabetContainer.SINGLETON, new SequenceAnnotation[]{
				                                                                                    new SequenceAnnotation( "ID", (String)list[i].getElement().getResultForName( "ID" ).getValue() ),
				                                                                                    new SequenceAnnotation( "Position", list[i].getElement().getResultForName( "Position" ).getValue().toString() ),
				                                                                                    new SequenceAnnotation( "Score", list[i].getElement().getResultForName( "Score" ).getValue().toString() )
				}, (String)list[i].getElement().getResultForName( "Sequence" ).getValue(), "" );
			}
			return new DataSet( "binding sites", seqs );
		}

		/**
		 * Returns an array containing all results.
		 * @return the array
		 */
		public ResultSet[] toArray(){
			ResultSet[] res = new ResultSet[getNumberOfResults()];
			for(int i=0;i<res.length;i++){
				res[i] = list[i].getElement();
			}
			return res;
		}

		/**
		 * Returns the score of the best element.
		 * @return the score
		 */
		public double getBestScore(){
			if(list[0] == null){
				return Double.NEGATIVE_INFINITY;
			}
			return -list[0].getWeight();
		}

		/**
		 * Returns the score of the worst element.
		 * @return the score
		 */
		public double getWorstScore(){
			if(curr == 0 && list[curr] == null){
				return Double.NEGATIVE_INFINITY;
			}
			if(list[curr] == null){
				return -list[curr-1].getWeight();
			}else{
				return -list[curr].getWeight();
			}
			
		}

		/**
		 * Returns the number of results.
		 * @return the number
		 */
		public int getNumberOfResults(){
			if(list[curr] == null){
				return curr;
			}else{
				return curr+1;
			}
		}

	}
	
	/**
	 * Scans the input sequences given in <code>ds</code> for TAL effector target sites given parameters <code>params</code> and a {@link TALgetterDiffSM} <code>model</code>.
	 * @param model the model
	 * @param params the parameters
	 * @param ds the input data
	 * @return the predicted target sites
	 * @throws Exception if something went wrong
	 */
	public static ResultList[] scan(TALgetterDiffSM model, TBSScannerParameterSet params, DataSet ds) throws Exception {
		if( !params.hasDefaultOrIsSet() ) {
			System.err.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}


		Pair<int[][],DataSet> pair = preprocess( ds );

		int[][] offsets = pair.getFirstElement();
		ds = pair.getSecondElement();

		//RVD-Alph
		String[] alph={"NI","NG","NN","NS","N*","ND","NK","NC","NV","NA","NH","HD","HG","HA","H*","HH","HI","HN","S*","SN","SS","IG","YG","NP","NT","IS"};
		//System.out.println("alph: "+alph.length);
		AlphabetContainer alphabetsRVD= new AlphabetContainer( new DiscreteAlphabet(true,alph) );

		Sequence tal = Sequence.create( alphabetsRVD, params.getTALSequence(), "-" );

		double size = getSize( ds, tal );
		if(size > 2E8){
			System.err.println("Data set too large. Currently at most 90 Mb allowed");
		}
		
		
		
		int firstPos = params.getFirstPosition();
		int downstreamOff = params.getDownstreamOffset();

		ResultList rl = new ResultList( params.getN() );
		ResultList rl2 = new ResultList( params.getN() );
		
		PVals pvals = params.computePValues();
		
		BackgroundDistribution bgd = null;
		if(pvals != PVals.NONE){
			bgd = new BackgroundDistribution( ds, model, tal, Math.max( (int)1E7, (int)size ), pvals );
		}
		
		double pValThresh = params.getPValue();
		
		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence seq = ds.getElementAt( i );

			String id = (String) seq.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue();

			id = id.trim();

			for(int j=0;j<seq.getLength()-tal.getLength();j++){
				Sequence sub = seq.getSubSequence( j, tal.getLength()+1 );
				sub = sub.annotate( true, seq.getAnnotation() );
				sub = sub.annotate( true, new ReferenceSequenceAnnotation( "seq", tal ) );
				double score = model.getLogScoreFor( sub );

				if(rl.better( score )){

					int off = j;

					double pVal = 0;
					if(bgd != null){
						pVal = bgd.getPValue( score );
					}
					if(pValThresh == 0 || pValThresh >= pVal){
						double eVal = pVal * size;

						int dist = (seq.getLength()-(j+tal.getLength()+1))+offsets[1][i];

						if(j >= firstPos && dist >= downstreamOff){

							String str = model.getMatchString( sub );

							Result[] tr = new Result[]{
							                           new CategoricalResult( "ID", "", id ),
							                           new NumericalResult( "Position", "", j+offsets[0][i] ),
							                           new NumericalResult( "Distance to end", "", (seq.getLength()-(j+tal.getLength()+1))+offsets[1][i] ),
							                           new CategoricalResult( "Sequence", "", sub.toString() ),
							                           new CategoricalResult( "Matches", "", str ),//TODO
							                           new NumericalResult( "Score", "", score ),
							                           pvals == PVals.NONE ? new CategoricalResult( "p-value", "", "NA" ) : new NumericalResult( "p-value", "", pVal ),
							                           pvals == PVals.NONE ? new CategoricalResult( "E-value", "", "NA" ) : new NumericalResult( "E-value", "", eVal)

							};					

							rl.add( new ResultSet( tr ),score );
							tr[0] = new GeneLinkResult( "ID", "", id );
							rl2.add( new ResultSet( tr ), score );
						}
					}
				}

			}

		}

		return new ResultList[]{rl,rl2};
	}


	public static Pair<int[][], DataSet> preprocess(DataSet data) throws Exception{

		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );

		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		IntList starts = new IntList();
		IntList ends = new IntList();
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt( i );
			String seqStr = seq.toString();
			Matcher match = acgt.matcher( seqStr );
			while(match.find()){
				int start = match.start();
				int end = match.end();
				seqs.add( Sequence.create( DNAAlphabetContainer.SINGLETON, seq.getAnnotation(), seqStr.substring( start, end ), "" ) );
				starts.add( start );
				ends.add( seq.getLength()-end );
			}
		}
		return new Pair<int[][],DataSet>(new int[][]{starts.toArray(),ends.toArray()},new DataSet("",seqs.toArray( new Sequence[0] )));
	}

	private static double getSize(DataSet data, Sequence tal){
		double size = 0;
		for(int i=0;i<data.getNumberOfElements();i++){
			size += data.getElementAt( i ).getLength()-tal.getLength();
		}
		return size;
	}


	private static HashSet<Sequence> makeHash(DataSet d){
		HashSet<Sequence> set = new HashSet<Sequence>();
		for(int i=0;i<d.getNumberOfElements();i++){
			set.add( d.getElementAt( i ) );
		}
		return set;
	}
	
}
