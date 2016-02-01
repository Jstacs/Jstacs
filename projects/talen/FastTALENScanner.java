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

package projects.talen;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.naming.OperationNotSupportedException;

import projects.talen.InfixTALENTargetFinder.TALENMatch;
import projects.talen.MatchFinder.Match;
import projects.tals.TALgetterDiffSM;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.GFFParser;
import de.jstacs.data.GFFParser.GFFEntry;
import de.jstacs.data.GFFParser.GFFEntry.GFFType;
import de.jstacs.data.GFFParser.GFFList;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.tools.ui.galaxy.MultilineSimpleParameter;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;


public class FastTALENScanner {

	public static String[] filters = new String[]{"Loose (q=0.35)", "Medium-Loose (q=0.375)", "Medium (q=0.4)", "Medium-Strict (q=0.45)", "Strict (q=0.5)", "Custom"};
	
	private static String[] archs = new String[]{"Cermak et al., NAR, 2011: 15-24 bp", "Li et al., NAR, 2011: 16-31 bp", "Miller et al., NBT, 2011, +28: 12-20 bp", "Miller et al., NBT, 2011, +63: 12-24 bp", "Mussolino et al., NAR, 2011: 12-20 bp", "Christian et al., PONE, 2012, +18: 13-17 bp", "Christian et al., PONE, 2012, full length: 15-24 bp", "Custom"};
	
	public enum Output{
		NONE,
		GFF3,
		GFF2
	}
	
	
	public static class FastTALENScannerParameterSet extends ParameterSet{

		public FastTALENScannerParameterSet() throws Exception {
			super();
			this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{"Use a previously uploaded file", "Paste sequences in FastA format"}, new ParameterSet[]{
			                                                                                                                                         new SimpleParameterSet( new FileParameter( "FastA file", "The sequences to scan for TAL effector target sites, FastA format", "fasta", true ) ),
			                                                                                                                                         new SimpleParameterSet( new MultilineSimpleParameter( "FastA sequences", "The sequences to scan for TAL effector target sites, FastA format", true ) )
			}, "Input sequences", "The DNA input sequences that are scanned for TALEN off targets. You can either use a previously uploaded file (see task &quot;GetData&quot; -&gt; &quot;Upload File&quot;), files from a data library (see &quot;Data Libraries&quot; of the menu item &quot;Shared Data&quot; on top of this page) or paste in sequences in FastA format", true ) );
			
			this.parameters.add( new FileParameter( "Annotation file", "A file containing genomic annotations (e.g., genes, mRNAs, exons) in GFF, GTF, or UCSC known genes BED format", "gff,gtf,txt", false ) );
			
			this.parameters.add( new MultilineSimpleParameter( "First RVD sequence", "The sequence of RVDs of the first TALEN monomer, seperated by '-'", true, "NI-HD-HD-NG-NN-NK-NK" ) );
			
			this.parameters.add( new MultilineSimpleParameter( "Second RVD sequence", "The sequence of RVDs of the second TALEN monomer, seperated by '-'", true, "NI-HD-HD-NG-NN-NK-NK" ) );
			
			this.parameters.add( new SimpleParameter( DataType.BOOLEAN, "N-Terminal first", "For the first RVD sequence, consider the architecture, where the endonuclease domain is used to the N-terminus instead of the standard C-terminal architecture", true, false ) );
			
			this.parameters.add( new SimpleParameter( DataType.BOOLEAN, "N-Terminal second", "For the second RVD sequence, consider the architecture, where the endonuclease domain is used to the N-terminus instead of the standard C-terminal architecture", true, false ) );
			
			this.parameters.add( new SimpleParameter( DataType.BOOLEAN, "Hetero-dimers only", "Consider only hetero-dimers of TALEN monomers instead of the standard search for TALEN hetero and homo-dimers", true, false ) );
			
			SimpleParameter min = new SimpleParameter( DataType.INT, "Minimum distance", "Minimum distance between TALEN monomer target sites", true, new NumberValidator<Comparable<Integer>>( 0, 100 ), 12 );
			
			SimpleParameter max = new SimpleParameter( DataType.INT, "Maximum distance", "Maximum distance between TALEN monomer target sites", true, new NumberValidator<Comparable<Integer>>( 0, 100 ), 24 );
			
			SelectionParameter architecture = new SelectionParameter( DataType.PARAMETERSET, archs, new Object[]{
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( min, max )
			}, "Architecture", "Select a pre-defined set of minimum and maximum distance between monomer target sites, i.e. spacer lengths, or enter custom values.", true );
			architecture.setDefault( "Custom" );
			
			this.parameters.add( architecture );
			
			SelectionParameter filter = new SelectionParameter( DataType.PARAMETERSET, filters, new Object[]{
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( ),
			                                                                                                 new SimpleParameterSet( new SimpleParameter( DataType.DOUBLE, "q", "The parameter q that defines the score that must be gained by an off-target site relative to the best-matching site", true, new NumberValidator<Double>( 0.35, 1.0 ), 0.4 ) )
			}, "Filter", "Filter off-targets using different thresholds on the score relative to the best-matching site", true );
			filter.setDefault( "Medium (q=0.4)" );
			
			this.parameters.add( filter );
			
			this.parameters.add( new MultilineSimpleParameter( "RVD specificities", "Set RVD specificities for additional RVDs or override existing RVD specificities. Details in help text.", false ) );
			
			
			this.parameters.add( new SimpleParameter( DataType.INT, "Maximum number of targets", "Limits the total number of reported targets in all input sequences, ranked by their score", true, new NumberValidator<Comparable<? extends Number>>( 1, 100000 ), 100 ) );
			
			this.parameters.add( new EnumParameter( Output.class, "Create an additional output file in GFF3 or GFF2 format that can be viewed in your favourite genome browser", true ) );
		
		}
		
		public double getFilter(){
			SelectionParameter filter = (SelectionParameter)parameters.get( "Filter" );
			int sel = filter.getSelected();
			switch(sel){
				case 0: return 0.35;
				case 1: return 0.375;
				case 2: return 0.4;
				case 3: return 0.45;
				case 4: return 0.5;
				case 5:
					return (Double)((SimpleParameterSet)filter.getValue()).getParameterAt( 0 ).getValue();
				default: throw new RuntimeException( "Filter not defined" );
			}
		}
		
		
		public Parameter[] getAllParameters(){
			return this.parameters.toArray( new Parameter[0] );
		}
		
		public void addParameter(int i, Parameter par){
			this.parameters.add( i, par );
		}

		public String getLeftTALSequence(){
			return (String)this.parameters.get( "First RVD sequence" ).getValue();
		}
		
		public String getRightTALSequence(){
			return (String)this.parameters.get( "Second RVD sequence" ).getValue();
		}
		
		public boolean getNTerm1(){
			return (Boolean)this.parameters.get( "N-Terminal first" ).getValue();
		}
		
		public boolean getNTerm2(){
			return (Boolean)this.parameters.get( "N-Terminal second" ).getValue();
		}
		
		public boolean getHeteroOnly(){
			return (Boolean)this.parameters.get( "Hetero-dimers only" ).getValue();
		}
		
		public int getMinimumDistance(){
			SelectionParameter arch = (SelectionParameter)parameters.get( "Architecture" );
			int sel = arch.getSelected();
			switch(sel){
				case 0 : return 15;
				case 1: return 16;
				case 2: return 12;
				case 3: return 12;
				case 4: return 12;
				case 5: return 13;
				case 6: return 15;
				case 7:
					return (Integer)((SimpleParameterSet)arch.getValue()).getParameterAt( 0 ).getValue();
				default: throw new RuntimeException( "Architecture not defined" );
			}
		}

		public int getMaximumDistance(){
			SelectionParameter arch = (SelectionParameter)parameters.get( "Architecture" );
			int sel = arch.getSelected();
			switch(sel){
				case 0 : return 24;
				case 1: return 31;
				case 2: return 20;
				case 3: return 24;
				case 4: return 20;
				case 5: return 17;
				case 6: return 24;
				case 7:
					return (Integer)((SimpleParameterSet)arch.getValue()).getParameterAt( 1 ).getValue();
				default: throw new RuntimeException( "Architecture not defined" );
			}
		}

		
		public void setInputPath(String path) throws Exception {
			((SelectionParameter)this.parameters.get( "Input sequences" )).setValue( "Use a previously uploaded file" );
			((FileParameter)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 )).setValue( new FileParameter.FileRepresentation( path ) );
		}
		
		public BufferedReader getInputSequences() throws Exception {
			if( ((SelectionParameter)this.parameters.get( "Input sequences" )).getSelected() == 0){
				String filename = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				return new BufferedReader( new FileReader( filename ) );
			}else{
				String content = (String)((SimpleParameterSet)((SelectionParameter)this.parameters.get( "Input sequences" )).getValue()).getParameterAt( 0 ).getValue();
				return new BufferedReader(new StringReader( content ));
			}
		}
		
		public GFFList getAnnotation() throws IOException{
			if(this.parameters.get( "Annotation file" ).isSet()){
				FileParameter par = (FileParameter)this.parameters.get( "Annotation file" );
				FileRepresentation rep = par.getFileContents();
				return GFFParser.parse( rep.getFilename(), rep.getExtension() );
			}else{
				return null;
			}
		}

		public int getN(){
			return (Integer) this.parameters.get( "Maximum number of targets" ).getValue();
		}

		public void setAnnotationPath( String value ) throws Exception {
			this.parameters.get( "Annotation file" ).setValue( new FileParameter.FileRepresentation( value ) );
		}

		public Output getOutput() {
			return (Output)parameters.get( Output.class.getSimpleName() ).getValue();
		}

		public void setRVDSpecificities( String filename ) throws IOException, IllegalValueException {
			BufferedReader reader = new BufferedReader( new FileReader( filename ) );
			String line = null;
			StringBuffer temp = new StringBuffer();
			while( (line = reader.readLine()) != null ){
				temp.append( line );
				temp.append( "\n" );
			}
			this.parameters.get( "RVD specificities" ).setValue( temp.toString() );
		}
		
		public Pair<String[], double[][][]> getSpecificities() throws WrongAlphabetException{
			String str = (String)this.parameters.get( "RVD specificities" ).getValue();
			if(str == null || str.trim().length() == 0){
				return null;
			}
			String[] lines = str.split( "\n" );
			LinkedList<String> rvdList = new LinkedList<String>();
			LinkedList<double[]> specList = new LinkedList<double[]>();
			double[] fp = null;
			DiscreteAlphabet dna = DNAAlphabet.SINGLETON;
			for(int i=0;i<lines.length;i++){
				String[] parts1 = lines[i].split( ":" );
				String[] rvds = parts1[0].split( "," );
				String[] specs = parts1[1].split( "," );
				for(int j=0;j<rvds.length;j++){
					rvds[j] = rvds[j].trim();
				}
				double[] specd = new double[4];
				for(int j=0;j<specs.length;j++){
					specs[j] = specs[j].trim();
					if(dna.isSymbol( specs[j] )){
						specd[ dna.getCode( specs[j] ) ] = 1.0;
					}else{
						specd[j] = Double.parseDouble( specs[j] );
					}
				}
				Normalisation.sumNormalisation( specd );
				for(int j=0;j<rvds.length;j++){
					if("0".equals( rvds[j] )){
						fp = specd;
					}else{
						rvdList.add( rvds[j] );
						specList.add( specd );
					}
				}
			}
			return new Pair<String[], double[][][]>( rvdList.toArray(new String[0]), new double[][][]{specList.toArray( new double[0][] ),{fp}} );
		}
		
	}

	public static class ResultList{

		private ComparableElement<Object, Double>[] list;
		private int curr;

		public ResultList(int n){
			this.list = new ComparableElement[n];
			this.curr = 0;
		}

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

		public void add(Object res, double value){
			if(list[curr] == null){
				list[curr] = new ComparableElement<Object, Double>( res, -value );
				if(curr < list.length-1){
					curr++;
					Arrays.sort( list, 0, curr );
				}else{
					Arrays.sort( list );
				}
				
			}else{
				if(-value < list[curr].getWeight()){
					list[curr] = new ComparableElement<Object, Double>( res, -value );
					Arrays.sort( list );
				}
			}
		}

		public ResultSet[] resultsToArray(){
			ResultSet[] res = new ResultSet[getNumberOfResults()];
			for(int i=0;i<res.length;i++){
				if(list[i].getElement() instanceof Object[]){
					res[i] = (ResultSet)((Object[])list[i].getElement())[0];
				}else{
					res[i] = (ResultSet)list[i].getElement();
				}
			}
			return res;
		}
		
		public GFFList gffToList(){
			LinkedList<GFFEntry> ens = new LinkedList<GFFParser.GFFEntry>();
			for(int i=0;i<getNumberOfResults();i++){
				if(list[i].getElement() instanceof Object[]){
					GFFEntry[] en = (GFFEntry[])((Object[])list[i].getElement())[1];
					for(int j=0;j<en.length;j++){
							ens.add(en[j]);
					}
				}else{
					return null;
				}
			}
			return new GFFList( ens );
		}

		public double getBestScore(){
			if(list[0] == null){
				return Double.NEGATIVE_INFINITY;
			}
			return -list[0].getWeight();
		}

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

		public int getNumberOfResults(){
			if(list[curr] == null){
				return curr;
			}else{
				return curr+1;
			}
		}

	}
	
	
	public ListResult[] scan(TALgetterDiffSM model, FastTALENScannerParameterSet params, BufferedReader dsRead, GFFList gff, int numThreads ) throws Exception {
		if( !params.hasDefaultOrIsSet() ) {
			System.err.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}

		//RVD-Alph
		//String[] alph={"NI","NG","NN","NS","N*","ND","NK","NC","NV","NA","NH","HD","HG","HA","H*","HH","HI","HN","S*","SN","SS","IG","YG","NP","NT","IS"};
		//System.out.println("alph: "+alph.length);
		AlphabetContainer alphabetsRVD= model.getRVDAlphabet();

		Sequence talLeft = Sequence.create( alphabetsRVD, params.getLeftTALSequence(), "-" );
		Sequence talRight = Sequence.create( alphabetsRVD, params.getRightTALSequence(), "-" );

		int minDist = params.getMinimumDistance();
		int maxDist = params.getMaximumDistance();

		ResultList rl = new ResultList( params.getN() );
		ListResult rl2 = new ListResult( "", "", null );
	
		StringBuffer lastHeader = new StringBuffer();
		
		Pair<int[][],DataSet> curr = null;
		
		InfixMatchFinder singleFind = new InfixMatchFinder( null, Math.min( 8, Math.min( talLeft.getLength(), talRight.getLength() )+1 ), model );//TODO length
		
		//SimpleMatchFinder singleFind = new SimpleMatchFinder(null, model);
		//TALENTargetFinder finder = new TALENTargetFinder( null, model, singleFind );
		InfixTALENTargetFinder finder = new InfixTALENTargetFinder( null, model, singleFind );
		
		double bestTotal = model.getBestPossibleScore( talLeft, null ) + model.getBestPossibleScore( talRight, null );
		
		//double totalThresh = (bestTotal)/(talLeft.getLength()+talRight.getLength()+2) + Math.log( 0.4 );
		//double singleThresh1 = model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1) + Math.log( 0.3 );
		//double singleThresh2 = model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1) + Math.log( 0.3 );
	
		double th = params.getFilter();
		
		double totalThresh = model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1) + model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1) + 2*Math.log( th );
//		double singleThresh1 = totalThresh - model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1);
//		double singleThresh2 = totalThresh - model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1);
		double singleThresh1 = model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1) + Math.log( th*0.9 );
		double singleThresh2 = model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1) + Math.log( th*0.9 );
		
		singleFind.getPreps( talLeft, singleThresh1*(talLeft.getLength()+1) );//TODO
		singleFind.getPreps( talRight, singleThresh2*(talRight.getLength()+1) );
		
		//System.out.println(singleThresh1+" "+(model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1) + Math.log( 0.4 ))+" "+(model.getBestPossibleScore( talLeft, null )/(talLeft.getLength()+1)));
		//System.out.println(singleThresh2+" "+(model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1) + Math.log( 0.4 ))+" "+(model.getBestPossibleScore( talRight, null )/(talRight.getLength()+1)));
		
	/*	Sequence test1 = Sequence.create( DNAAlphabetContainer.SINGLETON, "TACTAGTATTATTACT" );
		Sequence test2 = Sequence.create( DNAAlphabetContainer.SINGLETON, "TGTCCTAAACAAAGAT" );
		
		System.out.println(model.getMatchString( talLeft, test1 )+" "+(model.getPartialLogScoreFor( talLeft, test1, 0, 0, test1.getLength() )/(talLeft.getLength()+1))+" "+model.getBestPossibleScore( talLeft, null ));
		System.out.println(model.getMatchString( talRight, test2 )+" "+(model.getPartialLogScoreFor( talRight, test2, 0, 0, test2.getLength() )/(talRight.getLength()+1))+" "+model.getBestPossibleScore( talRight, null ));
		
		System.out.println((model.getPartialLogScoreFor( talLeft, test1, 0, 0, test1.getLength() )+model.getPartialLogScoreFor( talRight, test2, 0, 0, test2.getLength() ))+" "+
				(model.getBestPossibleScore( talLeft, null )+model.getBestPossibleScore( talRight, null )));
*/		
		
		if(numThreads == 1){

			Worker worker = new Worker(this, finder.clone(),talLeft,talRight,totalThresh,singleThresh1,singleThresh2,minDist,maxDist,params.getNTerm1(), params.getNTerm2(), params.getHeteroOnly(), params.getN(),model,gff,params.getOutput(),0);
			
			while( (curr = readNextSequences( numThreads, null, dsRead, lastHeader )) != null){
				
				worker.setData(curr);
				worker.compute();
				
			}
			
			ResultList rl3 = worker.rl;
			for(int j=0;j<rl3.getNumberOfResults();j++){
				rl.add( rl3.list[j].getElement(), -rl3.list[j].getWeight() );
			}
		}else{




			available = new LinkedList<Integer>();

			long size = 0;

			Worker[] workers = new Worker[numThreads <=2 ? 1 : numThreads-2];
			for(int i=0;i<workers.length;i++){
				workers[i] = new Worker(this, finder.clone(),talLeft,talRight,totalThresh,singleThresh1,singleThresh2,minDist,maxDist,params.getNTerm1(), params.getNTerm2(), params.getHeteroOnly(), params.getN(),model,gff,params.getOutput(),i);
				Thread temp = (new Thread(workers[i]));
				temp.setDaemon( true );
				temp.start();
				available.add( i );
			}

			SequenceParser parser = null;
			if(numThreads > 1){
				parser = new SequenceParser( this );
				Thread temp = (new Thread( parser ));
				temp.setDaemon( true );
				temp.start();
			}

			while( (curr = readNextSequences( numThreads, parser, dsRead, lastHeader )) != null){
				//System.out.println("next seq");
				int a = 0;
				synchronized( available ) {
					a = available.size();
				}
				while(a == 0){
					synchronized(this){
						wait(1000);
					}
					synchronized( available ) {
						a = available.size();
					}
				}
				//if(available.size() > 0){
				int idx = 0;
				synchronized( available ) {
					idx = available.pop();
					if(workers[idx].exception){
						throw new RuntimeException( "Exception in one of the threads." );
					}
				}
				//System.out.println("avail: "+idx);
				workers[idx].setData( curr );
				//System.out.println("set");
				if(numThreads <= 2){
					//System.out.println("waiting");
					synchronized(this){
						wait();
					}
				}
				/*}else{
				System.out.println("waiting");
				synchronized(this){
					wait();
				}
			}*/
				//System.out.println("reading next");
			}
			int a = 0;
			synchronized( available ) {
				a = available.size();
			}
			while(a < workers.length){
				synchronized(this){
					//System.out.println("main waiting "+a+" "+workers.length);
					wait(1000);
				}
				synchronized( available ) {
					a = available.size();
				}
			}
			//System.out.println("stopping");
			if(numThreads > 1){
				parser.stop();
			}
			for(int i=0;i<workers.length;i++){
				workers[i].stop();
				ResultList rl3 = workers[i].rl;
				for(int j=0;j<rl3.getNumberOfResults();j++){
					rl.add( rl3.list[j].getElement(), -rl3.list[j].getWeight() );
				}
				size += workers[i].size;
			}

		}
		
		//TODO GFF
		
		if(params.getOutput() != Output.NONE){
			rl2 = (rl.gffToList()).toListResult( "Off target predictions ("+params.getOutput().name()+")", "", params.getOutput() == Output.GFF3 );
		}
		
		//System.out.println("input size: "+size);
		
		Result[] res = new Result[params.getNumberOfParameters()-2];
		for(int i=0;i<res.length;i++){
			res[i] = new CategoricalResult( params.getParameterAt( i+2 ).getName(), params.getParameterAt( i+2 ).getComment(), params.getParameterAt( i+2 ).getValue() == null ? "NA" : 
				(i==5 ? "min="+params.getMinimumDistance()+", max="+params.getMaximumDistance() :
					(i==6 ? "q="+params.getFilter() :
						params.getParameterAt( i+2 ).getValue().toString().replaceAll( "\n", "; " )
					)
				));
		}
		
		return new ListResult[]{new ListResult( "Off target predictions", "The off targets predicted for TALENs "+params.getLeftTALSequence()+" and "+params.getRightTALSequence()+" with distance between "+params.getMinimumDistance()+" and "+params.getMaximumDistance(), new ResultSet( res ), rl.resultsToArray() ),rl2};
	}

	
	private static String getFeatString( String id, int start, int end, GFFList gff ) {
		String featString = "";
		if(gff != null){
			GFFList sub = gff.getSubList( id );
			if(sub != null){
				ArrayList<GFFEntry> feats = sub.getEntriesOverlapping( start, end );
				ArrayList<GFFEntry> feats2 = new ArrayList<GFFParser.GFFEntry>();
				if(feats != null){
					ArrayList<GFFEntry> parents = new ArrayList<GFFParser.GFFEntry>();
					LinkedList<GFFEntry> stack = new LinkedList<GFFParser.GFFEntry>();
					for(int j=0;j<feats.size();j++){
						LinkedList<GFFEntry> temp = feats.get( j ).getParents();
						if(temp != null){
							stack.addAll( temp );
						}
					}
					while(stack.size() > 0){
						GFFEntry temp = stack.pop();
						parents.add( temp );
						LinkedList<GFFEntry> temp2 = temp.getParents();
						if(temp2 != null){
							stack.addAll( temp2 );
						}
					}
					for(int j=0;j<feats.size();j++){
						if(!parents.contains( feats.get( j ) ) && feats.get( j ).getParsedType() != GFFType.chromosome){
							feats2.add( feats.get( j ) );
						}
					}
					if(feats2.size() > 0){
						featString = feats2.toString();
						featString = featString.substring( 1, featString.length()-1 );
					}
				}
			}
		}
		return featString;
	}


	public Pair<int[][],DataSet> readNextSequences(int numThreads, SequenceParser parser, BufferedReader read, StringBuffer lastHeader) throws Exception {
		//System.out.println("started reading");
		String str = null;
		
		StringBuffer line = new StringBuffer();
		
		IntList starts = new IntList();
		IntList ends = new IntList();
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		if(parser != null){
			parser.setList( seqs );
		}
		
		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		int size = 0;
		
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
						
						if(parser != null){
							parser.add( header, seqStr.substring( start, end ) );
						}else{

							SequenceAnnotation annotation = new SequenceAnnotation( "unparsed comment line", "unparsed comment line", new CategoricalResult( "unparsed comment", "unparsed comment", header ) );
							Sequence seq = new SparseSequence( DNAAlphabetContainer.SINGLETON, seqStr.substring( start, end ) );
							seq = seq.annotate( false, annotation );
							seqs.add( seq );
						}
						size += end-start;
						starts.add( start );
						ends.add( seqStr.length()-end );
					}
					if(size > 1E7 || str == null){
						int s = 0;
						if(parser != null){
							synchronized( seqs ) {
								s = starts.length()-seqs.size();
							}
							while(s != 0){
								if(parser.exception){
									throw new RuntimeException("Exception in parser thread");
								}
								synchronized( this ) {
									wait();
								}
								synchronized( seqs ) {
									s = starts.length()-seqs.size();
									//System.out.println(starts.length()+" <-> "+seqs.size());
								}

							};
						}
						//System.out.println("read: "+header);
						//System.out.println("finished reading "+starts.length()+" "+ends.length()+" "+seqs.size());
						return new Pair<int[][],DataSet>(new int[][]{starts.toArray(),ends.toArray()},new DataSet( "", seqs ));
					}
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
	
	private static LinkedList<Integer> available;
	
	private static class SequenceParser implements Runnable{
		
		private LinkedList<Sequence> list;
		private FastTALENScanner caller;

		
		private LinkedList<String> headers;
		private LinkedList<String> seqStrs;
		private boolean stop;
		
		private boolean exception;
		
		public SequenceParser(FastTALENScanner caller){
			this.headers = new LinkedList<String>();
			this.seqStrs = new LinkedList<String>();
			this.caller = caller;
		}
		
		public void setList(LinkedList<Sequence> list){
			this.list = list;
		}
		
		public synchronized void stop(){
			stop = true;
			notify();
		}
		
		public void add(String header, String seqStr){
			synchronized( seqStrs ) {
				seqStrs.add( seqStr );
			}
			synchronized( headers ) {
				headers.add( header );
			}
			synchronized( this ) {
				notify();
			}
		}
		
		public void run(){
			exception = false;
			while(!stop){
				int s = 0;
				synchronized( headers ) {
					s = headers.size();
				}
				if(s > 0){
					try{
						String header = null;
						synchronized( headers ) {
							header = headers.pop();
						}
						String seqStr = null;
						synchronized( seqStrs ) {
							seqStr = seqStrs.pop();
						}
						SequenceAnnotation annotation = new SequenceAnnotation( "unparsed comment line", "unparsed comment line", new CategoricalResult( "unparsed comment", "unparsed comment", header ) );

						Sequence seq = new SparseSequence( DNAAlphabetContainer.SINGLETON, seqStr );
						seq = seq.annotate( false, annotation );
						synchronized( list ) {
							list.add( seq );
						}
					}catch(Exception e){
						exception = true;
						e.printStackTrace();
						throw new RuntimeException();
					}
				}else{
					synchronized(this){
						synchronized( caller ) {
							caller.notify();
						}

						try {
							wait();
						} catch ( InterruptedException e ) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					
				}
			}
		}
		
	}
	
	private static class Worker implements Runnable{
		
		private FastTALENScanner caller;
		private long size;
		private Pair<int[][],DataSet> curr;
		private InfixTALENTargetFinder finder;
		private Sequence talLeft, talRight;
		private double totalThresh, singleThresh1, singleThresh2;
		private int minDist, maxDist, N;
		private ResultList rl;
		private TALgetterDiffSM model;
		private Output output;
		private GFFList gff;
		private boolean nTerm1Ext;
		private boolean nTerm2Ext;
		private boolean onlyhetero;
		
		private boolean stop;
		private boolean exception;
		
		private int index;
		
		
		
		/**
		 * @param size
		 * @param finder
		 * @param talLeft
		 * @param talRight
		 * @param totalThresh
		 * @param singleThresh1
		 * @param singleThresh2
		 * @param minDist
		 * @param maxDist
		 * @param n
		 * @param model
		 * @param output
		 * @param index
		 */
		public Worker( FastTALENScanner caller, InfixTALENTargetFinder finder, Sequence talLeft, Sequence talRight, double totalThresh,
						double singleThresh1, double singleThresh2, int minDist, int maxDist, boolean nTerm1, boolean nTerm2, boolean onlyhetero, int N, TALgetterDiffSM model,
						GFFList gff, Output output, int index ) {
			this.size = 0;
			this.finder = finder;
			this.talLeft = talLeft;
			this.talRight = talRight;
			this.totalThresh = totalThresh;
			this.singleThresh1 = singleThresh1;
			this.singleThresh2 = singleThresh2;
			this.minDist = minDist;
			this.maxDist = maxDist;
			this.nTerm1Ext = nTerm1;
			this.nTerm2Ext = nTerm2;
			this.onlyhetero = onlyhetero;
			this.N = N;
			this.model = model;
			this.gff = gff;
			this.index = index;
			this.rl = new ResultList( N );
			this.output = output;
			this.caller = caller;
		}

		public void setData(Pair<int[][],DataSet> data){
			synchronized(this){
				this.curr = data;
				//System.out.println("worker "+index+" is notified");
				notify();
			}
		}
		
		public void run(){
			try {
				synchronized(this){
					//System.out.println("worker "+index+" waits");
					wait();
				}
			} catch ( InterruptedException e ) {
				e.printStackTrace();
			}
			while(!stop){
				
				
				compute();
				if(exception){
					return;
				}
				synchronized(this){
					synchronized( available ) {
						available.add( index );
					}

					try{
						synchronized(caller){
							//System.out.println("worker "+index+" notifies");
							caller.notify();
						}
					}catch(IllegalMonitorStateException ex){
						ex.printStackTrace();
					}
					try {
						//synchronized(this){
							//System.out.println("worker "+index+" waits");
							wait();
						//}
					} catch ( InterruptedException e ) {
						e.printStackTrace();
					}
				}
			}
		}
		
		public void stop(){
			//System.out.println("stop worker "+index);
			stop = true;
			synchronized(this){
				notify();
			}
		}
		
		private void compute(){
			exception = false;
			int[][] offsets = curr.getFirstElement();
			DataSet dat = curr.getSecondElement();

			for(int i=0;i<dat.getNumberOfElements();i++){
				size += dat.getElementAt( i ).getLength();
			}

			finder.setDataSet( dat );	
			//thresh = Double.NEGATIVE_INFINITY;


			//System.out.println("started for: "+dat.getElementAt( 0 ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue());

			try{
				LimitedSortedList<TALENMatch> list = finder.getTALENMatches( talLeft, talRight, totalThresh, singleThresh1, singleThresh2, minDist, maxDist, N, true, nTerm1Ext, nTerm2Ext, onlyhetero );


				//System.out.println(list.getBestScore()+" <-> "+list.getWorstScore()+" ? "+rl.getWorstScore());

				ComparableElement<TALENMatch, Double>[] matches = list.getSortedList();



				int i=matches.length-1;
				while(i >= 0 && rl.better( matches[i].getWeight() )){

					TALENMatch match = matches[i].getElement();

					double score = matches[i].getWeight();


					Sequence tal1 = match.getTal1();
					Sequence tal2 = match.getTal2();
					
					

					Match match1 = match.getMatch1();
					Match match2 = match.getMatch2();


					Sequence seq1 = dat.getElementAt( match1.getSeqIdx() );
					Sequence seq2 = dat.getElementAt( match2.getSeqIdx() );

					int pos1 = match1.getSeqPos();
					int pos2 = match2.getSeqPos();

					String id = (String) seq1.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue();

					id = id.trim();
					int space = id.indexOf( ' ' );
					if(space >= 0){
						id = id.substring( 0, space );
					}

					int start = pos1+offsets[0][match1.getSeqIdx()]+1;
					int end = pos2+tal2.getLength()+offsets[0][match2.getSeqIdx()]+1;
					
					if(pos1 > pos2){
						start = pos2+offsets[0][match2.getSeqIdx()]+1;
						end = pos1+tal1.getLength()+offsets[0][match1.getSeqIdx()]+1;
					}

					byte cat = match.getCat();
					
					String featString = null;
					/*if(cat == 3){
						featString = getFeatString(id,pos2+offsets[0][match2.getSeqIdx()]+1,pos1+tal1.getLength()+offsets[0][match1.getSeqIdx()]+1,gff);
					}else if(cat == 2){
						featString = getFeatString( id, start, end, gff );//TODO
					}else if(cat == 3){
						featString = getFeatString( id, start, end, gff );//TODO
					}else{*/
						featString = getFeatString(id,start,end,gff);
					//}

					String archString = (tal1.equals( talLeft ) ? "first" : "second") + "-" + (tal2.equals( talLeft ) ? "first" : "second");

					int dist = pos2-(pos1+tal1.getLength()+1);
					if(cat == 3){
						dist = -(pos1 - pos2 - (tal2.getLength()+1)); 
					}

					Sequence sub1 = seq1.getSubSequence( pos1, tal1.getLength()+1 );
					Sequence sub2 = seq2.reverseComplement().getSubSequence( seq2.getLength()-pos2-(tal2.getLength()+1), tal2.getLength()+1 );
					Sequence between = null;
					if(cat == 3){
						between = seq1.getSubSequence( pos2+tal2.getLength()+1,-dist );
					}else if(cat == 1 || cat==2){
						if(match1.isRc()){
							sub1 = seq1.reverseComplement().getSubSequence( seq1.getLength()-pos1-tal1.getLength()-1, tal1.getLength()+1 );
							dist = -pos1+pos2-tal2.getLength()-1;
						}
						if(!match2.isRc()){
							sub2 = seq2.getSubSequence( pos2, tal2.getLength()+1 );
							dist = pos2-(pos1+tal1.getLength()+1);
						}
					//	System.out.println(pos1+" "+pos2+" "+dist+"; "+match1.getSeqPos()+" "+match2.getSeqPos()+"; "+match1.isRc()+" "+match2.isRc()+"; "+seq1.getLength());
						between = seq1.getSubSequence( pos1+tal1.getLength()+1,dist );
					}else{
						between = seq1.getSubSequence( pos1+tal1.getLength()+1,dist );
					}
					
					String site = null;
					if(cat == 3){
						site = sub2.reverseComplement().toString().toUpperCase()+between.toString().toLowerCase()+sub1.toString().toUpperCase();
					}else if(cat == 1 || cat == 2){
						if(match1.isRc()){
							site = sub1.reverseComplement().toString().toUpperCase()+between.toString().toLowerCase()+sub2.reverseComplement().toString().toUpperCase();
						}else{
							site = sub1.toString().toUpperCase()+between.toString().toLowerCase()+sub2.toString().toUpperCase();
						}
					}else{
						site = sub1.toString().toUpperCase()+between.toString().toLowerCase()+sub2.reverseComplement().toString().toUpperCase();
					}					
					
					String str1 = model.getMatchString( tal1, sub1 );
					String str2 = model.getMatchString( tal2, sub2 );

				//	System.out.println(score+" <-> "+((model.getPartialLogScoreFor( tal1, sub1, 0, 0, sub1.getLength() )/(sub1.getLength()))+(model.getPartialLogScoreFor( tal2, sub2, 0, 0, sub2.getLength() )/(sub2.getLength()))));
					
					Result[] tr = new Result[]{
					                           new CategoricalResult( "ID", "", id ),
					                           new NumericalResult( "Position 1", "", pos1+offsets[0][match1.getSeqIdx()] ),
					                           new NumericalResult( "Position 2", "", pos2+offsets[0][match2.getSeqIdx()] ),
					                           new NumericalResult( "Distance", "", dist ),
					                           new CategoricalResult( "Sequence 1", "", sub1.toString() ),
					                           new CategoricalResult( "Matches 1", "", str1 ),
					                           new CategoricalResult( "Sequence 2", "", sub2.toString() ),
					                           new CategoricalResult( "Matches 2", "", str2 ),
					                           new CategoricalResult( "Architecture", "", archString ),
					                           new CategoricalResult( "Full site", "", site ),
					                           new NumericalResult( "Score", "", score ),
					                           new CategoricalResult( "Features", "", featString )

					};

					if(output != Output.NONE){
						String predId = id+"_"+start+"_"+end+"_"+archString;
						GFFEntry[] ens = new GFFEntry[3];
						ens[0] = new GFFEntry( id, "TALoffer", "TALEN_dimer_target", start, end, score, Strand.UNKNOWN, -1, "ID="+predId+"; group="+predId );//TODO
						ens[1] = new GFFEntry( id, "TALoffer", "TALEN_monomer_target", pos1+offsets[0][match1.getSeqIdx()]+1, pos1+offsets[0][match1.getSeqIdx()]+tal1.getLength()+1, score, match1.isRc() ? Strand.REVERSE : Strand.FORWARD, -1, "Parent="+predId+"; group="+predId+"f" );
						ens[2] = new GFFEntry( id, "TALoffer", "TALEN_monomer_target", pos2+offsets[0][match2.getSeqIdx()]+1, pos2+offsets[0][match2.getSeqIdx()]+tal2.getLength()+1, score, match1.isRc() ? Strand.REVERSE : Strand.FORWARD, -1, "Parent="+predId+"; group="+predId+"r" );

						rl.add( new Object[]{new ResultSet( tr ),ens}, score );

					}else{
						rl.add( new ResultSet( tr ), score );
					}
					/*if(output != Output.NONE){
						String predId = id+"_"+start+"_"+end+"_"+archString;

						en.add( new GFFEntry( id, "TALoffer", "TALEN_dimer_target", start, end, score, Strand.UNKNOWN, -1, "ID="+predId+"; group="+predId ) );
						en.add( new GFFEntry( id, "TALoffer", "TALEN_monomer_target", pos1+offsets[0][match1.getSeqIdx()]+1, pos1+offsets[0][match1.getSeqIdx()]+tal1.getLength()+1, score, Strand.FORWARD, -1, "Parent="+predId+"; group="+predId+"f" ) );
						en.add( new GFFEntry( id, "TALoffer", "TALEN_monomer_target", pos2+offsets[0][match2.getSeqIdx()]+1, pos2+offsets[0][match2.getSeqIdx()]+tal2.getLength()+1, score, Strand.REVERSE, -1, "Parent="+predId+"; group="+predId+"r" ) );


					}*///TODO


					i--;
				}
				if(rl.getWorstScore() > totalThresh){
					//System.out.println("setting total, before "+totalThresh+" after: "+rl.getWorstScore()+" quot "+(totalThresh/rl.getWorstScore()));
					totalThresh = rl.getWorstScore();
				}
			}catch(Exception e){
				exception = true;
				e.printStackTrace();
				throw new RuntimeException();
			}
		}
		
	} 
	
}
