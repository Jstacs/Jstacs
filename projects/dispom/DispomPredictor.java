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
package projects.dispom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.PermutedSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder.RandomSeqType;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;

/**
 * Discriminative de-novo position distribution and motif finder.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class DispomPredictor {
	
	
	private static final String[] PREFIX = {
	                                       DispomPredictorParameterSet.HOME,DispomPredictorParameterSet.IGNORE_CHAR,DispomPredictorParameterSet.FG,DispomPredictorParameterSet.BG,
	                                       DispomPredictorParameterSet.XML_PATH,DispomPredictorParameterSet.P_VALUE,DispomPredictorParameterSet.ONEHIST
	        }; 
	
	private static DataSet getSample( AlphabetContainer con, String fileName, char ignore ) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		return new DataSet( con, 
				new SparseStringExtractor( fileName, ignore ) 
		);
	}
	
	/**
	 * This is the main of Dispom that starts the program. 
	 * 
	 * @param args the arguments for Dispom. Each argument has the form <code>name=value</code>.
	 * 
	 * @throws Exception if something went wrong.
	 */
	public static void main( String[] args ) throws Exception {
		
		//parse parameters
		ParameterSetTagger params = new ParameterSetTagger( PREFIX, new DispomPredictorParameterSet() );
		params.fillParameters( "=", args );
		System.out.println( "parameters:" );
		System.out.println( params );
		System.out.println("_________________________________");
		if( !params.hasDefaultOrIsSet() ) {
			System.out.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}
		
		//load data
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		String home = params.getValueFromTag( DispomPredictorParameterSet.HOME, String.class );
		char ignore = params.getValueFromTag( DispomPredictorParameterSet.IGNORE_CHAR, Character.class );

		int anz = 1;
		if( params.isSet( DispomPredictorParameterSet.BG ) ) {
			anz = 2;
		}
		
		DataSet[] data = new DataSet[anz];
		data[0] = getSample( con, home + File.separatorChar + params.getValueFromTag( DispomPredictorParameterSet.FG, String.class ), ignore );	
		if( anz > 1 ) {
			data[1] = getSample( con, home + File.separatorChar + params.getValueFromTag( DispomPredictorParameterSet.BG, String.class ), ignore );				
		}

	
		//create simple weights
		double[][] weights = new double[anz][];
		for( int i = 0; i < data.length; i++ ) {
			weights[i] = new double[data[i].getNumberOfElements()];
			Arrays.fill( weights[i], 1 );
			System.out.println( i + "\t# = " + data[i].getNumberOfElements() + "\tlength = " + data[i].getElementLength() + "\t" + data[i].getAnnotation() );
		}
		int sl = data[0].getElementLength();		
		
		String fName = params.getValueFromTag( DispomPredictorParameterSet.XML_PATH, String.class );
		if( !fName.endsWith( ".xml" ) ) {
			fName = fName + ".xml";
		}
		//restore classifier
		GenDisMixClassifier cl = new GenDisMixClassifier( XMLParser.extractForTag( FileManager.readFile( new File(fName) ),"classifier")  );
		DifferentiableStatisticalModel[] bestNSF = ArrayHandler.cast( DifferentiableStatisticalModel.class, cl.getDifferentiableSequenceScores() );

		// show
		System.out.println("_________________________________");
		MotifDiscoverer md = (MotifDiscoverer) bestNSF[0];
		System.out.println( md );
		System.out.println( "result: " + cl.getLastScore() );
		
		// predict BSs
		
		System.out.println("_________________________________");
		SignificantMotifOccurrencesFinder smof;
		if( anz > 1 ) {
			smof = new SignificantMotifOccurrencesFinder( md, data[1], null, params.getValueFromTag( DispomPredictorParameterSet.P_VALUE, Double.class ) );			
		} else if(params.getValueFromTag( DispomPredictorParameterSet.ONEHIST, Boolean.class )){
			Random r = new Random();
			Sequence[] seqs = new Sequence[1000];
			int numFg = data[0].getNumberOfElements();
			for( int n = 0; n < seqs.length; n++ ) {
				seqs[n] = new PermutedSequence( data[0].getElementAt( r.nextInt( numFg ) ) );
			}
			DataSet bg = new DataSet("permuted",seqs);
			smof = new SignificantMotifOccurrencesFinder( md, bg, null, params.getValueFromTag( DispomPredictorParameterSet.P_VALUE, Double.class ) );
		} else {
			smof = new SignificantMotifOccurrencesFinder( md, RandomSeqType.PERMUTED, false, 1000, params.getValueFromTag( DispomPredictorParameterSet.P_VALUE, Double.class ) );
		}
		Sequence seq, site, adjusted;
		MotifAnnotation ma;
		LinkedList<Sequence> list = new LinkedList<Sequence>();
		SequenceAnnotation[] seqAn;
		int start, end;
		for( int m = 0; m < md.getNumberOfMotifs(); m++ ) {
			list.clear();
			System.out.println();
			System.out.println( "predictions for motif " + m );
			System.out.println( "sequence\tposition\tstrand\tbinding site\tadjusted binding site\tp-value" );
			System.out.println( "------------------------------------------------------------------------");
			DataSet annotated = smof.annotateMotif(data[0], m);
			for( int j, i = 0; i < annotated.getNumberOfElements(); i++ ) {
				seqAn = data[0].getElementAt( i ).getAnnotation();
				start = seqAn == null ? 0 : seqAn.length;
				seq = annotated.getElementAt( i );
				seqAn = seq.getAnnotation();
				
				end = (seqAn == null ? 0 : seqAn.length);
				
				if( end - start > 0 ) {
					for( j = start; j < end; j++ ) {
						ma = (MotifAnnotation) seqAn[j];
						site = seq.getSubSequence( ma.getPosition(), ma.getLength() );
						adjusted = site;
						if(ma.getStrandedness() == Strand.REVERSE){
							adjusted = adjusted.reverseComplement();
						}
						list.add( adjusted );
						System.out.println( i + "\t" + ma.getPosition() + "\t" + ma.getStrandedness() + "\t" + site + "\t" + adjusted + "\t" + ma.getAnnotations()[1].getValue() );
					}
				}
			}
			if( list.size() >  0 ) {
				double[][] pfm = getPFM( new DataSet("",list.toArray(new Sequence[0])) );
				System.out.println( "------------------------------------------------------------------------");
				System.out.println("Position frequency matrix of sites: ");
				for(int i=0;i<pfm.length;i++){
					System.out.print("\t"+i);
				}
				System.out.println();
				for(int j=0;j<pfm[0].length;j++){
					System.out.print(DNAAlphabet.SINGLETON.getSymbolAt( j ));
					for(int i=0;i<pfm.length;i++){
						System.out.print("\t"+pfm[i][j]);
					}
					System.out.println();
				}
				System.out.println();
				System.out.println( "------------------------------------------------------------------------");
				System.out.println("Position weight matrix of sites: ");
				for(int i=0;i<pfm.length;i++){
					System.out.print("\t"+i);
				}
				System.out.println();
				double s = list.size();
				for(int j=0;j<pfm[0].length;j++){
					System.out.print(DNAAlphabet.SINGLETON.getSymbolAt( j ));
					for(int i=0;i<pfm.length;i++){
						System.out.print("\t"+(pfm[i][j]/s));
					}
					System.out.println();
				}
			}
		}
	}
	
	public static double[][] getPFM( DataSet data ) {
		if( data == null ) {
			return null;
		} else {
			double[][] pfm = new double[data.getElementLength()][];
			AlphabetContainer con = data.getAlphabetContainer();
			for( int l = 0; l < pfm.length; l++ ) {
				pfm[l] = new double[(int)con.getAlphabetLengthAt( l )];
			}
			Sequence seq;
			for( int n = 0; n < data.getNumberOfElements(); n++ ) {
				seq = data.getElementAt( n );
				for( int l = 0; l < pfm.length; l++ ) {
					pfm[l][seq.discreteVal( l )]++;
				}
			}
			return pfm;
		}
	}
}



/**
 * This class is a container for all parameters of Dispom. It also parses the parameter from Strings.
 *  
 * @author Jens Keilwagen
 */
class DispomPredictorParameterSet extends ParameterSet {

	public static final String HOME = "home";
	public static final String IGNORE_CHAR = "ignore";
	public static final String FG = "fg";
	public static final String BG = "bg";
	public static final String XML_PATH = "xml";
	public static final String P_VALUE = "p-val";
	public static final String ONEHIST = "one-histogram";
	
	public DispomPredictorParameterSet() throws Exception {
		super(
				new SimpleParameter( DataType.STRING, "home directory", "the path to the data directory", true, "./" ),
				new SimpleParameter( DataType.CHAR, "the ignore char for the data files", "the char that is used to mask comment lines in data files, e.g., '>' in a FASTA-file", true, '>' ),
				new SimpleParameter( DataType.STRING, "foreground file", "the file name of the foreground data file (the file containing sequences which are expected to contain binding sites of a common motif)", true ),
				new SimpleParameter( DataType.STRING, "background file", "the file name of the background data file", false ),
				new SimpleParameter( DataType.STRING, "classifier xml-file", "the file name of the xml file containing the classifier", true, "./classifier.xml" ),
				new SimpleParameter( DataType.DOUBLE, "p-value", "a p-value for predicting binding sites", true, new NumberValidator<Double>(0d,1d), new Double(1E-4) ),
				new SimpleParameter( DataType.BOOLEAN, "one histogram","If no background file is specificed, p-values are computed either using a joint histogram (true), or a sequence-wise histogram (false)",false,new Boolean( true ) )
		);
	}
}