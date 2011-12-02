package de.jstacs.models.hmm.transitions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.models.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class can be used in any {@link de.jstacs.models.hmm.AbstractHMM} allowing to use gradient based or sampling training algorithm.  
 * 
 * @author Jens Keilwagen
 */
public class HigherOrderTransition extends BasicHigherOrderTransition implements DifferentiableTransition, SamplingTransition {

	private static final String PREFIX = "samplingHOTransition-";
	
	private int offset;
	
	/**
	 * The files for saving the parameters during the sampling.
	 */
	protected File[] paramsFile;

	/**
	 * The counter for the sampling steps of each sampling.
	 */
	protected int[] counter;

	/**
	 * The index of the current sampling.
	 */
	protected int samplingIndex;

	/**
	 * The writer for the <code>paramsFile</code> in a sampling.
	 */
	protected BufferedWriter writer;

	/**
	 * The reader for the <code>paramsFile</code> after a sampling.
	 */
	protected BufferedReader reader;
	
	/**
	 * Reference used during parsing the parameters.
	 */
	private double[] params;
	
	/**
	 * The main constructor.
	 * 
	 * @param transitions the {@link TransitionElement}s for the internal use
	 * @param isSilent an array indicating for each state whether it is silent or not
	 * 
	 * @throws Exception if an error occurs during checking the {@link TransitionElement}s and creating internal fields 
	 */
	public HigherOrderTransition( boolean[] isSilent, TransitionElement... transitions ) throws Exception {
		super( isSilent, transitions );
		init();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link HigherOrderTransition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HigherOrderTransition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 *             
	 */
	public HigherOrderTransition( StringBuffer xml ) throws NonParsableException {
		super( xml );
		setParameterOffset();
	}
	
	private void init() {
		int n = 0;
		for( int t = 0; t < transitions.length; t++ ) {
			n += transitions[t].getNumberOfParameters();
		}
		params = new double[n];
	}
	
	private static final String XML_TAG = "HigherOrderTransition";
	
	protected String getXMLTag() {
		return XML_TAG;
	}
	
	@Override
	protected void appendFurtherInformation( StringBuffer xml ) {
		if( writer != null ) {
			throw new RuntimeException( "could not parse the model to XML while sampling" );
		}
		XMLParser.appendObjectWithTags( xml, offset, "offset" );
		if( paramsFile != null ) {
			try {
				XMLParser.appendObjectWithTags( xml, counter, "counter" );
				String content;
				for( int i = 0; i < paramsFile.length; i++ ) {
					if( paramsFile[i] != null ) {
						content = FileManager.readFile( paramsFile[i] ).toString();
					} else {
						content = "";
					}
					XMLParser.appendObjectWithTagsAndAttributes( xml, content, "fileContent", "pos=\"" + i + "\"" );
				}
			} catch ( IOException e ) {
				RuntimeException r = new RuntimeException( e.getMessage() );
				r.setStackTrace( e.getStackTrace() );
				throw r;
			}
		}
	}

	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		offset = (Integer) XMLParser.extractObjectForTags( xml, "offset" );
		if( XMLParser.hasTag( xml, "counter", null, null ) ) {
			counter = (int[]) XMLParser.extractObjectForTags( xml, "counter" );
			paramsFile = new File[counter.length];
			try {
				String content;
				Map<String,String> filter = new TreeMap<String, String>();
				for( int i = 0; i < paramsFile.length; i++ ) {
					filter.clear();
					filter.put( "pos", ""+i );
					content = XMLParser.extractObjectAndAttributesForTags( xml, "fileContent", null, filter, String.class );
					if( !content.equalsIgnoreCase( "" ) ) {
						paramsFile[i] = File.createTempFile( PREFIX, ".dat", null );
						FileManager.writeFile( paramsFile[i], new StringBuffer( content ) );
					}
				}
			} catch ( IOException e ) {
				NonParsableException n = new NonParsableException( e.getMessage() );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
		} else {
			counter = null;
			paramsFile = null;
		}
		writer = null;
		reader = null;
		init();
	}
	
	
	public HigherOrderTransition clone() throws CloneNotSupportedException {
		HigherOrderTransition clone = (HigherOrderTransition) super.clone();
		if( writer != null ) {
			throw new CloneNotSupportedException( "sampling was not stopped before" );
		} else {
			clone.writer = null;
		}
		clone.reader = null;
		clone.params = params.clone();
		
		if( paramsFile != null ) {
			try {
				clone.paramsFile = new File[paramsFile.length];
				clone.counter = new int[paramsFile.length];
				for( int i = 0; i < paramsFile.length; i++ ) {
					if( paramsFile[i] != null ) {
						clone.paramsFile[i] = File.createTempFile( PREFIX, ".dat", null );
						FileManager.copy( paramsFile[i].getAbsolutePath(), clone.paramsFile[i].getAbsolutePath() );
						clone.counter[i] = counter[i];
					}
				}
			} catch ( IOException e ) {
				CloneNotSupportedException c = new CloneNotSupportedException( e.getMessage() );
				c.setStackTrace( e.getStackTrace() );
				throw c;
			}

		}
		return clone;
	}

	@Override
	public void fillParameters( double[] params ) {
		fillParameters( params, offset );
	}
	
	/**
	 * This method allows to fill the current parameters using a specific offset.
	 * 
	 * @param params the parameters
	 * @param offset the offset indicating the start position
	 * 
	 * @see #fillParameters(double[])
	 * @see TransitionElement#fillParameters(double[], int)
	 */
	protected void fillParameters( double[] params, int offset ) {
		for( int o = offset, t = 0; t < transitions.length; t++ ) {
			o = ((TransitionElement)transitions[t]).fillParameters( params, o );
		}
	}
	
	@Override
	public int setParameterOffset( int offset ) {
		this.offset = offset;
		return setParameterOffset();
	}
	
	/**
	 * This method allows to set the parameter offset in each internally used {@link TransitionElement}.
	 * 
	 * @return the new parameter offset
	 * 
	 * @see #setParameterOffset(int)
	 * @see TransitionElement#setParameterOffset(int)
	 */
	protected int setParameterOffset() {
		int o = this.offset;
		for( int t = 0; t < transitions.length; t++ ) {
			o = ((TransitionElement)transitions[t]).setParameterOffset( o );
		}
		return o;
	}

	@Override
	public void setParameters( double[] params, int start ) {
		setParams( params, start+offset );
	}
	
	/**
	 * This method allows to set the new parameters using a specific offset.
	 * 
	 * @param params the parameters
	 * @param start the offset indicating the start position
	 * 
	 * @see #setParameters(double[], int)
	 * @see TransitionElement#setParameters(double[], int)
	 */
	protected void setParams( double[] params, int start ) {
		for( int s = start, t = 0; t < transitions.length; t++ ) {
			s = ((TransitionElement)transitions[t]).setParameters( params, s );
		}
	}
	
	@Override
	public void addGradientForLogPriorTerm( double[] gradient, int start ) {
		for( int t = 0; t < transitions.length; t++ ) {
			((TransitionElement)transitions[t]).addGradientForLogPriorTerm( gradient, start );
		}
	}

	@Override
	public double getLogScoreAndPartialDerivation( int layer, int index, int childIdx, IntList indices, DoubleList partDer, Sequence sequence, int sequencePosition ) {
		return ((TransitionElement)transitions[ getTransitionElementIndex( layer, index ) ]).getLogScoreAndPartialDerivation( childIdx, indices, partDer, sequence, sequencePosition );
	}

	@Override
	public void initForSampling( int starts ) throws IOException {
		for( int t = 0; t < transitions.length; t++ ) {
			if( ((TransitionElement)transitions[t]).getMinimalHyperparameter() <= 0 ) {
				throw new IllegalArgumentException( "All hyper-parameters must have a value > 0." );
			}
		}
		
		if( paramsFile != null && paramsFile.length == starts ) {
			FileOutputStream o;
			for( int i = 0; i < starts; i++ ) {
				if( paramsFile[i] != null ) {
					o = new FileOutputStream( paramsFile[i] );
					o.close();
				}
				counter[i] = 0;
			}
		} else {
			deleteParameterFiles();
			paramsFile = new File[starts];
			counter = new int[starts];
		}

	}
	
	private void deleteParameterFiles() {
		if( paramsFile != null ) {
			for( int i = 0; i < paramsFile.length; i++ ) {
				if( paramsFile[i] != null ) {
					paramsFile[i].delete();
				}
			}
		}
	}
	
	/* 
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		if( writer != null ) {
			writer.close();
		}
		if( reader != null ) {
			reader.close();
		}
		deleteParameterFiles();
		super.finalize();
	}
	
	@Override
	public void extendSampling( int sampling, boolean append ) throws IOException {
		if( paramsFile[sampling] == null ) {
			paramsFile[sampling] = File.createTempFile( PREFIX, ".dat", null );
			//System.out.println( paramsFile[start].getAbsolutePath() );
		} else {
			if( append ) {
				parseParameterSet( sampling, counter[sampling] - 1 );
				reader.close();
				reader = null;
			} else {
				counter[sampling] = 0;
			}
		}
		writer = new BufferedWriter( new FileWriter( paramsFile[sampling], append ) );
		samplingIndex = sampling;
	}

	@Override
	public boolean isInSamplingMode() {
		return writer != null;
	}

	@Override
	public boolean parseNextParameterSet() {
		if( writer != null ) {
			return false;
		}
		String str = null;
		try {
			str = reader.readLine();
		} catch ( IOException e ) {} finally {
			if( str == null ) {
				return false;
			}
		}

		parse( str );
		return true;
	}

	@Override
	public boolean parseParameterSet( int sampling, int n ) throws IOException {
		String str;
		if( reader != null ) {
			reader.close();
		}
		reader = new BufferedReader( new FileReader( paramsFile[sampling] ) );
		while( ( str = ( reader.readLine() ) ) != null ) {
			if( Integer.parseInt( str.substring( 0, str.indexOf( "\t" ) ) ) == n ) {
				parse( str );
				return true;
			}
		}
		return false;
	}
	
	private void parse( String str ) {
		String[] strArray = str.split( "\t" );
		int offset = 1;
		for( int p = 0; p < params.length; p++ ) {
			params[p] = Double.parseDouble(strArray[offset++]);
		}		
		setParams( params, 0 );
	}

	@Override
	public void samplingStopped() throws IOException {
		if( writer != null ) {
			writer.close();
			writer = null;
		}
	}

	@Override
	public void acceptParameters() throws IOException {
		writer.write( "" + ( counter[samplingIndex]++ ) );
		fillParameters( params, 0 );
		for( int p = 0; p < params.length; p++ ) {
			writer.write( "\t" + params[p] );
		}
		writer.newLine();
		writer.flush();
	}
	
	public double getLogPosteriorFromStatistic() {
		double logPost = 0;
		for( int t = 0; t < transitions.length; t++ ) {
			logPost += ((TransitionElement)transitions[t]).getLogPosteriorFromStatistic();
		}
		return logPost;
	}

	@Override
	public int getSizeOfEventSpace( int index ) {
		int off = this.offset;
		for(int i=0;i<transitions.length;i++){
			int num = ((TransitionElement)transitions[i]).getNumberOfParameters();
			if(index >= off && index < off + num){
				return num;
			}
			off += num;
		}
		return 0;
	}

	@Override
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list ) {
		int off = this.offset;
		for(int i=0;i<transitions.length;i++){
			int[] idxs = new int[((TransitionElement)transitions[i]).getNumberOfParameters()];
			for(int j=0;j<idxs.length;j++){
				idxs[j] = j + off + parameterOffset;
			}
			list.add( idxs );
			off += idxs.length;
		}
	}
	
	
}