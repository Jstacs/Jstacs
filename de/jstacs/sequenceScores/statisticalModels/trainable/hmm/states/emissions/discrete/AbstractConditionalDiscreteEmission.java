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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SamplingEmission;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.random.DiMRGParams;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jstacs.utils.random.FastDirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;


/**
 * The abstract super class of discrete emissions.
 * 
 * @author Jens Keilwagen, Michael Scharfe, Jan Grau
 */
public abstract class AbstractConditionalDiscreteEmission implements SamplingEmission, DifferentiableEmission {

	private double[] colors;
	
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
	 * The offset of the parameter indexes
	 */
	protected int offset;
	
	/**
	 * The alphabet of the emissions
	 */
	protected AlphabetContainer con;
	
	/**
	 * The parameters of the emission
	 */
	protected double[][] params;
	
	/**
	 * The parameters transformed to probabilites
	 */
	protected double[][] probs;
	
	/**
	 * The hyper-parameters for the prior on the parameters
	 */
	protected double[][] hyperParams;
	
	/**
	 * The array for storing the statistics for 
	 * each parameter
	 */
	protected double[][] statistic;
	
	/**
	 * The array for storing the gradients for
	 * each parameter
	 */
	protected double[][] grad;
	
	/**
	 * The log-normalization constants for each condition
	 */
	protected double[] logNorm;
	
	/**
	 * The equivalent sample sizes for each condition
	 */
	protected double[] ess;

	/**
	 * The hyper-parameters for initializing the parameters
	 */
	protected double[][] initHyperParams;
	
	private String shape;
	
	private boolean linear;
	
	private double logUniform;

	/**
	 * Returns the hyper-parameters for all parameters and a given ess.
	 * The equivalent sample size is distributed evenly across all parameters
	 * @param ess the equivalent sample size
	 * @param numConditions the number of conditions
	 * @param numEmissions the number of emissions, assumed to be equal for all conditions
	 * @return hyper-parameters for all parameters
	 */
	protected static double[][] getHyperParams(double ess, int numConditions, int numEmissions){
		double[] ess2 = new double[numConditions];
		Arrays.fill( ess2, ess/(double)numConditions );
		return getHyperParams( ess2, numEmissions );
	}
	
	private static double[][] getHyperParams( double[] ess, int number ){
		double[][] res = new double[ess.length][number];
		for(int i=0;i<res.length;i++){
			Arrays.fill( res[i], ess[i]/(double) number );
		}
		return res;
	}
	
	/**
	 * This is a simple constructor for a {@link AbstractConditionalDiscreteEmission} based on the equivalent sample size.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param numberOfConditions the number of conditions
	 * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 * 
	 * @see #AbstractConditionalDiscreteEmission(AlphabetContainer, double[][])
	 */
	protected AbstractConditionalDiscreteEmission( AlphabetContainer con, int numberOfConditions, double ess ) throws WrongAlphabetException {
		this( con, numberOfConditions, ess, ess );
	}

	protected AbstractConditionalDiscreteEmission( AlphabetContainer con, int numberOfConditions, double ess, double initEss ) throws WrongAlphabetException {
		this( con, getHyperParams( ess, numberOfConditions, (int) con.getAlphabetLengthAt( 0 )), getHyperParams( initEss, numberOfConditions, (int) con.getAlphabetLengthAt( 0 )) );
	}
	
	/**
	 * This is a simple constructor for a {@link AbstractConditionalDiscreteEmission} defining the individual hyper parameters.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param hyperParams the individual hyper parameters for each parameter
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 *  
	 * @see #AbstractConditionalDiscreteEmission(AlphabetContainer, double[][])
	 */
	protected AbstractConditionalDiscreteEmission( AlphabetContainer con, double[][] hyperParams ) throws WrongAlphabetException {
		this(con,hyperParams,hyperParams);
	}
	
	/**
	 * This constructor creates a {@link AbstractConditionalDiscreteEmission} defining the individual hyper parameters for the
	 * prior used during training and initialization.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param hyperParams the individual hyper parameters for each parameter (used during training)
	 * @param initHyperParams the individual hyper parameters for each parameter used in {@link #initializeFunctionRandomly()}
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	protected AbstractConditionalDiscreteEmission( AlphabetContainer con, double[][] hyperParams, double[][] initHyperParams ) throws WrongAlphabetException {
		if( !con.isDiscrete() || !con.isSimple() ) {
			throw new WrongAlphabetException("The AlphabetContainer has to be discrete and simple.");
		}
		this.con = con;
		ess = new double[hyperParams.length];
		this.hyperParams = new double[hyperParams.length][hyperParams[0].length];
		if(hyperParams == initHyperParams){
			this.initHyperParams = this.hyperParams;
		}else{
			this.initHyperParams = new double[initHyperParams.length][initHyperParams[0].length];
		}
		for( int i = 0; i < hyperParams.length; i++ ) {
			for( int j = 0; j < hyperParams[i].length; j++ ) {
				if( hyperParams[i][j] < 0 ) {
					throw new IllegalArgumentException( "Please check the hyper-parameter (" + i + ", " + j + ")." );
				}
				this.hyperParams[i][j] = hyperParams[i][j];
				if(this.hyperParams != this.initHyperParams){
					this.initHyperParams[i][j] = initHyperParams[i][j];
				}
				ess[i] += hyperParams[i][j];
			}
		}
		params = new double[hyperParams.length][hyperParams[0].length];
		probs = new double[hyperParams.length][hyperParams[0].length];
		statistic = new double[hyperParams.length][hyperParams[0].length];
		grad = new double[hyperParams.length][hyperParams[0].length];
		logNorm = new double[hyperParams.length];
		
		logUniform = - Math.log( con.getAlphabetLengthAt(0) );
		precompute();
	}
	
	/**
	 * Creates a {@link AbstractConditionalDiscreteEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	protected AbstractConditionalDiscreteEmission( StringBuffer xml ) throws NonParsableException {
		fromXML( xml );
	}
	
	public AbstractConditionalDiscreteEmission clone() throws CloneNotSupportedException {
		AbstractConditionalDiscreteEmission clone = (AbstractConditionalDiscreteEmission) super.clone();
		clone.colors = colors == null ? null : colors.clone();
		if(counter != null){
			clone.counter = counter.clone();
		}
		if(ess != null){
			clone.ess = ess.clone();
		}
		if(grad != null){
			clone.grad = ArrayHandler.clone( grad.clone() );
		}
		clone.hyperParams = ArrayHandler.clone( hyperParams );
		clone.initHyperParams = ArrayHandler.clone( initHyperParams );
		if(logNorm != null){
			clone.logNorm = logNorm.clone();
		}
		clone.params = ArrayHandler.clone( params );
		if( paramsFile != null ) {
			clone.paramsFile = new File[paramsFile.length];
			try {
				for( int i = 0; i < paramsFile.length; i++ ) {
					if( paramsFile[i] != null ) {
						clone.paramsFile[i] = createFile();
						FileManager.copy( paramsFile[i].getAbsolutePath(), clone.paramsFile[i].getAbsolutePath() );
					}
				}
			} catch( IOException e ) {
				CloneNotSupportedException c = new CloneNotSupportedException( e.getMessage() );
				c.setStackTrace( e.getStackTrace() );
				throw c;
			}
		}
		clone.probs = ArrayHandler.clone( probs );
		clone.reader = null; //XXX ??
		clone.statistic = ArrayHandler.clone( statistic );
		clone.writer = null; //XXX ??
		return clone;
	}
	
	/**
	 * Sets the graphviz shape of the node that uses this emission to some non-standard value
	 * (standard is &quot;house&quot;).
	 * @param shape the shape of the node
	 */
	public void setShape(String shape){
		this.shape = shape;
	}
	
	public void addGradientOfLogPriorTerm(double[] gradient, int offset) {
		for( int i = 0; i < params.length; i++ ) {
			for( int j = 0; j < params[i].length; j++, offset++ ) {
				gradient[offset+this.offset] += hyperParams[i][j] - ess[i]*probs[i][j];
			}
		}
	}

	public double getLogPriorTerm() {
		double res = 0;
		for( int i = 0; i < params.length; i++ ) {
			if( ess[i] > 0 ) {
				res += -ess[i]*logNorm[i];
				for( int j = 0; j < params[i].length; j++ ) {
					res += hyperParams[i][j] * params[i][j];
				}
			}
		}
		return res;
	}

	public double getLogProbAndPartialDerivationFor( boolean forward, int startPos, int endPos,
			IntList indices, DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
		int s, e;
		Sequence current;
		if( forward ) {
			current = seq;
			s = startPos;
			e = endPos;
		} else {
			current = seq.reverseComplement();
			int len = current.getLength();
			s = len - endPos -1;
			e = len - startPos -1;
		}
		int v = e-s+1;
		
		double res = 0;
		for(int i=0;i<grad.length;i++){
			Arrays.fill(grad[i],0);
		}
		while( s <= e ) {
			int condIdx = getConditionIndex( forward, s, seq );
			if(condIdx < 0){
				//condition value not available => uniform & no gradient
				res += logUniform;
			} else {
				v = getIndex(s, current);
				res -= logNorm[condIdx];
				for(int i=0;i<grad[condIdx].length;i++) {
					grad[condIdx][i] -= probs[condIdx][i];
				}
				res += params[condIdx][v];
				grad[condIdx][v]++;
			}
			s++;
		}
		int myOff = 0;
		for( int i = 0; i< grad.length; i++ ) {
			for(int j=0;j<grad[i].length;j++, myOff++){
				if( grad[i][j]!=0 ) {
					indices.add( offset + myOff );
					partDer.add( grad[i][j] );
				}
			}
		}
		return res;
	}

	public double getLogProbFor( boolean forward, int startPos, int endPos, Sequence seq) throws OperationNotSupportedException {
		int s, e;
		Sequence current;
		if( forward ) {
			current = seq;
			s = startPos;
			e = endPos;
		} else {
			current = seq.reverseComplement();
			int len = current.getLength();
			s = len - endPos -1;
			e = len - startPos -1;
		}
		
		double res = 0;
		while( s <= e ) {
			int condIdx = getConditionIndex( forward, s, seq );
			if(condIdx < 0){
				res += logUniform;
			} else {
				res -= logNorm[condIdx];
				res += params[condIdx][getIndex(s, current)];
			}
			s++;
		}
		return res;
	}

	public void initializeFunctionRandomly() {
		drawParameters( initHyperParams, true );
	}

	/**
	 * This method precomputes some normalization constant and probabilities.
	 * 
	 * @see #logNorm
	 * @see #probs
	 */
	protected void precompute() {
		for(int i = 0 ; i < params.length; i++ ) {
			logNorm[i] = Normalisation.logSumNormalisation(params[i], 0, params[i].length, probs[i], 0 );
			/*logNorm[i] = Normalisation.getLogSum( params[i] );
			for( int j = 0; j < params[i].length; j++ ) {
				probs[i][j] = Math.exp( params[i][j] - logNorm[i] );
			}*/
		}
	}

	private static final String XML_TAG = "ConditionalDiscreteEmission"; 
	
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, con, "alphabetContainer" );
		XMLParser.appendObjectWithTags( xml, params, "params" );
		XMLParser.appendObjectWithTags( xml, offset, "offset" );
		XMLParser.appendObjectWithTags( xml, hyperParams, "hyperParams" );
		XMLParser.appendObjectWithTags( xml, initHyperParams, "initHyperParams" );
		XMLParser.appendObjectWithTags( xml, statistic, "statistic" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.appendObjectWithTags( xml, shape, "shape" );
		XMLParser.appendObjectWithTags( xml, linear, "linear" );
		
		if( writer != null ) {
			throw new RuntimeException( "could not parse SamplingHigherOrderTransition to XML while sampling" );
		}
		
		XMLParser.appendObjectWithTags( xml, paramsFile != null, "hasParameters" );
		if( paramsFile != null ) {
			String content;
			try {
				XMLParser.appendObjectWithTags( xml, counter, "counter" );
				
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
		
		appendFurtherInformation( xml );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	/**
	 * This method appends further information to the XML representation. It allows subclasses to save further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 */
	protected void appendFurtherInformation( StringBuffer xml ) {
	}

	/**
	 * This method is internally used by the constructor {@link #AbstractConditionalDiscreteEmission(StringBuffer)}.
	 * 
	 * @param xml the {@link StringBuffer} containing the xml representation of an instance
	 * 
	 * @throws NonParsableException if the {@link StringBuffer} is not parsable
	 * 
	 * @see #AbstractConditionalDiscreteEmission(StringBuffer)
	 */
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		con = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "alphabetContainer" );
		logUniform = - Math.log( con.getAlphabetLengthAt(0) );
		params = (double[][]) XMLParser.extractObjectForTags( xml, "params" );
		probs = new double[params.length][params[0].length];
		grad = new double[params.length][params[0].length];
		logNorm = new double[params.length];
		precompute();
		
		offset = (Integer) XMLParser.extractObjectForTags( xml, "offset" );
		hyperParams = (double[][]) XMLParser.extractObjectForTags( xml, "hyperParams" );
		try{
			initHyperParams = (double[][]) XMLParser.extractObjectForTags( xml, "initHyperParams" );
		}catch(NonParsableException e){
			try{
				initHyperParams = ArrayHandler.clone( hyperParams );
			}catch(CloneNotSupportedException ex){}
		}
		statistic = (double[][]) XMLParser.extractObjectForTags( xml, "statistic" );
		ess = (double[]) XMLParser.extractObjectForTags( xml, "ess" );
		try{
			shape = XMLParser.extractObjectForTags( xml, "shape", String.class );
			linear = XMLParser.extractObjectForTags( xml, "linear",boolean.class );
		}catch(NonParsableException e){
			shape = null;
			linear = false;
		}
		
		if( XMLParser.extractObjectForTags( xml, "hasParameters", boolean.class ) ) {
			counter = XMLParser.extractObjectForTags( xml, "counter", int[].class );
			paramsFile = new File[counter.length];
			try {
				String content;
				Map<String,String> filter = new TreeMap<String, String>();
				for( int i = 0; i < paramsFile.length; i++ ) {
					filter.clear();
					filter.put( "pos", ""+i );
					content = XMLParser.extractObjectAndAttributesForTags( xml, "fileContent", null, filter, String.class );
					if( !content.equalsIgnoreCase( "" ) ) {
						paramsFile[i] = createFile();
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
		
		extractFurtherInformation( xml );
	}
	
	/**
	 * This method extracts further information from the XML representation. It allows subclasses to cast further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 *  
	 * @throws NonParsableException if the information could not be reconstructed out of the {@link StringBuffer} <code>xml</code>
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
	}

	public void joinStatistics(Emission... emissions){
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				for(int j=0;j<statistic.length;j++){
					for(int k=0;k<statistic[j].length;k++){
						statistic[j][k] += ((AbstractConditionalDiscreteEmission)emissions[i]).statistic[j][k];
					}
				}
			}
		}
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				for(int j=0;j<statistic.length;j++){
					System.arraycopy( this.statistic[j], 0, ((AbstractConditionalDiscreteEmission)emissions[i]).statistic[j], 0, this.statistic[j].length );
				}
			}
		}
	}
	
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {
		int s, e;
		Sequence current;
		
		if( forward ) {
			current = seq;
			s = startPos;
			e = endPos;
		} else {
			current = seq.reverseComplement();
			int len = current.getLength();
			s = len - endPos -1;
			e = len - startPos -1;
		}
		
		while( s <= e ) {
			int condIdx = getConditionIndex( forward, s, seq); //TODO error? seq vs. current
			if( condIdx>= 0 ) {
				statistic[condIdx][getIndex(s, current)] += weight;
			}
			s++;
		}
	}

	/**
	 * This method returns an index encoding the condition.
	 * 
	 * @param forward a switch to decide whether to use the forward or the reverse complementary strand (e.g. for DNA sequences)
	 * @param seqPos the position in the sequence <code>seq</code>
	 * @param seq the sequence
	 * 
	 * @return the index encoding the condition
	 */
	protected abstract int getConditionIndex( boolean forward, int seqPos, Sequence seq );
	
	/**
	 * Returns the index for position <code>seqPos</code> in sequence <code>seq</code>.
	 * @param seqPos the position
	 * @param seq the sequence
	 * @return the index
	 */
	protected int getIndex( int seqPos, Sequence seq ) {
		return seq.discreteVal( seqPos );
	}

	public void estimateFromStatistic() {
		for(int j=0;j<statistic.length;j++){
			double sum = 0;
			for( int i = 0; i < statistic[j].length; i++ ) {
				statistic[j][i] += hyperParams[j][i];
				sum += statistic[j][i];
			}
			if( sum == 0 ) {
				Arrays.fill( statistic[j], 1 );
				sum = statistic[j].length;
			}
			for( int i = 0; i < statistic[j].length; i++ ) {
				probs[j][i] = statistic[j][i] / sum;
				params[j][i] = Math.log( probs[j][i] );
			}
		}
		Arrays.fill( logNorm, 0 );
	}

	public void resetStatistic() {
		for(int i=0;i<hyperParams.length;i++){
			Arrays.fill( statistic[i], 0 );
		}
	}
	
	public void setParameter( double[] params, int offset ) {
		for( int i = 0; i < this.params.length; i++ ) {
			for( int j = 0; j < this.params[i].length; j++, offset++) {
				this.params[i][j] = params[this.offset+offset];
			}
		}
		precompute();
	}
	
	public AlphabetContainer getAlphabetContainer() {
		return con;
	}
	
	public void fillCurrentParameter( double[] params ) {
		int myOffset = offset;
		for( int i = 0; i < this.params.length; i++ ) {
			for(int j = 0; j < this.params[i].length; j++, myOffset++ ) {
				params[myOffset] = this.params[i][j];
			}
		}
	}
	
	public int setParameterOffset( int offset ) {
		this.offset = offset;
		for(int i=0;i<params.length;i++){
			offset += params[i].length;
		}
		return offset;
	}
	
	/**
	 * Draws the parameters of this {@link AbstractConditionalDiscreteEmission} from a Dirichlet distribution
	 * with given hyper-parameters. If the equivalent sample size (ess) according to the provided hyper-parameters
	 * is zero, parameters may be drawn from a uniform distribution on the simplex.
	 * @param hyper the hyper-parameters
	 * @param uniformBackup if a uniform distribution should be used in case of ess zero
	 */
	protected void drawParameters( double[][] hyper, boolean uniformBackup ) {
		double ess;
		DiMRGParams p;
    	for(int j=0;j<probs.length;j++){
    		ess = 0;
    		if( uniformBackup ) {
				for(int i=0;i<hyper[j].length;i++){
					ess += hyper[j][i];
				}
    		}
			if( uniformBackup && ess == 0 ) {
				p = new FastDirichletMRGParams(1d);
			}else{
				p = new DirichletMRGParams( hyper[j] );
			}	
    		DirichletMRG.DEFAULT_INSTANCE.generateLog( params[j], 0, hyper[j].length, p );
    	}
    	precompute();
	}
	
    public void drawParametersFromStatistic() {
    	for(int j=0;j<probs.length;j++){
    		for( int i = 0; i < statistic[j].length; i++ ) {
    			params[i][j] = statistic[i][j] + hyperParams[i][j];
    		}
		}
    	drawParameters( params, false );
	}

    public double getLogGammaScoreFromStatistic() {

    	double[][] hyper = getHyperParams( ess, (int) con.getAlphabetLengthAt( 0 ) );
        double res = Double.NEGATIVE_INFINITY;
        for(int j=0;j<ess.length;j++){
        	
        	double sum = 0;

        	for (double i : hyper[j]) sum += i;

        	res = Gamma.logOfGamma(sum);
        	for(int i = 0; i < hyper.length; i++)
        		res += Gamma.logOfGamma(statistic[j][i]) - Gamma.logOfGamma(hyper[j][i]);

        	sum = 0;
        	for (double i : statistic[j]) sum += i;

        	res -= Gamma.logOfGamma(sum);
        }
        return res;
    }
	
	
	public void acceptParameters() throws IOException {
		writer.write( "" + ( counter[samplingIndex]++ ) );
		for( int i = 0; i < params.length; i++ ) {
			for(int j=0;j<params[i].length;j++){
				writer.write( "\t" + params[i][j] );
			}
		}
		writer.write( "\n" );
		writer.flush();
	}
	
	public double getLogPosteriorFromStatistic() {
		double logPost =0;
		 for( int i = 0; i < params.length; i++ ) {
			 for(int j=0;j<params[i].length;j++){
				 logPost += statistic[i][j] * (params[i][j]-logNorm[i]);
			 }
         }
		return logPost;
	}

	public void extendSampling(int start, boolean append) throws IOException {
		
		if( paramsFile[start] == null ) {
			paramsFile[start] = createFile();
		} else {
			if( append ) {
				parseParameterSet( start, counter[start] - 1 );
				reader.close();
				reader = null;
			} else {
				counter[start] = 0;
			}
		}
		writer = new BufferedWriter( new FileWriter( paramsFile[start], append ) );
		samplingIndex = start;
	}

	public void initForSampling( int starts ) throws IOException {

		for(int i = 0; i < hyperParams.length; i++) {
			for(int j=0;j<hyperParams[i].length;j++){
				if(!Double.isNaN(hyperParams[i][j]) && hyperParams[i][j] <= 0){
					throw new IllegalArgumentException( "All (not NAN) hyper-parameters must have a value > 0. Please check the hyper-parameter " + i + "." );
				}
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

	public boolean isInSamplingMode() {
		return writer != null;
	}

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

	public boolean parseParameterSet(int start, int n) throws IOException {
		String str;
		if( reader != null ) {
			reader.close();
		}
		reader = new BufferedReader( new FileReader( paramsFile[start] ) );
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
		
		for( int i=0; i < params.length; i++ ) {
			for(int j=0;j<params[i].length;j++){
				params[i][j] = Double.parseDouble(strArray[offset++]);
				probs[i][j] = Math.exp(params[i][j]);
			}
			logNorm[i] = 0;
		}
		
	}

	public void samplingStopped() throws IOException {
		
		if( writer != null ) {
			writer.close();
			writer = null;
		}
	}
	
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
	public String getNodeShape(boolean forward) {
		String res;
		if( shape == null ) {
			res = "";
			if( getAlphabetContainer().isReverseComplementable() ) {
				res += "\"house\", orientation=";
				if( forward ) {
					res += "-";
				}
				res += "90";
			} else {
				res += "\"box\"";
			}
		} else {
			res = "\""+shape+"\"";
		}
		return res;
	}

	@Override
	public String getNodeLabel( double weight, String name, NumberFormat nf ) {
		if(weight < 0){
			return "\""+name+"\"";
		}else{
			StringBuffer buf = new StringBuffer();

			String namelabel = name;
			if(weight < 0.5){
				namelabel = "<font color=\"white\">"+namelabel+"</font>";
			}
			buf.append( "<<table border=\"0\" cellspacing=\"0\"><tr><td colspan=\""+((probs.length>1?1:0) + probs[0].length)+"\">"+namelabel+"</td></tr>" );
			
			DiscreteAlphabet abc = (DiscreteAlphabet) con.getAlphabetAt(0);
			buf.append( "<tr>" );
			if(probs.length>1){
				buf.append( "<td></td>" );
			}
			for(int j=0;j<probs[0].length;j++){
				buf.append( "<td border=\"0\">" );
				if(weight < 0.5){
					buf.append( "<font color=\"white\">"+abc.getSymbolAt(j)+"</font>" );
				} else {
					buf.append( abc.getSymbolAt(j) );
				}
				buf.append( "</td>" );
			}
			buf.append( "</tr>" );


			for(int i=0;i<probs.length;i++){
				buf.append( "<tr>" );
				if( probs.length > 1) {
					buf.append( "<td border=\"0\">" );
					if(weight < 0.5){
						buf.append( "<font color=\"white\">"+abc.getSymbolAt(i)+"</font>" );
					} else {
						buf.append( abc.getSymbolAt(i) );
					}
					buf.append( "</td>" );
				}
				double[] trans =  transformProbs( probs[i] );
				double en = (getInformationContent( probs[i] )+2.0)/3.0;
				for(int j=0;j<probs[i].length;j++){
					buf.append( "<td border=\"1\" width=\"25\" height=\"25\" bgcolor=\""+getColor(j)+" "+trans[j]+" "+en+"\">" );
					if( nf != null ) {
						buf.append( nf.format( probs[i][j]) );
					}
					buf.append( "</td>" );
				}
				buf.append( "</tr>" );
			}
			buf.append( "</table>>" );
			return buf.toString();
		}
		
	}
	
	private double getColor( int j ) {
		if(colors == null){
			colors = ToolBox.getUniqueHueValues( (int) con.getAlphabetLengthAt( 0 ) );
		}
		return colors[j];
	}

	/**
	 * If set to true, the probabilities are mapped to colors by directly, otherwise
	 * a logistic mapping is used to emphasize deviations from the uniform distribution.
	 * @param linear map probabilities linear
	 */
	public void setLinear( boolean linear ) {
		this.linear = linear;
	}
	
	private double getInformationContent(double[] probs){
		double max = Math.log( probs.length );
		double en = 0.0;
		for(int i=0;i<probs.length;i++){
			if(probs[i]>0){
				en -= probs[i] * Math.log( probs[i] );
			}
		}
		return (max-en)/max;
	}
	
	private double[] transformProbs(double[] probs){
		if(linear){
			return probs.clone();
		}else{
			double[] trans = new double[probs.length];

			double a=15,b=4;

			for(int i=0;i<probs.length;i++){
				trans[i] = 1.0/(1.0+Math.exp( -a*probs[i] + b ));
			}
			return trans;
		}
	}

	@Override
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list ) {
		int off = 0;
		for(int i=0;i<params.length;i++){
			int[] idxs = new int[params[i].length];
			for(int j=0;j<idxs.length;j++){
				idxs[j] = j + off + offset + parameterOffset;
			}
			list.add( idxs );
			off += idxs.length;
		}
	}

	@Override
	public int getNumberOfParameters() {
		return params.length*params[0].length;
	}

	@Override
	public int getSizeOfEventSpace() {
		return params.length*params[0].length;
	}
	
	private File createFile() throws IOException {
		return File.createTempFile( "samplingDEmission-", ".dat", null );
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) || ((AbstractConditionalDiscreteEmission)t).params.length != params.length ) {//TODO more?
			throw new IllegalArgumentException( "The transitions are not comparable." );
		}
		AbstractConditionalDiscreteEmission tt = (AbstractConditionalDiscreteEmission) t;
		for( int i = 0; i < params.length; i++ ) {
			System.arraycopy( tt.params[i], 0, params[i], 0, tt.params[i].length );
		}
		precompute();
	}
	
	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer();
		for( int i = 0; i < probs.length; i++ ) {
			for( int j = 0; j < probs[i].length; j++ ) {
				String cond;
				if( probs.length==1 ) {
					cond = "";
				} else {
					cond= getCondition(i);
				}
				sb.append( "P(X=" + con.getSymbol( 0, j ) + cond +") = " + nf.format( probs[i][j] ) + "\t");
				
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	/**
	 * Return the String representation of condition <code>i</code>.
	 * @param i the index of the condition
	 * 
	 * @return the String representation of condition <code>i</code>
	 */
	protected abstract String getCondition( int i );
}