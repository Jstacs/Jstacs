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

package de.jstacs.classifiers.differentiableSequenceScoreBased;

import java.io.OutputStream;
import java.util.Arrays;

import de.jstacs.NotTrainedException;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.StartDistanceForecaster;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.ClassDimensionException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox;
import de.jstacs.motifDiscovery.history.History;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.SafeOutputStream;

/**
 * This abstract class implements the main functionality of a {@link DifferentiableSequenceScore} based classifier.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public abstract class ScoreClassifier extends AbstractScoreBasedClassifier {

	/**
	 * The internally used scoring functions.
	 */
	protected DifferentiableSequenceScore[] score;

	/**
	 * The parameter set for the classifier.
	 */
	protected ScoreClassifierParameterSet params;

	/**
	 * This boolean indicates whether the classifier has been optimized with the
	 * method {@link de.jstacs.classifiers.AbstractClassifier#train(DataSet[])} or
	 * the weighted version.
	 */
	protected boolean hasBeenOptimized;

	/**
	 * This value is the score of the last optimization. 
	 */
	private double lastScore;

	/**
	 * This stream is used for comments, e.g. during the training, ... .
	 */
	protected SafeOutputStream sostream;

	/**
	 * This value should be used in {@link ScoreClassifier#getLastScore()} if the classifier is not trained. 
	 */
	public static final double NOT_TRAINED_VALUE = Double.NaN;

	/**
	 * The default history
	 */
	protected History template = null;//TODO new NoRevertHistory();
	
	/**
	 * Creates a new {@link ScoreClassifier} from a given
	 * {@link ScoreClassifierParameterSet} and {@link DifferentiableSequenceScore}s .
	 * 
	 * @param params
	 *            the parameter set for the classifier
	 * @param lastScore
	 *            the score of the last optimization, if no such value exists the
	 *            programmer should use {@link ScoreClassifier#NOT_TRAINED_VALUE}
	 * @param score
	 *            the {@link DifferentiableSequenceScore}s for the classes
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one {@link DifferentiableSequenceScore} could not be cloned
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(AlphabetContainer,
	 *      int, int)
	 */
	public ScoreClassifier( ScoreClassifierParameterSet params, double lastScore, DifferentiableSequenceScore... score ) throws CloneNotSupportedException {
		super( params.getAlphabetContainer(), params.getLength(), score.length );
		int i = 0, l, len = getLength();
		AlphabetContainer con = getAlphabetContainer();
		while( i < score.length ) {
			l = score[i].getLength();
			if( ( l == 0 || l == len ) && con.checkConsistency( score[i].getAlphabetContainer() ) ) {
				//everything is okay
				i++;
			} else {
				throw new IllegalArgumentException( "Please check the length (" + l + " vs. " + len + ")" +
						" and the AlphabetContainer of the DifferentiableSequenceScore with index " + i + "." );
			}
		}
		this.score = ArrayHandler.clone( score );
		hasBeenOptimized = false;
		if( isInitialized() ) {
			this.lastScore = lastScore;
		} else {
			this.lastScore = NOT_TRAINED_VALUE;
		}
		set( (ScoreClassifierParameterSet)params.clone() );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ScoreClassifier} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ScoreClassifier} could not be reconstructed out
	 *             of the XML representation (the {@link StringBuffer} could not
	 *             be parsed)
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public ScoreClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#clone()
	 */
	@Override
	public ScoreClassifier clone() throws CloneNotSupportedException {
		ScoreClassifier clone = (ScoreClassifier)super.clone();
		clone.params = (ScoreClassifierParameterSet)params.clone();
		clone.score = ArrayHandler.clone( score );
		clone.setOutputStream( this.sostream.doesNothing() ? null : SafeOutputStream.DEFAULT_STREAM );
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getClassifierAnnotation()
	 */
	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		CategoricalResult[] res = new CategoricalResult[score.length + 1];
		res[0] = new CategoricalResult( "classifier", "a <b>short</b> description of the classifier", getInstanceName() );
		int i = 0;
		while( i < score.length ) {
			res[i + 1] = new CategoricalResult( "class info " + i, "some information about the class", score[i++].getInstanceName() );
		}
		return res;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumericalCharacteristics()
	 */
	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {

		NumericalResult[] pars = new NumericalResult[score.length + ( hasBeenOptimized ? 1 : 0 )];
		if( hasBeenOptimized ) {
			pars[0] = new NumericalResult( "Last score", "The final score after the optimization", lastScore );
		}
		for( int i = 0; i < score.length; i++ ) {
			pars[i + ( hasBeenOptimized ? 1 : 0 )] = new NumericalResult( "Number of parameters " + ( i + 1 ),
					"The number of parameters for scoring function " + ( i + 1 ) + ", -1 indicates unknown number of parameters.",
					score[i].getNumberOfParameters() );
		}
		return new NumericalResultSet( pars );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		int i = 0;
		while( i < score.length && score[i].isInitialized() ) {
			i++;
		}
		return i == score.length;
	}

	/**
	 * This method indicates if the classifier has been optimized by a
	 * <code>train</code>-method.
	 * 
	 * @return <code>true</code> if the classifier has been optimized by a
	 *         <code>train</code>-method, <code>false</code> otherwise
	 */
	public boolean hasBeenOptimized() {
		return hasBeenOptimized;
	}

	/**
	 * Sets the {@link OutputStream} that is used e.g. for writing information
	 * during training. It is possible to set <code>o=null</code>, then nothing
	 * will be written.
	 * 
	 * @param o
	 *            the {@link OutputStream}
	 */
	public void setOutputStream( OutputStream o ) {
		sostream = SafeOutputStream.getSafeOutputStream( o );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#train(de.jstacs.data.DataSet[], double[][])
	 */
	@Override
	public void train( DataSet[] data, double[][] weights ) throws Exception {
		hasBeenOptimized = false;
		// check
		if( weights == null ) {
			weights = new double[data.length][];
		}
		if( data.length > 1  && data.length != weights.length ) {
			throw new IllegalArgumentException( "data and weights do not match" );
		}
		if( score.length != weights.length ) {
			throw new ClassDimensionException();
		}
		
		WeightedDataSetFactory wsf;
		DataSet[] reduced = new DataSet[data.length];
		double[][] newWeights = new double[weights.length][];
		AlphabetContainer abc = getAlphabetContainer();
		int j = 0;
		for( int l = getLength(), i = 0; i < score.length; i++ ) {
			if( weights[i] != null && data[j].getNumberOfElements() != weights[i].length ) {
				throw new IllegalArgumentException( "At least for one data set: The dimension of the data set and the weight do not match." );
			}
			if( i == 0 || data.length > 1 ) {
				if( !abc.checkConsistency( data[j].getAlphabetContainer() ) ) {
					throw new IllegalArgumentException( "At least one data set is not defined over the correct alphabets." );
				}
				if( data[i].getElementLength() != l ) {
					// throw new IllegalArgumentException( "At least one data set has not the correct length." );
					wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[i], weights[i], l );
				} else {
					wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[i], weights[i] );
				}
				reduced[i] = wsf.getDataSet();
				newWeights[i] = wsf.getWeights();
			} else {
				
				if( data[j].getElementLength() != l ) {
					// throw new IllegalArgumentException( "At least one data set has not the correct length." );
					wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[j], weights[i], l );
				} else {
					wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[j], weights[i] );
				}
				
				newWeights[i] = wsf.getWeights();
			}
			
			if( data.length > 1 ) {
				j++;
			}
		}
		lastScore = doOptimization( reduced, newWeights );
	}

	/**
	 * This method does the optimization of the <code>train</code>-method
	 * 
	 * @param reduced
	 *            the data sets
	 * @param newWeights
	 *            the weights
	 * 
	 * @return the value of the optimization
	 * 
	 * @throws Exception
	 *             if something went wrong during the optimization
	 */
	protected double doOptimization( DataSet[] reduced, double[][] newWeights ) throws Exception {
		// train
		byte algo = (Byte) params.getParameterForName( "algorithm" ).getValue();
		AbstractTerminationCondition tc = params.getTerminantionCondition();
		double linEps = (Double) params.getParameterForName( "line epsilon" ).getValue(),
			startDist = (Double) params.getParameterForName( "start distance" ).getValue();
		KindOfParameter plugIn;

		double[] best = null;
		double[][] res = new double[2][];		
		double max = Double.NEGATIVE_INFINITY, current;
		
		int iterations = AbstractDifferentiableStatisticalModel.getNumberOfStarts( score );
		sostream.writeln( getInstanceName() );
		DifferentiableSequenceScore[] bestSF = new DifferentiableSequenceScore[score.length], secure;
		if( iterations > 1 ) {
			secure = ArrayHandler.clone( score );
		} else {
			secure = null;
		}
		History[][] hist = MutableMotifDiscovererToolbox.createHistoryArray( score, template );
		int[][] minimalNewLength = MutableMotifDiscovererToolbox.createMinimalNewLengthArray( score );
		StartDistanceForecaster sd = new ConstantStartDistance( startDist );
		
		DiffSSBasedOptimizableFunction f = getFunction( reduced, newWeights );
		//old DifferentiableFunction g = new NegativeDifferentiableFunction( f );
		for( int i = 0; i < iterations; ) {
			// create structure
			createStructure( reduced, newWeights );
			f.reset( score );
			
			if( i == 0 ) {
				sostream.writeln( "optimizing " + f.getDimensionOfScope() + " parameters" );
			}
			sostream.writeln( "start " + ++i + ":" );

			//pre-optimize
			plugIn = preoptimize(f);	
			
			//new
			MutableMotifDiscovererToolbox.clearHistoryArray( hist );
			sd.reset();
			res = MutableMotifDiscovererToolbox.optimize( score, f, algo, tc, linEps, sd, sostream, false, hist, minimalNewLength, plugIn, true );//TODO
			current = res[0][0];

			/*
			//old
			res[1]= f.getParameters( plugIn );
			Optimizer.optimize( algo,
					g,
					res[1],
					Optimizer.TerminationCondition.SMALL_DIFFERENCE_OF_FUNCTION_EVALUATIONS,
					eps,
					linEps,
					new ConstantStartDistance( startDist ),
					sostream );
			current = f.evaluateFunction( res[1] );
			*/
			if( current > max ) {
				System.arraycopy( score, 0, bestSF, 0, score.length );
				//bestSF = score;
				best = res[1];
				max = current;
				System.gc();
			}
			if( iterations > 1 ) {
				score = ArrayHandler.clone( secure );
				if( sostream.doesNothing() ) {
					//System.out.println( "start " + i + ": " + res[0][0] );
				}
			}
		}
		sostream.writeln( "best = " + max );

		score = bestSF;
		setClassWeights( false, best );
		hasBeenOptimized = true;
		
		if( f instanceof AbstractMultiThreadedOptimizableFunction ) {
			((AbstractMultiThreadedOptimizableFunction)f).stopThreads();
		}
		return max;
	}
	
	/**
	 * This method allows to pre-optimize the parameter before the real optimization. It might be used useful, for instance,
	 * to find initial parameters on a part of the complete data.
	 * 
	 * @param f the function to be optimized
	 * 
	 * @return enumeration value that indicates how to retrieve the parameters from the function
	 * 
	 * @throws Exception if the pre-optimization fails
	 */
	protected KindOfParameter preoptimize( OptimizableFunction f ) throws Exception {
		return (KindOfParameter) params.getParameterForName( KindOfParameter.class.getSimpleName() ).getValue();
	}	
	
	/**
	 * Creates the structure that will be used in the optimization.
	 * 
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights of the data
	 * @param initRandomly
	 * 			  initialize the functions randomly
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected void createStructure( DataSet[] data, double[][] weights, boolean initRandomly ) throws Exception {
		boolean freeParams = params.useOnlyFreeParameter();
		DataSet[] d;
		if( !initRandomly && data.length == 1 && weights != null && weights.length > 1 ) {
			d = new DataSet[weights.length];
			Arrays.fill( d, data[0] );
		} else {
			d = data;
		}
		for( int i = 0; i < score.length; i++ ) {
			if(initRandomly){
				score[i].initializeFunctionRandomly( freeParams );
			}else{
				score[i].initializeFunction( i, freeParams, d, weights );
			}
		}
	}
	
	/**
	 * Sets the parameters of this classifier and the contained scoring functions
	 * to the supplied parameters.
	 * @param parameters the new parameters
	 * @throws Exception if the parameters could not be set
	 */
	public void initUsingParameters(double[] parameters) throws Exception{
		createStructure( null, null, true );
		double[] cw = new double[score.length];
		for(int i=0;i<cw.length - (params.useOnlyFreeParameter() ? 1 : 0);i++){
			cw[i] = parameters[i];
		}
		setClassWeights( false, cw, 0 );
		int off = score.length - (params.useOnlyFreeParameter() ? 1 : 0);
		for(int i=0;i<score.length;i++){
			score[i].setParameters( parameters, off );
			off += score[i].getNumberOfParameters();
		}
	}
	
	/**
	 * Creates the structure that will be used in the optimization.
	 * 
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights of the data	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected void createStructure( DataSet[] data, double[][] weights) throws Exception{
		createStructure( data, weights, false );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherClassifierInfosFromXML( xml );
		set( (ScoreClassifierParameterSet) XMLParser.extractObjectForTags( xml, "params" ) );
		hasBeenOptimized = XMLParser.extractObjectForTags( xml, "hasBeenOptimized", boolean.class );
		lastScore = XMLParser.extractObjectForTags( xml, "lastScore", double.class );
		score = XMLParser.extractObjectForTags( xml, "score", DifferentiableSequenceScore[].class );
	}
	
	/**
	 * Returns the function that should be optimized.
	 * 
	 * @param data
	 *            the data sets
	 * @param weights
	 *            the weights of the sequences of the data sets
	 * 
	 * @return the function that should be optimized
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected abstract DiffSSBasedOptimizableFunction getFunction( DataSet[] data, double[][] weights ) throws Exception;

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, params, "params" );
		XMLParser.appendObjectWithTags( xml, hasBeenOptimized, "hasBeenOptimized" );
		XMLParser.appendObjectWithTags( xml, lastScore, "lastScore" );
		XMLParser.appendObjectWithTags( xml, score, "score" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getScore(de.jstacs.data.Sequence, int, boolean)
	 */
	@Override
	protected double getScore( Sequence seq, int i, boolean check ) throws IllegalArgumentException, NotTrainedException, Exception {
		if( check ) {
			check( seq );
		}
		return getClassWeight( i ) + score[i].getLogScoreFor( seq, 0 );
	}

	/**
	 * Returns the score that was computed in the last optimization of the
	 * parameters.
	 * 
	 * @return score from the last parameter optimization
	 */
	public double getLastScore() {
		return lastScore;
	}

	/**
	 * Returns the internally used {@link DifferentiableSequenceScore} with index
	 * <code>i</code>.
	 * 
	 * @param i
	 *            the internal index of the {@link DifferentiableSequenceScore}
	 * 
	 * @return the internally used {@link DifferentiableSequenceScore} with index
	 *         <code>i</code>
	 * 
	 * @throws CloneNotSupportedException
	 *             if the {@link DifferentiableSequenceScore} could not be cloned
	 */
	public DifferentiableSequenceScore getDifferentiableSequenceScore( int i ) throws CloneNotSupportedException {
		return score[i].clone();
	}

	/**
	 * Returns all internally used {@link DifferentiableSequenceScore}s in the internal
	 * order.
	 * 
	 * @return the internally used {@link DifferentiableSequenceScore}s in the internal
	 *         order
	 * 
	 * @throws CloneNotSupportedException
	 *             if a {@link DifferentiableSequenceScore} could not be cloned
	 */
	public DifferentiableSequenceScore[] getDifferentiableSequenceScores() throws CloneNotSupportedException {
		return ArrayHandler.clone( score );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getXMLTag()
	 */
	@Override
	protected abstract String getXMLTag();

	private void set( ScoreClassifierParameterSet params ) {
		this.params = params;
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
	}
	
	/**
	 * This method returns the current {@link de.jstacs.parameters.ParameterSet} of the classifier.
	 * @return the current {@link de.jstacs.parameters.ParameterSet} of the classifier.
	 * @throws CloneNotSupportedException if the {@link de.jstacs.parameters.ParameterSet} could not be cloned
	 */
	public ScoreClassifierParameterSet getCurrentParameterSet() throws CloneNotSupportedException {
		return (ScoreClassifierParameterSet) params.clone();
	}
}
