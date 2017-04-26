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

package de.jstacs.classifiers.differentiableSequenceScoreBased.sampling;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NotTrainedException;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.Function;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.ClassDimensionException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sampling.SamplingComponent;
import de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel;
import de.jstacs.utils.Pair;

/**
 * A classifier that samples the parameters of {@link SamplingDifferentiableStatisticalModel}s by the Metropolis-Hastings algorithm.
 * The distribution the parameters are sampled from is the distribution {@latex.inline $P(\\vec{\\lambda}^{t})$} represented by the {@link DiffSSBasedOptimizableFunction} returned by
 * {@link SamplingScoreBasedClassifier#getFunction(DataSet[], double[][])}. As proposal distribution, a Gaussian distribution with given sampling
 * variance is used for each parameter.
 * Specifically, a new set of parameters {@latex.inline $\\vec{\\lambda}^{t}$} is drawn from a proposal distribution {@latex.inline $Q(\\vec{\\lambda}^{t} | \\vec{\\lambda}^{t-1})$},
 * where
 * {@latex.ilb \\[ Q(\\vec{\\lambda}^{t}|\\vec{\\lambda}^{t-1}) = \\prod_{i} \\mathcal{N}(\\lambda_i^{t}|\\lambda_i^{t-1},\\sigma_i^2)\\]}
 * and {@latex.inline $\\sigma_i^2$} is the sampling variance for parameter {@latex.inline $\\lambda_i$}. The sampling variances are adapted to the
 * size of the event space of each parameter based on a class-dependent variance provided to the constructor. This adaption depends on the correct
 * implementation of {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getSizeOfEventSpaceForRandomVariablesOfParameter(int)}. Let {@latex.inline $s_i$} be the size of the event space
 * of the random variable of parameter {@latex.inline $\\lambda_i$}, and let {@latex.inline $\\sigma^{2}$} be the class-dependent variance for the {@link SamplingDifferentiableStatisticalModel}
 * that {@latex.inline $\\lambda_i$} is a parameter of. Then {@latex.inline $\\sigma_i:=\\sigma^{2}*s_i$}. If {@latex.inline $s_i=0$} then {@latex.inline $\\sigma_i:=\\sigma^{2}$}.
 * 
 * After a new set of parameters {@latex.inline $\\vec{\\lambda}^{t}$} has been drawn, the sampling process decides if this new set of parameters is accepted
 * according to the distribution {@latex.inline $P(\\vec{\\lambda}^{t})$} that we want to sample from.
 * Specifically, the parameters are accepted, iff.
 * {@latex.ilb \\[ \\alpha < \\frac{ P(\\vec{\\lambda}^{t})Q(\\vec{\\lambda}^{t-1} | \\vec{\\lambda}^{t}) }{P(\\vec{\\lambda}^{t-1}) Q(\\vec{\\lambda}^{t} | \\vec{\\lambda}^{t-1})},\\]}
 * where {@latex.inline $\\alpha$} is drawn from a uniform distribution in {@latex.inline $[0,1]$}, i.e. {@link Random#nextDouble()}.
 * Otherwise, the parameters are rejected and {@latex.inline $\\vec{\\lambda}^{t}:=\\vec{\\lambda}^{t-1}$}.
 * 
 * Since the Gaussian distribution is symmetric around its mean, {@latex.inline $Q(\\vec{\\lambda}^{t-1} | \\vec{\\lambda}^{t})=Q(\\vec{\\lambda}^{t} | \\vec{\\lambda}^{t-1})$}, both terms
 * cancel, and the acceptance probability depends only on the current and previous values of {@latex.inline $P$}.
 * 
 * All sampled parameters are stored to separate temporary files for each concurrent sampling run by an internal {@link SamplingComponent}. The contents of these files
 * are stored together with the remaining representation of the {@link SamplingScoreBasedClassifier}, if {@link SamplingScoreBasedClassifier#toXML()} is called, and, hence,
 * can be stored to a monolithic file containing all information for, e.g., later classification procedures.
 * 
 * For determining the length of the burn-in phase and, as a consequence, the beginning of the stationary phase, a {@link BurnInTest} can be provided to the constructor of the classifier.
 * 
 * @author Jan Grau
 *
 */
public abstract class SamplingScoreBasedClassifier extends AbstractScoreBasedClassifier {

	private static Random r = new Random();
	
	/**
	 * Parameters
	 */
	protected SamplingScoreBasedClassifierParameterSet params;
	
	/**
	 * {@link SamplingDifferentiableStatisticalModel}s
	 */
	protected SamplingDifferentiableStatisticalModel[] scoringFunctions;
	
	/**
	 * the currently accepted parameters
	 */
	protected double[] currentParameters;
	
	/**
	 * The initial parameters if set by {@link SamplingScoreBasedClassifier#setInitParameters(double[])}, <code>null</code> otherwise
	 */
	protected double[] initParameters;
	
	/**
	 * The score achieved using {@link SamplingScoreBasedClassifier#currentParameters}
	 */
	protected double currentScore;
	
	/**
	 * The previously accepted parameters, backup for rollbacks
	 */
	protected double[] previousParameters;	
	
	/**
	 * The last accepted parameters for all samplings, backup for iterative
	 * sampling when checking for {@link BurnInTest}
	 */
	protected double[][] lastParameters;
	
	/**
	 * The scores yielded for the parameters in <code>lastParameters</code>
	 */
	protected double[] lastScore;
	
	/**
	 * The groups of parameters sampled using {@link SamplingScheme#GROUPED}
	 */
	private int[][] groupedParameters;
	
	/**
	 * The parameters offsets of all {@link SamplingScoreBasedClassifier#scoringFunctions}
	 */
	private int[] parameterOffsets;
	
	/**
	 * The sampling variances for sampling from Gaussian distribution
	 */
	private double[] classVariances;
	
	/**
	 * The standard deviations for sampling from the Gaussian distribution for each parameter
	 */
	private double[] samplingSds;
	
	/**
	 * true if sampling has been accomplished
	 */
	private boolean isTrained;
	
	/**
	 * The {@link BurnInTest}, may be null for no test
	 */
	protected BurnInTest burnInTest;
	
	/**
	 * The length of the burn-in phase as determined by {@link SamplingScoreBasedClassifier#burnInTest}
	 */
	protected Integer burnInLength;
	
	/**
	 * The sampling component that handles the (temporary) parameter files.
	 */
	private DiffSMSamplingComponent samplingComponent;
	
	/**
	 * The directory for temporary files used for storing sampled parameters.
	 */
	private File tempDir;
	
	/**
	 * If true, the temporary parameter files are deleted on exit of the program (default).
	 */
	private boolean deleteOnExit = true;
	
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, params, "parameters" );
		XMLParser.appendObjectWithTags( xml, scoringFunctions, "scoringFunctions" );
		XMLParser.appendObjectWithTags( xml, currentParameters, "currentParameters" );
		XMLParser.appendObjectWithTags( xml, initParameters, "initParameters" );
		XMLParser.appendObjectWithTags( xml, currentScore, "currentScore" );
		XMLParser.appendObjectWithTags( xml, previousParameters, "previousParameters" );
		XMLParser.appendObjectWithTags( xml, lastParameters, "lastParameters" );
		XMLParser.appendObjectWithTags( xml, lastScore, "lastScore" );
		XMLParser.appendObjectWithTags( xml, groupedParameters, "groupedParameters" );
		XMLParser.appendObjectWithTags( xml, parameterOffsets, "parametersOffsets" );
		XMLParser.appendObjectWithTags( xml, classVariances, "classVariances" );
		XMLParser.appendObjectWithTags( xml, samplingSds, "samplingSds" );
		XMLParser.appendObjectWithTags( xml, isTrained, "isTrained" );
		XMLParser.appendObjectWithTags( xml, burnInTest, "burnInTest" );
		XMLParser.appendObjectWithTags( xml, burnInLength, "burnInLength" );
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		try{
			StringBuffer sb = sfsc.saveParameters();
			XMLParser.addTags( sb, "sampledParameters" );
			xml.append( sb );
		}catch(Exception e){
			throw new RuntimeException( e );
		}
		return xml;
	}
	
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		params = XMLParser.extractObjectForTags( xml, "parameters", SamplingScoreBasedClassifierParameterSet.class );
		scoringFunctions = XMLParser.extractObjectForTags( xml, "scoringFunctions", SamplingDifferentiableStatisticalModel[].class );
		currentParameters = XMLParser.extractObjectForTags( xml, "currentParameters", double[].class );
		initParameters = XMLParser.extractObjectForTags( xml, "initParameters", double[].class );
		currentScore = XMLParser.extractObjectForTags( xml, "currentScore", double.class );
		previousParameters = XMLParser.extractObjectForTags( xml, "previousParameters", double[].class );
		lastParameters = XMLParser.extractObjectForTags( xml, "lastParameters", double[][].class );
		lastScore = XMLParser.extractObjectForTags( xml, "lastScore", double[].class );
		groupedParameters = XMLParser.extractObjectForTags( xml, "groupedParameters", int[][].class );
		parameterOffsets = XMLParser.extractObjectForTags( xml, "parameterOffsets", int[].class );
		classVariances = XMLParser.extractObjectForTags( xml, "classVariances", double[].class );
		samplingSds = XMLParser.extractObjectForTags( xml, "samplingSds", double[].class );
		isTrained = XMLParser.extractObjectForTags( xml, "isTrained", boolean.class );
		burnInTest = XMLParser.extractObjectForTags( xml, "burnInTest", BurnInTest.class );
		burnInLength = XMLParser.extractObjectForTags( xml, "burnInLength", Integer.class );
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		try{
			sfsc.initForSampling( params.getNumberOfStarts() );
			StringBuffer sb = XMLParser.extractForTag( xml, "sampledParameters" );
			sfsc.createFiles( sb );
		}catch(Exception e){
			NonParsableException ex = new NonParsableException( e.getMessage() );
			ex.setStackTrace( e.getStackTrace() );
		}
	}
	
	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public SamplingScoreBasedClassifier(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	/**
	 * Creates a new {@link SamplingScoreBasedClassifier} using the parameters in <code>params</code>,
	 * a specified {@link BurnInTest} (or <code>null</code> for no burn-in test), a set of sampling variances,
	 * which may be different for each of the classes (in analogy to equivalent sample size for the Dirichlet distribution),
	 * and set set of {@link SamplingDifferentiableStatisticalModel}s for each of the classes.
	 * @param params the external parameters of this classifier
	 * @param burnInTest the burn-in test (or <code>null</code> for no burn-in test)
	 * @param classVariances the variances used for sampling for the parameters of each class
	 * @param scoringFunctions the scoring functions for each of the classes
	 * @throws CloneNotSupportedException if the scoring functions or the burn-in test could not be cloned
	 * @see de.jstacs.sampling.VarianceRatioBurnInTest
	 */
	protected SamplingScoreBasedClassifier(SamplingScoreBasedClassifierParameterSet params, BurnInTest burnInTest, double[] classVariances, SamplingDifferentiableStatisticalModel... scoringFunctions) throws CloneNotSupportedException{
		super(params.getAlphabetContainer(), params.getLength(), scoringFunctions.length);
		this.params = params.clone();
		this.scoringFunctions = ArrayHandler.clone(scoringFunctions);
		this.burnInLength = null;
		this.classVariances = classVariances.clone();
		if(burnInTest != null){
			this.burnInTest = burnInTest.clone();
		}
	}
	
	
	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		CategoricalResult[] res = new CategoricalResult[scoringFunctions.length + 1];
		res[0] = new CategoricalResult( "classifier", "a <b>short</b> description of the classifier", getInstanceName() );
		int i = 0;
		while( i < scoringFunctions.length ) {
			res[i + 1] = new CategoricalResult( "class info " + i, "some information about the scoring function for class "+i, scoringFunctions[i++].getInstanceName() );
		}
		return res;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumericalCharacteristics()
	 */
	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {

		NumericalResult[] pars = new NumericalResult[scoringFunctions.length];
		
		for( int i = 0; i < scoringFunctions.length; i++ ) {
			pars[i] = new NumericalResult( "Number of parameters " + ( i + 1 ),
					"The number of parameters for scoring function " + ( i + 1 ) + ", -1 indicates unknown number of parameters.",
					scoringFunctions[i].getNumberOfParameters() );
		}
		return new NumericalResultSet( pars );
	}
	
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}
	
	/**
	 * Returns the function that should be sampled from.
	 * 
	 * @param data
	 *            the samples
	 * @param weights
	 *            the weights of the sequences of the samples
	 * 
	 * @return the function that should be sampled from
	 * 
	 * @throws Exception
	 *             if the function could not be created
	 */
	protected abstract Function getFunction( DataSet[] data, double[][] weights ) throws Exception;
	
	/**
	 * Allows for a modification of the value returned by the function
	 * obtained by {@link SamplingScoreBasedClassifier#getFunction(DataSet[], double[][])}.
	 * This is for instance necessary in case of {@link de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction} to
	 * obtain a proper posterior or supervised posterior.
	 * @param value the original value
	 * @return the modified value
	 */
	protected double modifyFunctionValue(double value){
		return value;
	}
	
	/**
	 * Returns a sampling component suited for this {@link SamplingScoreBasedClassifier}
	 * @return the sampling component
	 */
	protected DiffSMSamplingComponent getSamplingComponent() {
		if(samplingComponent == null){
			samplingComponent = new DiffSMSamplingComponent( params.getOutfilePrefix() );
		}
		return samplingComponent;
	}
	
	/**
	 * Returns the directory for parameter files set in this {@link SamplingScoreBasedClassifier}.
	 * If this value is <code>null</code>, the default directory of the executing OS is used for the parameter
	 * files.
	 * @return the temp directory
	 */
	public File getTempDir() {
		return tempDir;
	}

	/**
	 * Sets the directory for parameter files set in this {@link SamplingScoreBasedClassifier}.
	 * If <code>tempDir</code> is <code>null</code>, the default directory of the executing OS is used for the parameter
	 * files. If this value is reset after training, all sampled parameters will be lost.
	 * The value set by this method is not stored in the XML-representation.
	 * @param tempDir the temp directory
	 */
	public void setTempDir( File tempDir ) {
		this.tempDir = tempDir;
		samplingComponent = null;
		isTrained = false;
	}

	
	/**
	 * Returns <code>true</code> if the temporary parameter files shall
	 * be deleted on exit of the program.
	 * @return if temp files are deleted
	 */
	public boolean getDeleteOnExit() {
		return deleteOnExit;
	}

	/**
	 * If set to <code>true</code> (which is the default), the temporary files for storing sampled parameter
	 * values are deleted on exit of the program. If this value is set to <code>true</code> it cannot be
	 * reset to <code>false</code>, again, after sampling started due to the restrictions of {@link File#deleteOnExit()}.
	 * If you want to retain those
	 * parameters, nonetheless, you can call {@link SamplingScoreBasedClassifier#toXML()}
	 * and save this {@link StringBuffer}, which also contains the sampled
	 * parameter values, somewhere.
	 * The value set by this method is not stored in the XML-representation.
	 * @param deleteOnExit if temp files shall be deleted on exit
	 * @throws Exception if set to <code>false</code> after sampling started
	 */
	public void setDeleteOnExit( boolean deleteOnExit ) throws Exception {
		if(samplingComponent != null && deleteOnExit == false){
			throw new Exception("Cannot revoke delete on exit after creating files.");
		}
		this.deleteOnExit = deleteOnExit;
	}

	/**
	 * Initializes all internal fields and initializes the {@link SamplingScoreBasedClassifier#scoringFunctions}s randomly
	 * @param starts number of starts
	 * @param adaptVariance if true, variance is adapted to size of event space
	 * @param outfilePrefix the prefix of the outfiles
	 * @throws Exception if the scoring functions could not be initialized
	 */
	protected void init(int starts, boolean adaptVariance, String outfilePrefix) throws Exception {
		boolean freeParams = params.getFreeParameters();
		int numParams = 0;
		numParams += scoringFunctions.length - (freeParams ? 1 : 0);
		LinkedList<int[]> list = new LinkedList<int[]>();
		
		int[] temp2 = new int[scoringFunctions.length-(freeParams ? 1 : 0)+1];
		for(int i=0;i<temp2.length;i++){
			temp2[i] = i-1;
		}
		list.add( temp2 );
		parameterOffsets = new int[scoringFunctions.length+1];
		for(int i=0;i<scoringFunctions.length;i++){
			parameterOffsets[i] = numParams;
			scoringFunctions[i].initializeFunctionRandomly( freeParams );
			int[][] temp = scoringFunctions[i].getSamplingGroups( numParams );
			for(int j=0;j<temp.length;j++){
				temp2 = new int[temp[j].length+1];
				temp2[0] = i;
				for(int k=0;k<temp[j].length;k++){
					temp2[k+1] = temp[j][k];
				}
				list.add( temp2 );
			}
			numParams += scoringFunctions[i].getNumberOfParameters();
		}
		parameterOffsets[parameterOffsets.length-1] = numParams;
		groupedParameters = list.toArray( new int[0][0] );

		currentParameters = new double[numParams];
		samplingSds = new double[numParams];
		double[] temp = getClassWeights();//TODO init?
		int j=0;
		for(;j<temp.length-(freeParams ? 1 : 0);j++){
			if(freeParams){
				currentParameters[j] = temp[j] - temp[temp.length-1];
			}else{
				currentParameters[j] = temp[j];
			}
			samplingSds[j] = Math.sqrt( classVariances[j] );
		}
		for(int i=0;i<scoringFunctions.length;i++){
			temp = scoringFunctions[i].getCurrentParameterValues();
			System.arraycopy( temp, 0, currentParameters, j, temp.length );
			for(int k=0;k<temp.length;k++){
				if(adaptVariance){
					try{
						samplingSds[j+k] = Math.sqrt( classVariances[i]*scoringFunctions[i].getSizeOfEventSpaceForRandomVariablesOfParameter( k ) );
						if(Double.isNaN( samplingSds[j+k] )){
							samplingSds[j+k] = Math.sqrt( classVariances[i] );
						}
					}catch(Exception e){
						samplingSds[j+k] = Math.sqrt( classVariances[i] );
					}
				}else{
					samplingSds[j+k] = Math.sqrt( classVariances[i] );
				}
			}
			j += temp.length;
		}
		if(initParameters != null && initParameters.length == currentParameters.length){
			currentParameters = initParameters.clone();
		}
		previousParameters = currentParameters.clone();
		
		lastScore = new double[starts];
		Arrays.fill( lastScore, Double.NEGATIVE_INFINITY );
		lastParameters = new double[starts][];
		lastParameters[0] = currentParameters.clone();
		for(int i=1;i<starts;i++){
			lastParameters[i] = currentParameters.clone();
			for(int k=0;k<scoringFunctions.length;k++){
				scoringFunctions[k].initializeFunctionRandomly( freeParams );
				System.arraycopy( scoringFunctions[k].getCurrentParameterValues(), 0, lastParameters[i], parameterOffsets[k], parameterOffsets[k+1] - parameterOffsets[k] );
			}
		}
		
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		sfsc.initForSampling( starts );
	}
	
	
	/**
	 * Samples a predefined number of steps appended to the current sampling
	 * @param function the objective function
	 * @param component the sampling component with selected sampling
	 * @param test the burn-in test
	 * @param numSteps the number of steps
	 * @param scheme the {@link SamplingScheme}
	 * @return the value of the last accepted parameters
	 * @throws Exception if either the function could not be evaluated on the current parameters or the 
	 * 					sampled parameters could not be stored
	 */
	protected double sampleNSteps( Function function, DiffSMSamplingComponent component, BurnInTest test, int numSteps, SamplingScheme scheme ) throws Exception{
		double previousValue, newValue;
		if(currentScore == Double.NEGATIVE_INFINITY){
			newValue = modifyFunctionValue( function.evaluateFunction( currentParameters ) );
		}else{
			newValue = currentScore;
		}
		for(int i=0;i<numSteps;i++){
			previousValue = newValue;
			newValue = doOneSamplingStep( function, scheme, previousValue );
			if(Double.isNaN( newValue )){
				newValue = previousValue;
			}
			currentScore = newValue;
			component.acceptParameters();
			if(test != null){
				test.setValue( newValue );
			}
		}
		return newValue;
	}
	
	/**
	 * Samples as many steps as needed to get into the stationary phase according to
	 * {@link SamplingScoreBasedClassifier#burnInTest} and then samples the number of
	 * stationary steps as set in {@link SamplingScoreBasedClassifier#params}.
	 * @param sfsc the current sampling component
	 * @param function the objective function
	 * @throws Exception if the sampling could not be extended, e.g. due to evaluation errors
	 */
	protected void sample(DiffSMSamplingComponent sfsc, Function function) throws Exception{
		boolean afterBurnIn = false;
		int numIterations = 0;
		int starts = params.getNumberOfStarts(), numberOfTestIterations = params.getNumberOfTestSamplings(), numberOfStationaryIterations = params.getNumberOfStationarySamplings();
		SamplingScheme scheme = params.getSamplingScheme();
		boolean first = true;
		while(!afterBurnIn){
			for(int i=0;i<starts;i++){
				sfsc.extendSampling( i, !first );
				burnInTest.setCurrentSamplingIndex( i );
				sampleNSteps( function, sfsc, burnInTest, numberOfTestIterations, scheme );
			}
			numIterations += numberOfTestIterations;
			burnInLength = burnInTest.getLengthOfBurnIn();
			if(numIterations > burnInLength){
				afterBurnIn = true;
			}
			first=false;
		}
		for(int i=0;i<starts;i++){
			sfsc.extendSampling( i, true );
			sampleNSteps( function, sfsc, burnInTest, numberOfStationaryIterations-(numIterations-burnInLength), scheme );
		}
		
	}
	

	/**
	 * Performs one sampling step, i.e., one sampling of all parameter values.
	 * @param function the objective function
	 * @param scheme the {@link SamplingScheme}
	 * @param previousValue the value of the last sampling or minus infinity
	 * 			for the first sampling run
	 * @return the value of the last accepted parameter values or {@link Double#NaN} if none
	 * 			of the sampled parameters where accepted
	 * @throws Exception if the function could not be evaluated or an unknown {@link SamplingScheme} was provided
	 */
	protected double doOneSamplingStep( Function function, SamplingScheme scheme, double previousValue ) throws Exception{
		double returnValue = Double.NaN, temp;
		switchPars( 0, currentParameters.length, false );
		for(int i=0;i<scoringFunctions.length;i++){
			currentParameters[i] = r.nextGaussian()*samplingSds[i] + previousParameters[i];
		}
		switch( scheme ) {
			case ALL_PARAMETERS:
				for(int i=scoringFunctions.length;i<currentParameters.length;i++){
					currentParameters[i] = r.nextGaussian()*samplingSds[i] + previousParameters[i];
				}
				returnValue = testParameters( function, previousValue );
				if(Double.isNaN( returnValue )){
					switchPars( 0, currentParameters.length, true );
				}
				break;
			case FUNCTION_WISE:
				temp = testParameters( function, previousValue );
				if(!Double.isNaN( temp )){
					previousValue = temp;
					returnValue = temp;
				}
				switchPars( 0, scoringFunctions.length, Double.isNaN( temp ) );
				for(int i=0;i<scoringFunctions.length;i++){
					for(int j=parameterOffsets[i];j<parameterOffsets[i+1];j++){
						currentParameters[i] = r.nextGaussian()*samplingSds[j] + previousParameters[j];
					}
					temp = testParameters( function, previousValue );
					if(!Double.isNaN( temp )){
						previousValue = temp;
						returnValue = temp;
					}
					switchPars( parameterOffsets[i], parameterOffsets[i+1], Double.isNaN( temp ) );
				}
				break;
			case GROUPED:
				temp = testParameters( function, previousValue );
				if(!Double.isNaN( temp )){
					previousValue = temp;
					returnValue = temp;
				}
				switchPars( 0, scoringFunctions.length, Double.isNaN( temp ) );
				int idx;
				for(int i=0;i<groupedParameters.length;i++){
					for(int j=1;j<groupedParameters[i].length;j++){
						idx = groupedParameters[i][j];
						currentParameters[idx] = r.nextGaussian()*samplingSds[idx] + previousParameters[idx];
					}
					temp = testParameters( function, previousValue );
					if(!Double.isNaN( temp )){
						previousValue = temp;
						returnValue = temp;
					}
					switchPars( groupedParameters[i], Double.isNaN( temp ) );
				}
				break;
			case INDIVIDUAL:
				temp = testParameters( function, previousValue );
				if(!Double.isNaN( temp )){
					previousValue = temp;
					returnValue = temp;
				}
				for(int j=0;j<currentParameters.length;j++){
					currentParameters[j] = r.nextGaussian()*samplingSds[j] + previousParameters[j];
					temp = testParameters( function, previousValue );
					if(!Double.isNaN( temp )){
						previousValue = temp;
						returnValue = temp;
					} else {
						currentParameters[j] = previousParameters[j];
					}
				}
				break;
			default:
				throw new Exception( "Sampling scheme not implemented." );
		}
		return returnValue;
	}
	
	/**
	 * Test the current parameter values for acceptance.
	 * @param function the objective function
	 * @param previousValue the value of the previously accepted parameters
	 * @return the value of the newly accepted parameters or {@link Double#NaN} if the parameters where
	 * 			not accepted
	 * @throws DimensionException if the current parameters do not fit the function
	 * @throws EvaluationException if the function could not be evaluated due to other problems, e.g. infinite values
	 */
	private double testParameters( Function function, double previousValue) throws DimensionException, EvaluationException{
		double newValue = modifyFunctionValue( function.evaluateFunction( currentParameters ) );
		//System.out.println(newValue+" "+previousValue+" "+Math.exp(newValue - previousValue));
		if(Math.log( r.nextDouble() ) < newValue - previousValue){
			//System.out.println("a");
			return newValue;
		}else{
			return Double.NaN;
		}
	}
	
	/**
	 * Exchanges the parameter values from <code>start</code> to <code> end
	 * between {@link SamplingScoreBasedClassifier#currentParameters} and {@link SamplingScoreBasedClassifier#previousParameters}. If
	 * <code>rollback</code> is <code>true</code>, the values are copied form the previous parameter to the current parameters, and vice versa
	 * otherwise.
	 * @param start the first position (inclusive)
	 * @param end the last position (exclusive)
	 * @param rollback do a rollback or commit the currently sampled parameters
	 */
	private void switchPars(int start, int end, boolean rollback){
		if(rollback){
			System.arraycopy( previousParameters, start, currentParameters, start, end-start );
		}else{
			System.arraycopy( currentParameters, start, previousParameters, start, end-start );
		}
	}
	
	/**
	 * Exchanges the parameter values from <code>start</code> to <code> end
	 * between {@link SamplingScoreBasedClassifier#currentParameters} and {@link SamplingScoreBasedClassifier#previousParameters}. If
	 * <code>rollback</code> is <code>true</code>, the values are copied form the previous parameter to the current parameters, and vice versa
	 * otherwise.
	 * @param idxs the indexes to be copied
	 * @param rollback do a rollback or commit the currently sampled parameters
	 */
	private void switchPars(int[] idxs, boolean rollback){
		if(rollback){
			for(int i=1;i<idxs.length;i++){
				currentParameters[idxs[i]] = previousParameters[idxs[i]];
			}
		}else{
			for(int i=1;i<idxs.length;i++){
				previousParameters[idxs[i]] = currentParameters[idxs[i]];
			}
		}
	}
	
	@Override
	protected double getScore( Sequence seq, int cls, boolean check ) throws IllegalArgumentException, NotTrainedException, Exception {
		if(check){
			super.check( seq );
		}
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		int starts = params.getNumberOfStarts();
		sfsc.samplingStopped();
		if(burnInLength == null){
			precomputeBurnInLength( sfsc );
		}
		double score = 0, n = 0;
		for(int i=0;i<starts;i++){
			sfsc.parseParameterSet( i, burnInLength );
			while(sfsc.parseNextParameterSet()){
				setParameters(currentParameters);
				score += getClassWeight( cls ) + scoringFunctions[cls].getLogScoreFor( seq );
				n++;
			}
		}
		return score/n;
	}
	
	@Override
	public double[] getScores( DataSet s ) throws Exception {
		if( scoringFunctions.length != 2 ) {
			throw new OperationNotSupportedException( "This method is only for 2-class-classifiers." );
		}
		if( s == null ) {
			return new double[0];
		}
		check( s );
		int starts = params.getNumberOfStarts();
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		sfsc.samplingStopped();
		if(burnInLength == null){
			precomputeBurnInLength( sfsc );
		}
		double[] scores = new double[s.getNumberOfElements()];
		double n = 0;
		for(int i=0;i<starts;i++){
			sfsc.parseParameterSet( i, burnInLength );
			while(sfsc.parseNextParameterSet()){
				setParameters(currentParameters);
				for(int j=0;j<scores.length;j++){
					Sequence seq = s.getElementAt( j );
					scores[j] += getClassWeight( 0 ) - getClassWeight( 1 ) + scoringFunctions[0].getLogScoreFor( seq ) - scoringFunctions[1].getLogScoreFor( seq );
				}
				n++;
			}
		}
		for(int i=0;i<scores.length;i++){
			scores[i] /= n;
		}
		return scores;
	}


	/**
	 * Sets the initial parameters of the sampling to <code>parameters</code>.
	 * @param parameters the initial parameters
	 */
	public void setInitParameters(double[] parameters){
		this.initParameters = parameters.clone();
	}
	
	/**
	 * Sets the current parameters for the class weights and in all scoring functions 
	 * @param currentParameters
	 */
	protected void setParameters( double[] currentParameters ) {
		this.setClassWeights( false, currentParameters, 0 );
		for(int i=0;i<scoringFunctions.length;i++){
			scoringFunctions[i].setParameters( currentParameters, parameterOffsets[i] );
		}
	}

	@Override
	public boolean isInitialized() {
		return isTrained;
	}
	
	/**
	 * Checks the data and the weight and returns a reduced variant of the data, where
	 * multiple occuring sequences are represented by their weight
	 * @param data the data
	 * @param weights the external weights
	 * @return the weighted daa
	 * @throws ClassDimensionException the data do not fit the number of scoring functions
	 * @throws WrongAlphabetException the alphabet of the data does not fit that of the scoring functions
	 * @throws WrongLengthException the lenght of the sequences does not match that of the scoring functions
	 */
	private Pair<DataSet[], double[][]> check(DataSet[] data, double[][] weights) throws ClassDimensionException, WrongAlphabetException, WrongLengthException{
		if( weights != null && data.length != weights.length ) {
			throw new IllegalArgumentException( "data and weights do not match" );
		}
		if( scoringFunctions.length != data.length ) {
			throw new ClassDimensionException();
		}
		if( weights == null ) {
			weights = new double[data.length][];
		}
		WeightedDataSetFactory wsf;
		DataSet[] reduced = new DataSet[data.length];
		double[][] newWeights = new double[data.length][];
		AlphabetContainer abc = getAlphabetContainer();
		for( int l = getLength(), i = 0; i < scoringFunctions.length; i++ ) {
			if( weights[i] != null && data[i].getNumberOfElements() != weights[i].length ) {
				throw new IllegalArgumentException( "At least for one sample: The dimension of the sample and the weight do not match." );
			}
			if( !abc.checkConsistency( data[i].getAlphabetContainer() ) ) {
				throw new IllegalArgumentException( "At least one sample is not defined over the correct alphabets." );
			}
			if( data[i].getElementLength() != l ) {
				// throw new IllegalArgumentException( "At least one sample has not the correct length." );
				wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[i], weights[i], l );
			} else {
				wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data[i], weights[i] );
			}
			reduced[i] = wsf.getDataSet();
			newWeights[i] = wsf.getWeights();
		}
		return new Pair<DataSet[], double[][]>( reduced, newWeights );
	}
	
	/**
	 * Does a single sampling run for a predefined number of steps.
	 * @param s the data
	 * @param weights the weights for the data
	 * @param numSteps the number of sampling steps
	 * @param outfilePrefix the prefix of the outfile where the parameter values
	 * 			are stored
	 * @throws Exception if the scoring functions could not be initialized or the sampling could not be extended, e.g. due to evaluation errors
	 */
	public void doSingleSampling( DataSet[] s, double[][] weights, int numSteps, String outfilePrefix ) throws Exception{
		Pair<DataSet[], double[][]> pair= check( s, weights );
		s = pair.getFirstElement();
		weights = pair.getSecondElement();
		init( 1, params.getAdaptVariance(), outfilePrefix );
		Function function = getFunction( s, weights );
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		sfsc.extendSampling( 0, false );
		sampleNSteps( function, sfsc, null, numSteps, params.getSamplingScheme() );
		sfsc.samplingStopped();
	}

	protected int getNumberOfParameters() {
		int res = scoringFunctions.length;
		for( int i = 0; i < scoringFunctions.length; i++ ) {
			res += scoringFunctions[i].getNumberOfParameters();
		}
		return res;
	}
	
	@Override
	public void train( DataSet[] s, double[][] weights ) throws Exception {
		Pair<DataSet[], double[][]> pair= check( s, weights );
		s = pair.getFirstElement();
		weights = pair.getSecondElement();
		init( params.getNumberOfStarts(), params.getAdaptVariance(), params.getOutfilePrefix() );
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		if(burnInTest != null){
			burnInTest.resetAllValues();
		}
		Function function = getFunction( s, weights );
		sample( sfsc, function );
		sfsc.samplingStopped();
		isTrained = true;
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		try {
			int starts = params.getNumberOfStarts();
			DiffSMSamplingComponent sfsc = getSamplingComponent();
			sfsc.samplingStopped();
			if(burnInLength == null){
				precomputeBurnInLength( sfsc );
			}
			System.out.println("burnInLength= " + burnInLength);
			out=sb;
			for(int i=0;i<starts;i++){
				sb.append( "start = " + i +"\n"); 
				sfsc.parseParameterSet( i, burnInLength );
				while(sfsc.parseNextParameterSet()){
					sb.append( currentScore + "\t" + Arrays.toString(currentParameters) + "\n" );
				}
				sb.append("\n");
			}
			out=null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return sb.toString();
	}
	
	static StringBuffer out = null;
		
	/**
	 * Precomputes the length of the burn-in phase, e.g. useful for computing scores of
	 * multiple sequences
	 * @param sfsc the current sampling component
	 * @throws Exception if the parameters values could not be parsed
	 */
	protected void precomputeBurnInLength(DiffSMSamplingComponent sfsc) throws Exception{
		if(burnInTest == null){
			burnInLength = 0;
			return;
		}
		int starts = params.getNumberOfStarts();
		for(int i=0;i<starts;i++){
			sfsc.parseParameterSet( i, 0 );
			burnInTest.setCurrentSamplingIndex( i );
			while(sfsc.parseNextParameterSet()){
				burnInTest.setValue( currentScore );
			}
		}
		burnInLength = burnInTest.getLengthOfBurnIn();
	}

	
	/**
	 * Returns the sampled parameter values with the maximum value of the objective function
	 * @return the best parameters
	 * @throws Exception if the parameters values could not be parsed
	 */
	public double[] getBestParameters() throws Exception {
		int starts = params.getNumberOfStarts();
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		
		double best = Double.NEGATIVE_INFINITY;
		double[] bestParameters = null;
		for(int i=0;i<starts;i++){
			sfsc.parseParameterSet( i, 0 );
			while(sfsc.parseNextParameterSet()){
				if(currentScore > best){
					if(bestParameters == null){
						bestParameters = currentParameters.clone();
					}else{
						System.arraycopy( currentParameters, 0, bestParameters, 0, currentParameters.length );
					}
					best = currentScore;
				}
			}
		}
		System.out.println(best);
		return bestParameters;
	}
	
	
	/**
	 * Returns the mean parameters over all samplings of all stationary phases.
	 * @param testBurnIn true if the length of the burn-in phase shall be computed
	 * @param minBurnInSteps minimum number of steps considered as burn-in
	 * @return the mean parameters
	 * @throws Exception if the parameters values could not be parsed
	 */
	protected double[] getMeanParameters(boolean testBurnIn, int minBurnInSteps) throws Exception {
		int starts = params.getNumberOfStarts();
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		
		if(testBurnIn && burnInLength == null){
			precomputeBurnInLength( sfsc );
		}
		double[] meanParameters = null;
		double n=0,k;
		for(int i=0;i<starts;i++){
			k=0;
			sfsc.parseParameterSet( i, 0 );
			while(sfsc.parseNextParameterSet()){
				if(k >= minBurnInSteps && (!testBurnIn || k > burnInLength) ){
					if(meanParameters == null){
						meanParameters = currentParameters.clone();
					}else{
						for(int j=0;j<currentParameters.length;j++){
							meanParameters[j] += currentParameters[j];
						}
					}
					n++;
				}
				k++;
			}
		}
		for(int i=0;i<meanParameters.length;i++){
			meanParameters[i] /= n;
		}
		return meanParameters;
	}
	
	/**
	 * Combines parameter files such that they are accepted as parameter files
	 * of this {@link SamplingScoreBasedClassifier}
	 * @param add if <code>true, parameter files are appended to the current ones, i.e., the number
	 * 			of samplings is augmented by these files
	 * @param files the parameter files
	 * @throws Exception if the parameter files could not be joined
	 */
	public void joinAndSetParameterFiles(boolean add, File... files ) throws Exception{
		
		if(lastParameters == null){
			init( params.getNumberOfStarts(), params.getAdaptVariance(), params.getOutfilePrefix() );
		}
		
		DiffSMSamplingComponent sfsc = getSamplingComponent();
		
		boolean b = sfsc.joinAndSetParameterFiles(add, files);
		if(b){
			if(add){
				params.setNumberOfStarts(params.getNumberOfStarts()+files.length);
			}else{
				params.setNumberOfStarts(files.length);
			}
		}
	}

	/**
	 * Sampling scheme for sampling the parameters of the scoring functions.
	 * While sampling all parameters in each step may lead to low acceptance rates,
	 * it may be suited for scoring functions with low number of parameters, or in cases
	 * where an intuitive grouping of parameters is not possible. On the other hand,
	 * grouped sampling may escape local maxima only after a large number of sampling steps.
	 * Sampling function-wise sampled parameters grouped by scoring functions.
	 * @author Jan Grau
	 * @see SamplingDifferentiableStatisticalModel#getSamplingGroups(int)
	 */
	public enum SamplingScheme{
		/**
		 * All parameters are drawn before decision of acceptance
		 */
		ALL_PARAMETERS,
		/**
		 * The parameters of each function are drawn before decision of acceptance
		 */
		FUNCTION_WISE,
		/**
		 * The parameters of each group are drawn before decision of acceptance
		 */
		GROUPED, 
		//for testing Jens
		INDIVIDUAL
	}
	
	/**
	 * The {@link SamplingComponent} that handles storing and loading sampled parameters values
	 * to and from files.
	 * 
	 * @author Jan Grau
	 *
	 */
	protected class DiffSMSamplingComponent implements SamplingComponent{

		/**
		 * The files for parameter output
		 */
		private File[] outfiles;
		/**
		 * The {@link PrintWriter} for the current sampling
		 */
		private PrintWriter curr;
		/**
		 * The {@link SparseStringExtractor} for parsing the parameters of
		 * the current sampling
		 */
		private SparseStringExtractor extract;
		/**
		 * The prefix of the outfiles
		 */
		private String outfilePrefix;
		
		/**
		 * If this {@link DiffSMSamplingComponent} is ready for storing sampled parameters
		 */
		private boolean inSamplingMode;
		/**
		 * the index of the current sampling
		 */
		private int currSampling;
		
		/**
		 * Creates a new {@link DiffSMSamplingComponent} that uses temporary files
		 * with name prefix <code>outfilePrefix</code> to store sampled parameters.
		 * These files are either stored in the default temp directory of the OS, or,
		 * if specified, in the directory set by {@link SamplingScoreBasedClassifier#setTempDir(File)} before
		 * starting the sampling process.
		 * @param outfilePrefix the prefix of the parameter files
		 */
		public DiffSMSamplingComponent(String outfilePrefix){
			this.inSamplingMode = false;
			this.outfilePrefix = outfilePrefix;
			this.currSampling = -1;
		}
		
		
		/**
		 * Combines parameter files such that they are accepted as parameter files
		 * of the calling {@link SamplingScoreBasedClassifier}
		 * @param add if true, parameter files are appended to the current ones, i.e., the number
		 * 			of samplings is augmented by these files
		 * @param files the parameter files
		 * @throws Exception if the parameter files could not be loaded
		 * @return <code>true</code> if the parameters could be joined
		 */
		public boolean joinAndSetParameterFiles( boolean add, File[] files ) throws Exception {
			int l=-1;
			double[][] newLastParameters = new double[files.length][];
			for(int i=0;i<files.length;i++){
				SparseStringExtractor ex = new SparseStringExtractor( files[i] );
				while(ex.hasMoreElements()){
					String str = ex.nextElement();
					String[] splt1 = str.split( "\t" );
					splt1[1] = splt1[1].substring( 1, splt1[1].length()-1 );
					String[] splt = splt1[1].split( ", " );
					if(l==-1){
						l = splt.length;
						if(currentParameters != null && currentParameters.length != l){
							throw new Exception( "Number of parameters does not match scoring functions: expected "+currentParameters.length+" but found "+splt.length );
						}
					}else if(splt.length != 0 && l != splt.length){
							throw new Exception( "Numbers of parameters between different files do not match: expected "+l+" but found "+splt.length );
					}
					if(!ex.hasMoreElements()){
						newLastParameters[i] = new double[splt.length];
						for(int j=0;j<newLastParameters[i].length;j++){
							newLastParameters[i][j] = Double.parseDouble( splt[j] );
						}
					}
				}
			}
			int off = 0;
			File[] temp;
			double[][] temp2;

			off = add ? outfiles.length : 0;
			temp = new File[off + files.length];
			temp2 = new double[temp.length][];
			
			if( add ) {
				System.arraycopy( outfiles, 0, temp, 0, off );
				System.arraycopy( lastParameters, 0, temp2, 0, off );
			}
			
			for(int i=0;i<files.length;i++){
				temp[i+off] = getOutfile( i+off );
				FileManager.copy( files[i].getAbsolutePath(), temp[i+off].getAbsolutePath() );
				temp2[i+off] = newLastParameters[i];
			}
			outfiles = temp;
			lastParameters = temp2;
			

			if(!add){
				System.arraycopy( lastParameters[0], 0, currentParameters, 0, currentParameters.length );
			}
			return true;
		}
		
		//TODO new
		public int getNumberOfParameterSets( int sampling ) throws Exception {
			currSampling = sampling;
			extract = new SparseStringExtractor( outfiles[sampling] );
			int i=0;
			while(extract.hasMoreElements() ){
				extract.nextElement();
				i++;
			}
			return i;
		}

		@Override
		public boolean parseParameterSet( int sampling, int n ) throws Exception {
			currSampling = sampling;
			extract = new SparseStringExtractor( outfiles[sampling] );
			int i=0;
			while(extract.hasMoreElements() && ++i < n){
				String s = extract.nextElement();
				/*if( out!=null ) {
					out.append(s+"\n");
				}/**/
			}
			/*if( out!=null ) {
				out.append("\n");
			}/**/
			return parseNextParameterSet();
		}

		@Override
		public boolean parseNextParameterSet() {
			if(!extract.hasMoreElements()){
				return false;
			}
			String str = extract.nextElement();
			String[] splt1 = str.split( "\t" );
			splt1[1] = splt1[1].substring( 1, splt1[1].length()-1 );
			String[] splt = splt1[1].split( ", " );
			if(currentParameters.length != splt.length){
				return false;
			}else{
				for(int j=0;j<splt.length;j++){
					currentParameters[j] = Double.parseDouble( splt[j] );
				}
				currentScore = Double.parseDouble( splt1[0] );
			}
			return true;
		}

		@Override
		public void initForSampling( int starts ) throws IOException {
			outfiles = new File[starts];
			for(int i=0;i<starts;i++){
				outfiles[i] = getOutfile( i );
			}
			inSamplingMode = true;
		}

		/**
		 * Returns the output file for sampling <code>idx</code>.
		 * @param idx the index of the sampling
		 * @return the output file
		 * @throws IOException if the file could not be created
		 */
		private File getOutfile(int idx) throws IOException {
			File f =  File.createTempFile( outfilePrefix, idx+".sam",tempDir );
			if(deleteOnExit){
				f.deleteOnExit();
			}
			return f;
			
		}
		
		@Override
		public void extendSampling( int sampling, boolean append ) throws IOException {
			if(currSampling >= 0){
				System.arraycopy( currentParameters, 0, lastParameters[currSampling], 0, currentParameters.length );
				lastScore[currSampling] = currentScore;
			}
			System.arraycopy( lastParameters[sampling], 0, currentParameters, 0, currentParameters.length );
			currentScore = lastScore[sampling];
			currSampling = sampling;
			curr = new PrintWriter( new FileOutputStream( outfiles[sampling], append ) );
		}

		@Override
		public void samplingStopped() throws IOException {
			if(curr != null){
				curr.close();
			}
			inSamplingMode = false;
		}

		@Override
		public boolean isInSamplingMode() {
			return inSamplingMode;
		}

		@Override
		public void acceptParameters() throws IOException {
			curr.println(currentScore+"\t"+Arrays.toString( currentParameters ));
			curr.flush();
		}

		
		/**
		 * Saves the parameter values of all parameter files to
		 * a {@link StringBuffer} representing these as XML.
		 * @return the {@link StringBuffer} containing the parameters
		 * @throws IOException if the files could not be read
		 */
		protected StringBuffer saveParameters() throws IOException{
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<outfiles.length;i++){
				StringBuffer temp = FileManager.readFile( outfiles[i] );
				XMLParser.addTagsAndAttributes( temp,"outfile","pos=\"i\"");
				sb.append( temp );
			}
			return sb;
		}
		
		/**
		 * Creates files out of file contents saved as XML.
		 * @param contents the parameter values to be stored in the files
		 * @throws NonParsableException if the parameters could not be parsed
		 * @throws IOException if the files could not be written
		 * @see SamplingScoreBasedClassifier.DiffSMSamplingComponent#saveParameters()
		 */
		protected void createFiles(StringBuffer contents) throws NonParsableException, IOException{
			HashMap<String, String> posFilter = new HashMap<String, String>();
			for(int i=0;i<outfiles.length;i++){
				posFilter.put( "pos", i+"" );
				StringBuffer temp = XMLParser.extractForTag( contents, "outfile", null, posFilter );
				FileManager.writeFile( outfiles[i], temp );
			}
		}
		
	}
	
}
