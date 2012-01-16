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

package de.jstacs.classifier.differentiableSequenceScoreBased.sampling;

import de.jstacs.NonParsableException;
import de.jstacs.classifier.differentiableSequenceScoreBased.DiffSSBasedOptimizableFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.io.XMLParser;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel;

/**
 * A classifier that samples its parameters from a {@link LogGenDisMixFunction} using the
 * Metropolis-Hastings algorithm. For details on the algorithm see {@link SamplingScoreBasedClassifier}.
 * 
 * The {@link LogGenDisMixFunction} includes several known posterior distributions, including the posterior ({@link LearningPrinciple#MAP})
 * and the supervised posterior ({@link LearningPrinciple#MSP}). For non-uniform values of the mixture parameters <code>beta</code> the distribution
 * we sample from is less well defined, although sampling is possible in general.
 * 
 * @author Jan Grau
 *
 */
public class SamplingGenDisMixClassifier extends SamplingScoreBasedClassifier {

	/**
	 * The prior on the parameter values
	 */
	private LogPrior prior;
	
	/**
	 * The weights of the three components of the
	 * {@link LogGenDisMixFunction}, i.e., likelihood,
	 * conditional likelihood, and prior
	 */
	private double[] beta;
	
	/**
	 * The factor the function value is multiplied with.
	 * For posterior and supervised posterior, this is 2, while
	 * for PGDT it is 3. If <code>beta</code> contains non-uniform
	 * values for the three components, this defaults to 1.
	 */
	private double factor;
	
	/**
	 * Creates a new {@link SamplingGenDisMixClassifier} using the external parameters
	 * <code>params</code>, a burn-in test, a set of sampling variances for the different classes,
	 * a prior on the parameters, weights <code>beta</code> for the three components of the
	 * {@link LogGenDisMixFunction}, i.e., likelihood, conditional likelihood, and prior,
	 * and scoring functions that model the distribution for each of the classes.
	 * @param params the external parameters
	 * @param burnInTest the burn-in test, or <code>null</code> for no burn-in test
	 * @param classVariances the sampling variances for the parameters in the different classes
	 * @param prior the prior on the parameters
	 * @param beta The weights of the three components of the {@link LogGenDisMixFunction}
	 * @param scoringFunctions the scoring functions for the different classes
	 * @throws CloneNotSupportedException if the scoring functions could not be cloned
	 */
	public SamplingGenDisMixClassifier( SamplingGenDisMixClassifierParameterSet params, BurnInTest burnInTest, double[] classVariances, LogPrior prior, double[] beta, SamplingDifferentiableStatisticalModel... scoringFunctions ) throws CloneNotSupportedException {
		super( params, burnInTest, classVariances, scoringFunctions );
		this.prior = prior;
		this.beta = beta.clone();
		this.factor = 0;
		for(int i=0;i<beta.length;i++){
			if(beta[i] > 0 && factor == 0){
				factor = 1.0/beta[i];
			}else if(beta[i] > 0 && factor > 0 && 1.0/factor != beta[i]){
				factor = 1.0;
			}
		}
	}
	
	/**
	 * Creates a new {@link SamplingGenDisMixClassifier} using the external parameters
	 * <code>params</code>, a burn-in test, a set of sampling variances for the different classes,
	 * a prior on the parameters, a learning principle,
	 * and scoring functions that model the distribution for each of the classes.
	 * @param params the external parameters
	 * @param burnInTest the burn-in test, or <code>null</code> for no burn-in test
	 * @param classVariances the sampling variances for the parameters in the different classes
	 * @param prior the prior on the parameters
	 * @param principle the learning principle, i.e., the objective function we sample from
	 * @param scoringFunctions the scoring functions for the different classes
	 * @throws CloneNotSupportedException if the scoring functions could not be cloned
	 */
	public SamplingGenDisMixClassifier( SamplingGenDisMixClassifierParameterSet params, BurnInTest burnInTest, double[] classVariances, LogPrior prior, LearningPrinciple principle, SamplingDifferentiableStatisticalModel... scoringFunctions ) throws CloneNotSupportedException {
		this(params,burnInTest,classVariances,prior,LearningPrinciple.getBeta( principle ), scoringFunctions);
	}

	/**
	 * Creates a new {@link SamplingGenDisMixClassifier} from its XML-representation
	 * @param xml the XML-representation
	 * @throws NonParsableException if <code>xml</code> could not be parsed
	 */
	public SamplingGenDisMixClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected DiffSSBasedOptimizableFunction getFunction( DataSet[] data, double[][] weights ) throws Exception {
		
		DiffSSBasedOptimizableFunction fun = new LogGenDisMixFunction( ((SamplingGenDisMixClassifierParameterSet) params).getNumberOfThreads(), scoringFunctions,
				data, weights, prior, beta, false, params.getFreeParameters() );
		fun.reset();
		return fun;
	}
	
	@Override
	protected double modifyFunctionValue( double value ){
		return factor*value;
	}


	@Override
	protected String getXMLTag() {
		return "SamplingGenDisMixClassifier";
	}
	
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, prior, "prior" );
		XMLParser.appendObjectWithTags( xml, beta, "beta" );
		return xml;
	}
	
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		prior = XMLParser.extractObjectForTags( xml, "prior", LogPrior.class );
		beta = XMLParser.extractObjectForTags( xml, "beta", double[].class );
	}
	
	/**
	 * Returns a standard, i.e., non-sampling, {@link GenDisMixClassifier}, where the parameters
	 * are set to those that yielded the maximum value of the objective functions among all sampled
	 * parameter values.
	 * @param params the external parameters of the {@link GenDisMixClassifier}
	 * @return the {@link GenDisMixClassifier} with set parameter values
	 * @throws Exception if the {@link GenDisMixClassifier} could not be created
	 */
	public GenDisMixClassifier getClassifierForBestParameters(GenDisMixClassifierParameterSet params) throws Exception{
		GenDisMixClassifier cls = new GenDisMixClassifier( params, prior, beta, scoringFunctions );
		cls.initUsingParameters( getBestParameters() );
		return cls;
	}
	
	/**
	 * Returns a standard, i.e., non-sampling, {@link GenDisMixClassifier}, where the parameters
	 * are set to the mean values over all sampled
	 * parameter values in the stationary phase.
	 * @param params the external parameters of the {@link GenDisMixClassifier}
	 * @param testBurnIn if burn-in phase is tested, otherwise parameters starting from index <code>minBurnInSteps</code> are considered
	 * @param minBurnInSteps the minimum number of steps before the stationary phase
	 * @return the {@link GenDisMixClassifier} with set parameter values
	 * @throws Exception if the {@link GenDisMixClassifier} could not be created
	 */
	public GenDisMixClassifier getClassifierForMeanParameters(GenDisMixClassifierParameterSet params, boolean testBurnIn, int minBurnInSteps) throws Exception{
		GenDisMixClassifier cls = new GenDisMixClassifier( params, prior, beta, scoringFunctions );
		cls.initUsingParameters( getMeanParameters( testBurnIn, minBurnInSteps ) );
		return cls;
	}

}
