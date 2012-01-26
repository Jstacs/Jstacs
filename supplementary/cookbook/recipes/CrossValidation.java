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
package supplementary.cookbook.recipes;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.assessment.KFoldCrossValidation;
import de.jstacs.classifiers.assessment.KFoldCrossValidationAssessParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.MixtureTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization;

/**
 * This class exemplarily shows how to perform a cross validation.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class CrossValidation {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the foreground data set</li>
	 * <li>args[1] contains the path to the background data set</li>
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		
		//the length of our input sequences
		int length = data[0].getElementLength();

		data[1] = new DataSet( new DNADataSet( args[1] ), length );
		 
		AlphabetContainer container = data[0].getAlphabetContainer();
		
		//create a new PWM
		BayesianNetworkTrainSM pwm = new BayesianNetworkTrainSM( new BayesianNetworkTrainSMParameterSet(
				//the alphabet and the length of the model:
				container, length, 
				//the equivalent sample size to compute hyper-parameters
				4, 
				//some identifier for the model
				"my PWM", 
				//we want a PWM, which is an inhomogeneous Markov model (IMM) of order 0
				ModelType.IMM, (byte) 0, 
				//we want to estimate the MAP-parameters
				LearningType.ML_OR_MAP ) );
		 
		//create a new mixture model using 2 PWMs
		MixtureTrainSM mixPwms = new MixtureTrainSM(
				//the length of the mixture model
				length, 
				//the two components, which are PWMs
				new TrainableStatisticalModel[]{pwm,pwm},
				//the number of starts of the EM
				10,
				//the equivalent sample sizes
				new double[]{pwm.getESS(),pwm.getESS()},
				//the hyper-parameters to draw the initial sequence-specific component weights (hidden variables)
				1,
				//stopping criterion
				new SmallDifferenceOfFunctionEvaluationsCondition(1E-6),
				//parameterization of the model, LAMBDA complies with the
				//parameterization by log-probabilities
				Parameterization.LAMBDA);
		 
		//create a new inhomogeneous Markov model of order 3
		BayesianNetworkTrainSM mm = new BayesianNetworkTrainSM( 
				new BayesianNetworkTrainSMParameterSet( container, length, 256, "my iMM(3)", ModelType.IMM, (byte) 3, LearningType.ML_OR_MAP ) );
		 
		//create a new PWM scoring function
		BayesianNetworkDiffSM dPwm = new BayesianNetworkDiffSM(
				//the alphabet and the length of the scoring function
				container, length, 
				//the equivalent sample size for the plug-in parameters
				4, 
				//we use plug-in parameters
				true, 
				//a PWM is an inhomogeneous Markov model of order 0
				new InhomogeneousMarkov(0));
		 
		//create a new mixture scoring function
		MixtureDiffSM dMixPwms = new MixtureDiffSM(
				//the number of starts
				2,
				//we use plug-in parameters
				true,
				//the two components, which are PWMs
				dPwm,dPwm);
		 
		//create a new scoring function that is an inhomogeneous Markov model of order 3
		BayesianNetworkDiffSM dMm = new BayesianNetworkDiffSM(container, length, 4, true, new InhomogeneousMarkov(3));
		 
		//create the classifiers
		int threads = AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors();
		AbstractScoreBasedClassifier[] classifiers = new AbstractScoreBasedClassifier[]{
									   //model based with mixture model and Markov model
									   new TrainSMBasedClassifier( mixPwms, mm ),
									   //conditional likelihood based classifier
									   new MSPClassifier( new GenDisMixClassifierParameterSet(container, length, 
											   //method for optimizing the conditional likelihood and 
											   //other parameters of the numerical optimization
											   Optimizer.QUASI_NEWTON_BFGS, 1E-2, 1E-2, 1, true, KindOfParameter.PLUGIN, false, threads),
											   //mixture scoring function and Markov model scoring function
											   dMixPwms,dMm )
		};
		 
		//create an new k-fold cross validation using above classifiers
		KFoldCrossValidation cv = new KFoldCrossValidation( classifiers );
		 
		//we use a specificity of 0.999 to compute the sensitivity and a sensitivity of 0.95 to compute FPR and PPV
		NumericalPerformanceMeasureParameterSet mp = (NumericalPerformanceMeasureParameterSet) PerformanceMeasureParameterSet.createFilledParameters();
		//we do a 10-fold cross validation and partition the data by means of the number of symbols
		KFoldCrossValidationAssessParameterSet cvpars = new KFoldCrossValidationAssessParameterSet(PartitionMethod.PARTITION_BY_NUMBER_OF_SYMBOLS, length, true, 2);
		 
		//compute the result of the cross validation and print them to System.out
		System.out.println( cv.assess( mp, cvpars, data ) );
	}

}
