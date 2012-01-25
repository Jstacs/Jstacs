package supplementary.cookbook.recipes;

import java.util.Arrays;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;

/**
 * This class exemplarily shows how to train and evaluate a {@link GenDisMixClassifier}.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see TrainSMBasedClassifierTest
 */
public class GenDisMixClassifierTest {

	/**
	 * @param args
	 * <ul>
	 * <li>args[0] contains the path to the foreground data set</li>
	 * <li>args[1] contains the path to the background data set</li>
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		//read FastA-files
		DataSet[] data = {
		         new DNADataSet( args[0] ),
		         new DNADataSet( args[1] )
		};
		AlphabetContainer container = data[0].getAlphabetContainer();
		int length = data[0].getElementLength();
		
		//equivalent sample size =^= ESS
		double essFg = 4, essBg = 4;
		//create DifferentiableSequenceScore, here PWM
		DifferentiableStatisticalModel pwmFg = new BayesianNetworkDiffSM( container, length, essFg, true, new InhomogeneousMarkov(0) );
		DifferentiableStatisticalModel pwmBg = new BayesianNetworkDiffSM( container, length, essBg, true, new InhomogeneousMarkov(0) );
		
		//create parameters of the classifier
		GenDisMixClassifierParameterSet cps = new GenDisMixClassifierParameterSet(
				container,//the used alphabets
				length,//sequence length that can be modeled/classified
				Optimizer.QUASI_NEWTON_BFGS, 1E-9, 1E-11, 1,//optimization parameter
				false,//use free parameters or all
				KindOfParameter.PLUGIN,//how to start the numerical optimization
				true,//use a normalized objective function
				AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors()//number of compute threads		
		);
		
		//create classifiers
		LearningPrinciple[] lp = LearningPrinciple.values();
		GenDisMixClassifier[] cl = new GenDisMixClassifier[lp.length+1];
		//elementary learning principles
		int i = 0;
		for( ; i < cl.length-1; i++ ){
			System.out.println( "classifier " + i + " uses " + lp[i] );
			cl[i] = new GenDisMixClassifier( cps, new CompositeLogPrior(), lp[i], pwmFg, pwmBg );
		}
		
		//use some weighted version of log conditional likelihood, log likelihood, and log prior
		double[] beta = {0.3,0.3,0.4};
		System.out.println( "classifier " + i + " uses the weights " + Arrays.toString( beta ) );
		cl[i] = new GenDisMixClassifier( cps, new CompositeLogPrior(), beta, pwmFg, pwmBg );
		
		//do what ever you like
		
		//e.g., train
		for( i = 0; i < cl.length; i++ ){
			cl[i].train( data );
		}
		
		//e.g., evaluate (normally done on a test data set)
		PerformanceMeasureParameterSet mp = PerformanceMeasureParameterSet.createFilledParameters();
		for( i = 0; i < cl.length; i++ ){
			System.out.println( cl[i].evaluate( mp, true, data ) );
		}
	}

}
