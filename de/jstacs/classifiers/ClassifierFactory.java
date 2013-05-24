package de.jstacs.classifiers;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModelFactory;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;

/**
 * This class allows to easily create classifiers from {@link TrainableStatisticalModel}s, {@link DifferentiableStatisticalModel}s, and {@link DifferentiableSequenceScore}s.
 * Most parameters of the classifiers are set to default values. If you like to set some of these additional parameters to non-standard values, directly use the constructors of 
 * {@link TrainSMBasedClassifier}, {@link GenDisMixClassifier}, or {@link MSPClassifier}.
 * 
 * @author Jan Grau
 *
 */
public class ClassifierFactory {

	/**
	 * Creates a classifier that is based on at least two {@link TrainableStatisticalModel}s. Such models can be created using the 
	 * {@link TrainableStatisticalModelFactory} or by directly using constructors of sub-classes of {@link TrainableStatisticalModel}.
	 * After the classifier has been created, it can be trained using the {@link AbstractClassifier#train(de.jstacs.data.DataSet...)} method using
	 * one data set for each {@link TrainableStatisticalModel} provided.
	 * 
	 * @param models the models for the individual classes
	 * @return the classifier
	 * @throws IllegalArgumentException
	 *             if the {@link TrainableStatisticalModel}s do not describe a common domain of
	 *             sequences
	 * @throws CloneNotSupportedException
	 *             if at least one {@link TrainableStatisticalModel} could not be cloned
	 * @throws ClassDimensionException
	 *             if the number of classes is below 2
	 */
	public static AbstractClassifier createGenerativeClassifier(TrainableStatisticalModel... models) throws IllegalArgumentException, CloneNotSupportedException, ClassDimensionException{
		return new TrainSMBasedClassifier( models );
	}
	
	
	/**
	 * Creates a classifier that is based on at least two {@link DifferentiableStatisticalModel}s. Such models can be created using the 
	 * {@link DifferentiableStatisticalModelFactory} or by directly using constructors of sub-classes of {@link DifferentiableStatisticalModel}.
	 * After the classifier has been created, it can be trained using the {@link AbstractClassifier#train(de.jstacs.data.DataSet...)} method using
	 * one data set for each {@link DifferentiableStatisticalModel} provided. The {@link LearningPrinciple} (generative ML or MAP, 
	 * or discriminative MCL or MSP) must be provided as first parameter.
	 *  
	 * @param principle the learning principle
	 * @param models the models for the individual classes
	 * @return the classifier
	 * @throws IllegalArgumentException if the lengths of the {@link DifferentiableStatisticalModel}s (see {@link DifferentiableStatisticalModel#getLength()}) are incompatible
	 * @throws Exception if something else went wrong
	 */
	public static AbstractClassifier createClassifier(LearningPrinciple principle, DifferentiableStatisticalModel... models) throws IllegalArgumentException, Exception{
		int length = models[0].getLength();
		for(int i=1;i<models.length;i++){
			int l = models[i].getLength();
			if(l != length){
				if(length == 0){
					length = l;
				}else if(l != 0){
					throw new IllegalArgumentException( "Model "+i+" has different but fixed length than a previous model." );
				}
			}
		}
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( models[0].getAlphabetContainer(), length, Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, 1 );
		return new GenDisMixClassifier( params, new CompositeLogPrior(), principle, models );
	}
	
	/**
	 * Creates a classifier that is based on at least two {@link DifferentiableSequenceScore}s. Such scores can be created using the 
	 * {@link DifferentiableStatisticalModelFactory} (since {@link DifferentiableStatisticalModel}s are also {@link DifferentiableSequenceScore}s) 
	 * or by directly using constructors of sub-classes of {@link DifferentiableSequenceScore}.
	 * After the classifier has been created, it can be trained using the discriminative MCL principle using the {@link AbstractClassifier#train(de.jstacs.data.DataSet...)} method using
	 * one data set for each {@link DifferentiableSequenceScore} provided.
	 *  
	 * @param scores the scores for the individual classes
	 * @return the classifier
	 * @throws IllegalArgumentException if the lengths of the {@link DifferentiableSequenceScore}s (see {@link DifferentiableSequenceScore#getLength()}) are incompatible
	 * @throws Exception if something else went wrong
	 */
	public static AbstractClassifier createClassifier(DifferentiableSequenceScore... scores) throws IllegalArgumentException, Exception{
		int length = scores[0].getLength();
		for(int i=1;i<scores.length;i++){
			int l = scores[i].getLength();
			if(l != length){
				if(length == 0){
					length = l;
				}else if(l != 0){
					throw new IllegalArgumentException( "Model "+i+" has different but fixed length than a previous model." );
				}
			}
		}
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( scores[0].getAlphabetContainer(), length, Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, 1 );
		return new MSPClassifier( params, scores );
	}
	
	
}
