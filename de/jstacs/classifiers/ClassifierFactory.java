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

package de.jstacs.classifiers;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.sequenceScores.SequenceScore;
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
	 * @return the sequence length
	 * @throws IllegalArgumentException if the lengths of the {@link DifferentiableSequenceScore}s (see {@link DifferentiableSequenceScore#getLength()}) are incompatible
	 */
	private static int getLength( SequenceScore...models ) {
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
		return length;
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
		return createClassifier(LearningPrinciple.getBeta(principle), models);
	}
	
	/**
	 * Creates a classifier that is based on at least two {@link DifferentiableStatisticalModel}s. Such models can be created using the 
	 * {@link DifferentiableStatisticalModelFactory} or by directly using constructors of sub-classes of {@link DifferentiableStatisticalModel}.
	 * After the classifier has been created, it can be trained using the {@link AbstractClassifier#train(de.jstacs.data.DataSet...)} method using
	 * one data set for each {@link DifferentiableStatisticalModel} provided. The <code>beta</code> array determines the employed learning principle
	 * and must be provided as first parameter.
	 *  
	 * @param beta an array specifying the learning principle by weighting factors for conditional likelihood, likelihood and prior
	 * @param models the models for the individual classes
	 * @return the classifier
	 * @throws IllegalArgumentException if the lengths of the {@link DifferentiableStatisticalModel}s (see {@link DifferentiableStatisticalModel#getLength()}) are incompatible
	 * @throws Exception if something else went wrong
	 * 
	 * @see LearningPrinciple
	 * @see LogGenDisMixFunction
	 */
	public static AbstractClassifier createClassifier(double[] beta, DifferentiableStatisticalModel... models) throws IllegalArgumentException, Exception{
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( models[0].getAlphabetContainer(), getLength(models), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, 1 );
		return new GenDisMixClassifier( params, new CompositeLogPrior(), beta, models );
	}
	
	/**
	 * Creates a classifier that is based on at least two {@link DifferentiableSequenceScore}s. Such scores can be created using the 
	 * {@link DifferentiableStatisticalModelFactory} (since {@link DifferentiableStatisticalModel}s are also {@link DifferentiableSequenceScore}s) 
	 * or by directly using constructors of sub-classes of {@link DifferentiableSequenceScore}.
	 * After the classifier has been created, it can be trained using the discriminative MCL principle using the {@link AbstractClassifier#train(de.jstacs.data.DataSet...)} method using
	 * one data set for each {@link DifferentiableSequenceScore} provided.
	 *  
	 * @param models the models for the individual classes
	 * @return the classifier
	 * @throws IllegalArgumentException if the lengths of the {@link DifferentiableSequenceScore}s (see {@link DifferentiableSequenceScore#getLength()}) are incompatible
	 * @throws Exception if something else went wrong
	 */
	public static AbstractClassifier createClassifier(DifferentiableSequenceScore... models) throws IllegalArgumentException, Exception{
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( models[0].getAlphabetContainer(), getLength(models), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, 1 );
		return new MSPClassifier( params, models );
	}
}
