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

package de.jstacs.classifier;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.classifier.MeasureParameters.Measure;
import de.jstacs.classifier.ScoreBasedPerformanceMeasureDefinitions.ThresholdMeasurePair;
import de.jstacs.classifier.utils.PValueComputation;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sample.ElementEnumerator;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.RangeParameter;
import de.jstacs.results.ImageResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.utils.REnvironment;

/**
 * This class is the main class for all score based classifiers. Score based
 * classifiers enable you to compute many different measures easily. For
 * instance one can use the package &quot;ROCR&quot; in R to compute or plot
 * many of them.
 * 
 * @author Jens Keilwagen, Jan Grau, Andre Gohr
 */
public abstract class AbstractScoreBasedClassifier extends AbstractClassifier {

	/**
	 * The weights of the classes, do not have to be probabilities.
	 */
	private double[] classWeights;

	/**
	 * The constructor for a homogeneous classifier. Such a classifier can
	 * handle sequences of arbitrary length. The class weights are set initially
	 * to <code>0</code>.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * @param classes
	 *            the number of different classes
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(AlphabetContainer,
	 *      int, int, double)
	 */
	public AbstractScoreBasedClassifier( AlphabetContainer abc, int classes ) {
		this( abc, 0, classes, 0d );
	}

	/**
	 * The constructor for a homogeneous classifier. Such a classifier can
	 * handle sequences of arbitrary length. The class weights are set initially
	 * to <code>classWeight</code>.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * @param classes
	 *            the number of different classes
	 * @param classWeight
	 *            the value of all class weights
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(AlphabetContainer,
	 *      int, int, double)
	 */
	public AbstractScoreBasedClassifier( AlphabetContainer abc, int classes, double classWeight ) {
		this( abc, 0, classes, classWeight );
	}

	/**
	 * The constructor for an inhomogeneous classifier. Such a classifier can
	 * handle sequences of fixed length. The class weights are set initially to
	 * <code>0</code>.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * @param length
	 *            the length of the sequences that can be classified
	 * @param classes
	 *            the number of different classes
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(AlphabetContainer,
	 *      int, int, double)
	 */
	public AbstractScoreBasedClassifier( AlphabetContainer abc, int length, int classes ) {
		this( abc, length, classes, 0d );
	}

	/**
	 * The constructor for an inhomogeneous classifier. Such a classifier can
	 * handle sequences of fixed length. The class weights are set initially to
	 * <code>classWeight</code>.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * @param length
	 *            the length of the sequences that can be classified
	 * @param classes
	 *            the number of different classes
	 * @param classWeight
	 *            the value of all class weights
	 * 
	 * @throws IllegalArgumentException
	 *             if the length and the possible length of the
	 *             {@link AlphabetContainer} does not match or the number of
	 *             classes is less than 2
	 * 
	 * @see AbstractClassifier#AbstractClassifier(AlphabetContainer, int)
	 */
	public AbstractScoreBasedClassifier( AlphabetContainer abc, int length, int classes, double classWeight )
																												throws IllegalArgumentException {
		super( abc, length );
		if( classes < 2 ) {
			throw new IllegalArgumentException( "You should have at least 2 classes." );
		}
		createDefaultClassWeights( classes, classWeight );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AbstractScoreBasedClassifier} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractScoreBasedClassifier} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see AbstractClassifier#AbstractClassifier(StringBuffer)
	 */
	public AbstractScoreBasedClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#clone()
	 */
	@Override
	public AbstractScoreBasedClassifier clone() throws CloneNotSupportedException {
		AbstractScoreBasedClassifier erg = (AbstractScoreBasedClassifier)super.clone();
		erg.classWeights = classWeights.clone();
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#classify(de.jstacs.data.Sequence)
	 */
	@Override
	public byte classify( Sequence seq ) throws Exception {
		return classify( seq, true );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#getResults(de.jstacs.data.Sample[], de.jstacs.classifier.MeasureParameters, boolean, boolean)
	 */
	@Override
	protected LinkedList<? extends Result> getResults( Sample[] s, MeasureParameters params, boolean exceptionIfNotComputeable ) throws Exception {
		if( s.length != getNumberOfClasses() ) {
			throw new ClassDimensionException();
		}
		if( s.length != 2 ) {
			return super.getResults( s, params, exceptionIfNotComputeable );
		} else {
			LinkedList<Result> list = new LinkedList<Result>();
			double[][] sortedScores = null;

			int numSelected = params.getNumberOfValues();
			if( numSelected > 0 ) {
				sortedScores = getSortedScore( s );
			}
			if( params.isSelected( Measure.ClassificationRate ) ) {
				list.add( new NumericalResult( Measure.ClassificationRate.getNameString(),
						Measure.ClassificationRate.getCommentString(),
						ScoreBasedPerformanceMeasureDefinitions.getClassificationRateFor2Classes( sortedScores[0], sortedScores[1] ) ) );
				//list.add( getClassificationRate( s ) );
				numSelected--;
			}
			if( params.isSelected( Measure.Sensitivity ) ) {
				double sp = (Double)( (ParameterSet)params.getValueFor( Measure.Sensitivity.getNameString() ) ).getParameterAt( 0 )
						.getValue();
				ThresholdMeasurePair sn = ScoreBasedPerformanceMeasureDefinitions.getSensitivityForSpecificity( sortedScores[0],
						sortedScores[1],
						sp );
				list.add( new NumericalResult( Measure.Sensitivity.getNameString(),
						Measure.Sensitivity.getCommentString() + " = " + sp,
						sn.getMeasure() ) );
				list.add( new NumericalResult( "Threshold for sensitivity",
						"The threshold for the sensitivity for a specificity of " + sp,
						sn.getThreshold() ) );
				numSelected--;
			}
			if( params.isSelected( Measure.FalsePositiveRate ) ) {
				double sn = (Double)( (ParameterSet)params.getValueFor( Measure.FalsePositiveRate.getNameString() ) ).getParameterAt( 0 )
						.getValue();
				ThresholdMeasurePair fpr = ScoreBasedPerformanceMeasureDefinitions.getFPRForSensitivity( sortedScores[0],
						sortedScores[1],
						sn );
				list.add( new NumericalResult( Measure.FalsePositiveRate.getNameString(),
						Measure.FalsePositiveRate.getCommentString() + " = " + sn,
						fpr.getMeasure() ) );
				list.add( new NumericalResult( "Threshold for false positive rate",
						"The threshold for the false positive rate for a fixed sensitivity of " + sn,
						fpr.getThreshold() ) );
				numSelected--;
			}
			if( params.isSelected( Measure.PositivePredictiveValue ) ) {
				double sn = (Double)( (ParameterSet)params.getValueFor( Measure.PositivePredictiveValue.getNameString() ) ).getParameterAt( 0 )
						.getValue();
				ThresholdMeasurePair ppv = ScoreBasedPerformanceMeasureDefinitions.getPPVForSensitivity( sortedScores[0],
						sortedScores[1],
						sn );
				list.add( new NumericalResult( Measure.PositivePredictiveValue.getNameString(),
						Measure.PositivePredictiveValue.getCommentString() + " = " + sn,
						ppv.getMeasure() ) );
				list.add( new NumericalResult( "Threshold for positive predictive value",
						"The threshold for the positive predictive value for a fixed sensitivity of " + sn,
						ppv.getThreshold() ) );
				numSelected--;
			}
			ArrayList<double[]> list2 = null;
			if( params.isSelected( Measure.AreaUnderROCCurve ) ) {
				if( params.isSelected( Measure.AreaUnderROCCurve ) ) {
					list2 = new ArrayList<double[]>( s[0].getNumberOfElements() / 2 );
				}
				list.add( new NumericalResult( Measure.AreaUnderROCCurve.getNameString(),
						Measure.AreaUnderROCCurve.getCommentString(),
						ScoreBasedPerformanceMeasureDefinitions.getAUC_ROC( sortedScores[0], sortedScores[1], list2 ) ) );
				numSelected--;
			}
			ArrayList<double[]> list3 = null;
			if( params.isSelected( Measure.AreaUnderPrecisionRecallCurve ) ) {
				if( params.isSelected( Measure.AreaUnderPrecisionRecallCurve ) ) {
					list3 = new ArrayList<double[]>( s[0].getNumberOfElements() / 2 );
				}
				list.add( new NumericalResult( Measure.AreaUnderPrecisionRecallCurve.getNameString(),
						Measure.AreaUnderPrecisionRecallCurve.getCommentString(),
						ScoreBasedPerformanceMeasureDefinitions.getAUC_PR( sortedScores[0], sortedScores[1], list3 ) ) );
				numSelected--;
			}
			if( params.isSelected( Measure.MaximumCorrelationCoefficient ) ) {
				double cc, t;
				try {
					ThresholdMeasurePair maxCC = ScoreBasedPerformanceMeasureDefinitions.getMaxOfCC( sortedScores[0], sortedScores[1] );
					cc = maxCC.getMeasure();
					t = maxCC.getThreshold();
				} catch( IllegalArgumentException e ) {
					if( exceptionIfNotComputeable ) {
						throw e;
					} else {
						cc = Double.NaN;
						t = Double.NaN;
					}
				}
				list.add( new NumericalResult( Measure.MaximumCorrelationCoefficient.getNameString(),
						Measure.MaximumCorrelationCoefficient.getCommentString(), cc ) );
				list.add( new NumericalResult( "Threshold for maximum correlation coefficient",
						"The threshold for the maximum correlation coefficient of all possible correlation coefficients", t ) );
				numSelected--;
			}
			if( params.isSelected( Measure.PartialROCCurve ) ) {
				RangeParameter specs = (RangeParameter)( (ParameterSet)params.getValueFor( Measure.PartialROCCurve.getNameString() ) ).getParameterAt( 0 );
				double[][] val = ScoreBasedPerformanceMeasureDefinitions.getPartialROC( sortedScores[0], sortedScores[1], specs );
				specs.resetToFirst();
				int i = 0;
				do {
					list.add( new NumericalResult( "" + val[0][i], "The sensitivity for a specificity of " + val[0][i], val[1][i] ) );
					list.add( new NumericalResult( "Threshold for a specificity of " + val[0][i],
							"The threshold for the sensitivity for a specificity of " + val[0][i],
							val[2][i] ) );
					i++;
				} while( specs.next() );
				numSelected--;
			}
			if( params.isSelected( Measure.ReceiverOperatingCharacteristicCurve ) ) {
				if( list2 == null ) {
					list2 = new ArrayList<double[]>( s[0].getNumberOfElements() / 2 );
					ScoreBasedPerformanceMeasureDefinitions.getAUC_ROC( sortedScores[0], sortedScores[1], list2 );
				}
				list.add( new DoubleTableResult( Measure.ReceiverOperatingCharacteristicCurve.getNameString(),
						Measure.ReceiverOperatingCharacteristicCurve.getCommentString(),
						list2 ) );
				numSelected--;
			}
			if( params.isSelected( Measure.PrecisionRecallCurve ) ) {
				if( list3 == null ) {
					list3 = new ArrayList<double[]>( s[0].getNumberOfElements() / 2 );
					ScoreBasedPerformanceMeasureDefinitions.getAUC_PR( sortedScores[0], sortedScores[1], list3 );
				}
				list.add( new DoubleTableResult( Measure.PrecisionRecallCurve.getNameString(),
						Measure.PrecisionRecallCurve.getCommentString(),
						list3 ) );
				numSelected--;
			}

			if( exceptionIfNotComputeable && numSelected > 0 ) {
				throw new IllegalArgumentException( "There are measures that could not be evaluate with this classifier (" + this.getClass() + ")." );
			}

			return list;
		}
	}

	/**
	 * Returns the specific class weights of a
	 * {@link AbstractScoreBasedClassifier}.
	 * 
	 * @return the class weights of the classifier
	 */
	public double[] getClassWeights() {
		return classWeights.clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#getNumberOfClasses()
	 */
	@Override
	public int getNumberOfClasses() {
		return classWeights.length;
	}

	/**
	 * This method returns the score for a given {@link Sequence} and a given
	 * class.
	 * 
	 * @param seq
	 *            the given {@link Sequence}
	 * @param i
	 *            the index of the class
	 * 
	 * @return the score for a given {@link Sequence} and a given class
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractScoreBasedClassifier#getScore(Sequence, int, boolean)
	 */
	public double getScore( Sequence seq, int i ) throws Exception {
		return getScore( seq, i, true );
	}

	/**
	 * Sets new class weights. <br>
	 * The logarithmic probabilities of an item <code>i</code> given class 0 to
	 * class <code>n</code> are computed to classify this item into a class. The
	 * class weights are added to each of these logarithmic probabilities. As
	 * higher (relational) the class weight of class <code>j</code>, as more
	 * probable it becomes, that any item is classified into this class. <br>
	 * Class weights do not have to be logarithmic probabilities. If
	 * {@latex.inline $\\sum_{j=0}^n \\exp(classWeight_j) \\stackrel{!}{=}1$},
	 * the class weights may be interpreted as logarithmic class-a-priori-probabilities.
	 * 
	 * @param add
	 *            indicates if the class weights are added to the current class
	 *            weights
	 * @param weights
	 *            the array of weights, for each class the weight that is added
	 *            in a classification
	 * 
	 * @throws ClassDimensionException
	 *             if something is wrong with the number of classes
	 */
	public final void setClassWeights( boolean add, double... weights ) throws ClassDimensionException {
		int c = getNumberOfClasses();
		if( weights == null || c != weights.length ) {
			throw new ClassDimensionException();
		}
		setClassWeights( add, weights, 0 );
	}
	
	/**
	 * Sets new class weights.
	 * 
	 * <br><br>
	 * 
	 * Only for internal use.
	 * 
	 * @param add
	 *            indicates if the class weights are added to the current class
	 *            weights
	 * @param weights
	 *            an array of weights that might have more entries than the classifier has classes
	 * @param start
	 * 			  the start index
	 * 
	 * @see #setClassWeights(boolean, double...)
	 */
	protected final void setClassWeights( boolean add, double[] weights, int start ) {
		if( add ) {
			for( int i = 0; i < classWeights.length; i++ ) {
				classWeights[i] += weights[start+i];
			}
		} else {
			for( int i = 0; i < classWeights.length; i++ ) {
				classWeights[i] = weights[start+i];
			}
		}
	}

	/**
	 * Sets a new threshold for 2-class-classifiers.<br>
	 * Only available if this {@link AbstractScoreBasedClassifier} distinguishes
	 * between 2 classes 0 and 1. In this case, <code>t</code> will be interpreted as
	 * {@latex.inline $\\log\\left(\\frac{P(class1)}{P(class0)}\\right)$}. A large
	 * <code>t</code> (greater than 0) makes the classifier to decide more often for
	 * class 1. A small <code>t</code> (smaller than 0) makes the classifier to
	 * decide more often for class 0.
	 * 
	 * @param add
	 *            indicates if the class weights are added to the current class
	 *            weights
	 * @param t
	 *            the new threshold
	 * 
	 * @throws OperationNotSupportedException
	 *             if the classifier is no 2-class-classifier
	 */
	public final void setThresholdClassWeights( boolean add, double t ) throws OperationNotSupportedException {
		int c = getNumberOfClasses();
		if( c != 2 ) {
			throw new OperationNotSupportedException();
		}
		if( classWeights == null ) {
			classWeights = new double[2];
		}
		// logarithmic bound
		// => t = log( p(1)/p(0) ) = log( (1-p)/p )
		// <=> exp(t) = (1-p)/p
		// <=> exp(t)p = (1-p)
		// <=> (exp(t) + 1)p = 1
		// <=> p = 1/(exp(t)+1)
		// <=> log(p) = log( 1/(exp(t)+1) )
		// => log(1-p) = t + log(p)
		double logP = -Math.log1p( Math.exp( t ) );
		if( add ) {
			classWeights[0] += logP;
			classWeights[1] += t + logP;
		} else {
			classWeights[0] = logP;
			classWeights[1] = t + logP;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#test(de.jstacs.data.Sample[])
	 */
	@Override
	public ConfusionMatrix test( Sample... testData ) throws Exception {
		if( testData.length != getNumberOfClasses() ) {
			throw new ClassDimensionException();
		}

		ConfusionMatrix matrix = new ConfusionMatrix( testData.length );
		ElementEnumerator ei;
		for( int i = 0; i < testData.length; i++ ) {
			if( testData[i] != null ) {
				check( testData[i] );
				ei = new ElementEnumerator( testData[i] );
				while( ei.hasMoreElements() ) {
					matrix.add( i, classify( ei.nextElement(), false ) );
				}
			}
		}
		return matrix;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = new StringBuffer( 300 );
		XMLParser.appendObjectWithTags( xml, classWeights, "classWeights" );
		return xml;
	}

	/**
	 * This method checks if the given {@link Sample} can be used.
	 * 
	 * @param s
	 *            the {@link Sample} to be checked
	 * 
	 * @throws NotTrainedException
	 *             if the classifier is not trained
	 * @throws IllegalArgumentException
	 *             if something is wrong with the {@link Sample} <code>s</code>
	 * 
	 * @see AbstractClassifier#setNewAlphabetContainerInstance(AlphabetContainer)
	 */
	protected void check( Sample s ) throws NotTrainedException, IllegalArgumentException {
		if( !isTrained() ) {
			throw new NotTrainedException( "The classifier is not trained yet." );
		}
		int length = getLength();
		if( length != 0 && s.getElementLength() != length ) {
			throw new IllegalArgumentException( "The sequences have not the correct length." );
		}
		if( !setNewAlphabetContainerInstance( s.getAlphabetContainer() ) ) {
			throw new IllegalArgumentException( "The sequences are not defined over the correct alphabets." );
		}
	}

	/**
	 * This method checks if the given {@link Sequence} can be used.
	 * 
	 * @param seq
	 *            the {@link Sequence} to be checked
	 * 
	 * @throws NotTrainedException
	 *             if the classifier is not trained
	 * @throws IllegalArgumentException
	 *             if something is wrong with the {@link Sequence}
	 *             <code>seq</code>
	 */
	protected void check( Sequence seq ) throws NotTrainedException, IllegalArgumentException {
		if( !isTrained() ) {
			throw new NotTrainedException( "The classifier is not trained yet." );
		}
		int length = getLength();
		if( length != 0 && seq.getLength() != length ) {
			throw new IllegalArgumentException( "The sequence has not the correct length." );
		}
		if( !getAlphabetContainer().checkConsistency( seq.getAlphabetContainer() ) ) {
			throw new IllegalArgumentException( "The sequence is not defined over the correct alphabets." );
		}
	}

	/**
	 * This method classifies a {@link Sequence}. It enables you to check the
	 * constraints (<code>alphabets</code>, <code>length</code>,
	 * {@link #isTrained()} ).
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param check
	 *            indicates if the constraints will be checked
	 * 
	 * @return the index of the class the {@link Sequence} is assigned to
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractScoreBasedClassifier#check(Sequence)
	 */
	protected byte classify( Sequence seq, boolean check ) throws Exception {
		if( check ) {
			check( seq );
		}
		byte clazz = 0;
		double max = getScore( seq, clazz, false ), current;
		for( byte i = 1; i < getNumberOfClasses(); i++ ) {
			current = getScore( seq, i, false );
			if( current > max ) {
				max = current;
				clazz = i;
			}
		}
		return clazz;
	}

	/**
	 * This method creates new class weights. Each class weight has the same
	 * value <code>val</code>. So the class weights do not have any influence on
	 * the classification.
	 * 
	 * @param classes
	 *            the number of different classes
	 * @param val
	 *            the value that is used for all classes
	 * 
	 * @throws IllegalArgumentException
	 *             if the number of classes is below 2
	 */
	protected void createDefaultClassWeights( int classes, double val ) throws IllegalArgumentException {
		if( classes < 2 ) {
			throw new IllegalArgumentException();
		}
		classWeights = new double[classes];
		Arrays.fill( classWeights, val );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		classWeights = XMLParser.extractObjectForTags( xml, "classWeights", double[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

	/**
	 * Returns the class weight for the class with a given <code>index</code>.
	 * 
	 * @param index
	 *            the given index of the class
	 * 
	 * @return the class weight
	 */
	protected double getClassWeight( int index ) {
		return classWeights[index];
	}

	/**
	 * This method returns the score for a given {@link Sequence} and a given
	 * class.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param i
	 *            the index of the class
	 * @param check
	 *            the switch to decide whether to check
	 *            {@link AlphabetContainer} and the length of the
	 *            {@link Sequence} or not
	 * 
	 * @return the score for a given {@link Sequence} and a given class
	 * 
	 * @throws NotTrainedException
	 *             if the classifier is not trained
	 * @throws IllegalArgumentException
	 *             if something is wrong with the {@link Sequence}
	 *             <code>seq</code>
	 * @throws Exception
	 *             if something went wrong
	 */
	protected abstract double getScore( Sequence seq, int i, boolean check ) throws IllegalArgumentException,
			NotTrainedException,
			Exception;

	/**
	 * This method returns the scores of the classifier for any {@link Sequence}
	 * in the {@link Sample}. The scores are stored in the array according to
	 * the index of the {@link Sequence} in the {@link Sample}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * <b>Only for 2-class-classifiers.</b>
	 * 
	 * @param s
	 *            the {@link Sample}
	 * 
	 * @return the array of scores
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public double[] getScores( Sample s ) throws Exception {
		if( classWeights.length != 2 ) {
			throw new OperationNotSupportedException( "This method is only for 2-class-classifiers." );
		}
		if( s == null ) {
			return new double[0];
		}
		check( s );
		double[] score = new double[s.getNumberOfElements()];
		ElementEnumerator ei = new ElementEnumerator( s );
		Sequence seq;
		for( int i = 0; i < score.length; i++ ) {
			seq = ei.nextElement();
			// for probabilistical models this should be:
			// ln p(x,fg|params) - ln p(x,bg|params)
			score[i] = getScore( seq, 0, false ) - getScore( seq, 1, false );
			if( Double.isNaN( score[i] ) ) {
				throw new IllegalArgumentException( "Could not classify sequence " + i
													+ ": "
													+ seq
													+ "\nfg: "
													+ getScore( seq, 0, false )
													+ "\nbg: "
													+ getScore( seq, 1, false ) );
			}
		}
		return score;
	}

	/**
	 * Returns the p-value for a {@link Sequence} <code>candidate</code> with
	 * respect to a given background {@link Sample}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The p-value is the percentage of background {@link Sequence}s that have a
	 * score that is at least as high as the score for the
	 * <code>candidate</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * It is not recommended to use this method in a <code>for</code>-loop. In
	 * such cases one should use the method that works on two {@link Sample}s.
	 * 
	 * @param candidate
	 *            the candidate {@link Sequence}
	 * @param bg
	 *            the background {@link Sample}
	 * 
	 * @return the p-value for the {@link Sequence} <code>candidate</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractScoreBasedClassifier#getPValue(Sample, Sample)
	 * @see PValueComputation#getPValue(double[], double)
	 */
	public double getPValue( Sequence candidate, Sample bg ) throws Exception {
		double[] scores = createStatistic( bg );
		return PValueComputation.getPValue( scores, getScore( candidate, 0 ) - getScore( candidate, 1 ) );
	}

	/**
	 * Returns the p-values for all {@link Sequence}s in the {@link Sample}
	 * <code>candidates</code> with respect to a given background {@link Sample}
	 * .
	 * 
	 * <br>
	 * <br>
	 * 
	 * The p-value is the percentage of background {@link Sequence}s that have a
	 * score that is at least as high as the score for a {@link Sequence} in
	 * <code>candidates</code>.
	 * 
	 * @param candidates
	 *            the {@link Sample} with candidate sequences
	 * @param bg
	 *            the background sample
	 * 
	 * @return the p-values for all sequences in <code>candidates</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractScoreBasedClassifier#getPValue(Sequence, Sample)
	 * @see PValueComputation#getPValue(double[], double)
	 */
	public double[] getPValue( Sample candidates, Sample bg ) throws Exception {
		double[] scores = createStatistic( bg ), pVal = new double[candidates.getNumberOfElements()];
		Sequence candidate;
		for( int i = 0; i < pVal.length; i++ ) {
			candidate = candidates.getElementAt( i );
			pVal[i] = PValueComputation.getPValue( scores, getScore( candidate, 0 ) - getScore( candidate, 1 ) );
		}
		return pVal;
	}

	// just put the scores in the array -> sort
	private double[] createStatistic( Sample bg ) throws Exception {
		double[] scores = getScores( bg );
		Arrays.sort( scores );
		// System.out.println( scores[0] + " " + scores[scores.length-1] );
		return scores;
	}

	private double[][] getSortedScore( Sample[] s ) throws Exception {
		double[][] scores = new double[][]{ getScores( s[0] ), getScores( s[1] ) };
		Arrays.sort( scores[0] );
		Arrays.sort( scores[1] );
		//System.out.println( "class 0: " + scores[0][0] + " to " + scores[0][scores[0].length-1] );
		//System.out.println( "class 1: " + scores[1][0] + " to " + scores[1][scores[1].length-1] );
		return scores;
	}

	/**
	 * This class is for {@link Result}s given as a table of <code>double</code>
	 * s.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see Result
	 */
	public static class DoubleTableResult extends Result {

		private double[][] content;

		private static final String XML_TAG = "DoubleTableResult";

		/**
		 * This is the default constructor that creates an instance based on the results
		 * given in <code>list</code>
		 * 
		 * @param name the name of the result
		 * @param comment the comment for the result
		 * @param list the list of values
		 */
		public DoubleTableResult( String name, String comment, AbstractList<double[]> list ) {
			super( name, comment, DataType.LIST );
			content = new double[list.size()][];
			for( int i = 0; i < content.length; i++ ) {
				content[i] = list.get( i ).clone();
			}
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link DoubleTableResult} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link DoubleTableResult} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see Result#Result(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public DoubleTableResult( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.results.Result#fromXML(java.lang.StringBuffer)
		 */
		@Override
		protected void fromXML( StringBuffer representation ) throws NonParsableException {
			StringBuffer xml = XMLParser.extractForTag( representation, XML_TAG );
			extractMainInfo( xml );
			content = XMLParser.extractObjectForTags( xml, "content", double[][].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		}

		/**
		 * Return the line with a given <code>index</code> from the table.
		 * 
		 * @param index
		 *            the given index
		 * 
		 * @return the line with the given <code>index</code>
		 */
		public double[] getLine( int index ) {
			return content[index].clone();
		}

		/**
		 * Returns the number of lines in this table.
		 * 
		 * @return the number of lines in this table
		 */
		public int getNumberOfLines() {
			return content.length;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString() {
			return "[table] \t " + name + " \t(" + comment + ")";
		}

		/* (non-Javadoc)
		 * @see de.jstacs.results.Result#getResult()
		 */
		@Override
		public double[][] getResult() {
			double[][] res = new double[content.length][];
			for( int i = 0; i < res.length; i++ ) {
				res[i] = content[i].clone();
			}
			return res;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.Storable#toXML()
		 */
		public StringBuffer toXML() {
			int l = content.length == 0 ? 0 : content[0].length;
			StringBuffer xml = new StringBuffer( 500 + content.length * l * 10 );
			appendMainInfo( xml );
			XMLParser.appendObjectWithTags( xml, content, "content" );
			XMLParser.addTags( xml, XML_TAG );
			return xml;
		}

		/**
		 * This method plots an array of {@link DoubleTableResult}s in one
		 * image.
		 * 
		 * @param dtr
		 *            the array of {@link DoubleTableResult}s
		 * @param e
		 *            the R environment
		 * 
		 * @return the image as an {@link ImageResult}
		 * 
		 * @throws Exception
		 *             if something went wrong
		 * 
		 * @see ImageResult#ImageResult(String, String, java.awt.image.BufferedImage)
		 */
		public static final ImageResult plot( REnvironment e, DoubleTableResult... dtr ) throws Exception {
			String opt = dtr[0].name;
			int i = 1;
			while( i < dtr.length && dtr[i].name.equalsIgnoreCase( opt ) ) {
				i++;
			}
			if( i != dtr.length ) {
				opt = null;
			}
			return new ImageResult( opt, "This plot shows the " + opt + ".", e.plot( getPlotCommands( e, opt, dtr ).toString() ) );
		}

		/**
		 * This method copies the data to the server side and creates a
		 * {@link StringBuffer} containing the plot commands.
		 * 
		 * @param dtr
		 *            the array of {@link DoubleTableResult}s
		 * @param e
		 *            the R environment
		 * @param plotOptions
		 *            <ol>
		 *            <li>recommended for {@link Measure#ReceiverOperatingCharacteristicCurve}
		 *            or {@link Measure#PrecisionRecallCurve}:
		 *            {@link Measure#getNameString()}</li>
		 *            <li>any String that can be parsed to R plot options</li>
		 *            <li>
		 * 
		 * @return a {@link StringBuffer} containing the plot commands
		 * 
		 * @throws Exception
		 *             if something went wrong
		 */
		public static final StringBuffer getPlotCommands( REnvironment e, String plotOptions, DoubleTableResult... dtr ) throws Exception {
			return getPlotCommands( e, plotOptions, (String[]) null, dtr );
		}
		
		/**
		 * This method copies the data to the server side and creates a
		 * {@link StringBuffer} containing the plot commands.
		 * 
		 * @param dtr
		 *            the array of {@link DoubleTableResult}s
		 * @param e
		 *            the R environment
		 * @param plotOptions
		 *            <ol>
		 *            <li>recommended for {@link Measure#ReceiverOperatingCharacteristicCurve}
		 *            or {@link Measure#PrecisionRecallCurve}:
		 *            {@link Measure#getNameString()}</li>
		 *            <li>any String that can be parsed to R plot options</li>
		 *            <li>
		 * @param colors array of colors for the dtrs
		 * 
		 * @return a {@link StringBuffer} containing the plot commands
		 * 
		 * @throws Exception
		 *             if something went wrong
		 */
		public static final StringBuffer getPlotCommands( REnvironment e, String plotOptions, int[] colors, DoubleTableResult... dtr ) throws Exception {
			String[] col = new String[colors.length];
			for( int i = 0; i < col.length; i++ ) {
				col[i] = "" + colors[i];
			}
			return getPlotCommands( e, plotOptions, col, dtr );
		}
		
		/**
		 * This method copies the data to the server side and creates a
		 * {@link StringBuffer} containing the plot commands.
		 * 
		 * @param dtr
		 *            the array of {@link DoubleTableResult}s
		 * @param e
		 *            the R environment
		 * @param plotOptions
		 *            <ol>
		 *            <li>recommended for {@link Measure#ReceiverOperatingCharacteristicCurve}
		 *            or {@link Measure#PrecisionRecallCurve}:
		 *            {@link Measure#getNameString()}</li>
		 *            <li>any String that can be parsed to R plot options</li>
		 *            <li>
		 * @param colors array of colors for the dtrs
		 * 
		 * @return a {@link StringBuffer} containing the plot commands
		 * 
		 * @throws Exception
		 *             if something went wrong
		 */
		public static final StringBuffer getPlotCommands( REnvironment e, String plotOptions, String[] colors, DoubleTableResult... dtr ) throws Exception {
			int i = 0;
			while( i < dtr.length ) {
				// this is not really secure,
				// but as long as nobody changes something in REnviroment this is good and fast
				e.createMatrix( "dtr" + i, dtr[i].content );
				i++;
			}

			if( plotOptions == null ) {
				plotOptions = dtr[0].name == null ? "" : dtr[0].name;
			}
			if( plotOptions.equals( Measure.ReceiverOperatingCharacteristicCurve.getNameString() ) ) {
				plotOptions = ", xlim=c(0, 1), ylim=c(0, 1), xlab=\"false positive rate\", ylab=\"sensitivity\", main=\"ROC curve\", lwd=3";
			} else if( plotOptions.equals( Measure.PrecisionRecallCurve.getNameString() ) ) {
				plotOptions = ", xlim=c(0, 1), ylim=c(0, 1), xlab=\"recall\", ylab=\"precision\", main=\"PR curve\", lwd=3";
			} else {
				plotOptions = plotOptions.trim();
				if( plotOptions.charAt( 0 ) != ',' ) {
					plotOptions = ", " + plotOptions;
				}
			}

			StringBuffer p = new StringBuffer( dtr.length * 200 );
			p.append( "plot( dtr0[,1], dtr0[,2], col="+((colors == null || colors.length==0) ? 1 : ("\"" + colors[0] + "\"") )+", type=\"l\"" + plotOptions + " );" );
			for( i = 1; i < dtr.length; ) {
				p.append( "\nlines( dtr" + i + "[,1], dtr" + i + "[,2], col=" + ((colors == null|| colors.length==0) ? ++i : ("\"" + colors[i++] + "\"") ) + ", lwd=3 );" );
			}
			return p;
		}
	}
}
