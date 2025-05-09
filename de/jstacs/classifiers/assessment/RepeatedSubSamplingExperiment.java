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

package de.jstacs.classifiers.assessment;

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.ClassDimensionException;
import de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ProgressUpdater;

/**
 * This class implements a repeated subsampling experiment. A repeated
 * subsampling experiment uses the following procedure to assess classifiers.<br>
 * The user supplies data sets for each class the classifiers are capable to
 * distinguish. In each step of the repeated subsampling these data sets are used to generate test and train
 * data sets by subsampling of these data sets. Hence, the subsampled train and test data set of 
 * each class do not have to be disjoint. The user defines how many elements
 * the subsampled test and train data sets should contain. After subampling the
 * train data sets are used to train the classifiers and the test data sets are
 * used to assess the performance of the classifiers to predict the elements
 * therein. Additionally the user defines how often these procedure should be
 * repeated and which assessment measures are used to assess the classifiers.
 * 
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com)
 * 
 */
public class RepeatedSubSamplingExperiment extends ClassifierAssessment<RepeatedSubSamplingAssessParameterSet> {

	//	**********************
	//	class variables
	//	**********************

	//	**********************
	//	class methods
	//	**********************

	//	**********************
	//	member variables
	//	**********************

	//	**********************
	//	constructors
	//	**********************	

	/**
	 * Creates a new {@link RepeatedSubSamplingExperiment} from an array of
	 * {@link AbstractClassifier}s and a two-dimensional array of {@link TrainableStatisticalModel}
	 * s, which are combined to additional classifiers. If
	 * <code>buildClassifiersByCrossProduct</code> is <code>true</code>, the
	 * cross-product of all {@link TrainableStatisticalModel}s in <code>aMs</code> is built to
	 * obtain these classifiers.
	 * 
	 * @param aCs
	 *            the predefined classifiers
	 * @param aMs
	 *            the {@link TrainableStatisticalModel}s that are used to build additional
	 *            classifiers
	 * @param buildClassifiersByCrossProduct
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of S are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param checkAlphabetConsistencyAndLength
	 *            indicates if alphabets and lengths shall be checked for
	 *            consistency
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if the classifiers use different alphabets
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      TrainableStatisticalModel[][], boolean, boolean)
	 */
	protected RepeatedSubSamplingExperiment( AbstractClassifier[] aCs, TrainableStatisticalModel[][] aMs, boolean buildClassifiersByCrossProduct,
												boolean checkAlphabetConsistencyAndLength ) throws IllegalArgumentException,
																							WrongAlphabetException,
																							CloneNotSupportedException,
																							ClassDimensionException {
		super( aCs, aMs, buildClassifiersByCrossProduct, checkAlphabetConsistencyAndLength );
	}

	/**
	 * Creates a new {@link RepeatedSubSamplingExperiment} from a set of
	 * {@link AbstractClassifier}s.
	 * 
	 * @param aCs
	 *            contains the classifiers to be assessed.<br>
	 *            If model based classifiers are trained, the order of the
	 *            models in the classifiers determines, which model will be
	 *            trained using which data set in method
	 *            <code>assess( ... )</code>.<br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate a assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link de.jstacs.data.AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier...)
	 */
	public RepeatedSubSamplingExperiment( AbstractClassifier... aCs ) throws IllegalArgumentException, WrongAlphabetException,
																		CloneNotSupportedException, ClassDimensionException {
		super( aCs );
	}

	/**
	 * Creates a new {@link RepeatedSubSamplingExperiment} from a set of
	 * {@link TrainableStatisticalModel}s. The argument <code>buildClassifiersByCrossProduct</code>
	 * determines how these {@link TrainableStatisticalModel}s are combined to classifiers.
	 * 
	 * @param buildClassifiersByCrossProduct
	 * <br>
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param aMs
	 * <br>
	 *            Contains the models in the following way (suppose a k-class
	 *            problem): the first dimension encodes the class (here it is
	 *            k), the second dimension (<code>aMs[i]</code>) contains the
	 *            models according to class <code>i</code>.<br>
	 *            If the models are trained directly (during assessment), the
	 *            order of given models during initiation of this assessment
	 *            object determines, which data set will be used for training
	 *            which model. In general the first model will be trained using
	 *            the first data set in <code>s</code>... . <br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate an assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link de.jstacs.data.AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(boolean, TrainableStatisticalModel[][])
	 */
	public RepeatedSubSamplingExperiment( boolean buildClassifiersByCrossProduct, TrainableStatisticalModel[]... aMs ) throws IllegalArgumentException,
																									WrongAlphabetException,
																									CloneNotSupportedException,
																									ClassDimensionException {
		super( buildClassifiersByCrossProduct, aMs );
	}

	/**
	 * This constructor allows to assess a collection of given
	 * {@link AbstractClassifier}s and those constructed using the given
	 * {@link de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel}s by a
	 * {@link RepeatedSubSamplingExperiment}. <br>
	 * 
	 * @param aCs
	 *            contains some {@link AbstractClassifier} that should be
	 *            assessed in addition to the {@link AbstractClassifier}
	 *            constructed using the given
	 *            {@link de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel}s
	 * @param buildClassifiersByCrossProduct
	 * <br>
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param aMs
	 * <br>
	 *            Contains the models in the following way (suppose a k-class
	 *            problem): the first dimension encodes the class (here it is
	 *            k), the second dimension (<code>aMs[i]</code>) contains the
	 *            models according to class <code>i</code>.<br>
	 *            If models are trained directly (during assessment), the order
	 *            of given models during initiation of this assessment object
	 *            determines, which data set will be used for training which
	 *            model. In general the first model will be trained using the
	 *            first data set in <code>s</code>... . <br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate an assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link de.jstacs.data.AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      boolean, TrainableStatisticalModel[][])
	 */
	public RepeatedSubSamplingExperiment( AbstractClassifier[] aCs, boolean buildClassifiersByCrossProduct, TrainableStatisticalModel[]... aMs )
																															throws IllegalArgumentException,
																															WrongAlphabetException,
																															CloneNotSupportedException,
																															ClassDimensionException {
		super( aCs, buildClassifiersByCrossProduct, aMs );
	}

	//	**********************
	//	member methods
	//	**********************

	/**
	 * Evaluates the classifier.
	 * 
	 * @param mp
	 *            defines which performance measures are used to assess
	 *            classifiers
	 * @param pU
	 *            allows to abort this assessment by setting
	 *            <code>pU.isCancelled()=true</code>. The last step of this
	 *            assessment is continued afterwards this assessment stops.
	 * @param s
	 *            contains the data to be used for assessment. The order of the
	 *            data sets is important. <br>
	 *            If model based classifiers are trained, the order of the
	 *            models in the classifiers determines, which model will be
	 *            trained using which data set. The first model in the classifier
	 *            will be trained using the first data set in <code>s</code>. If
	 *            the models are trained directly, the order of given models
	 *            during initiation of this assessment object determines, which
	 *            data set will be used for training which model. In general the
	 *            first model will be trained using the first data set in
	 *            <code>s</code>... . <br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate a assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * @param assessPS
	 *            contains parameters for a run of this
	 *            {@link RepeatedSubSamplingExperiment}. Must be of type
	 *            {@link RepeatedSubSamplingExperiment}.
	 * 
	 * @throws IllegalArgumentException
	 *             if the given <code>assessPS</code> is not of type
	 *             {@link RepeatedSubSamplingExperiment}
	 * @throws Exception
	 *             if something went wrong
	 */
	@Override
	protected void evaluateClassifier( NumericalPerformanceMeasureParameterSet mp, RepeatedSubSamplingAssessParameterSet assessPS, DataSet[] s, double[][] weights, ProgressUpdater pU ) throws IllegalArgumentException,
			Exception {

		int eL = assessPS.getElementLength();
		boolean exceptionIfMPNotComputable = assessPS.getExceptionIfMPNotComputable();
		int repeats = assessPS.getRepeats();
		double[] trainN = assessPS.getTrain_TestNumbers( true );
		double[] testN = assessPS.getTrain_TestNumbers( false );
		int subSeqL = assessPS.getElementLength();

		if( trainN.length != this.myAbstractClassifier[0].getNumberOfClasses() ) {
			throw new IllegalArgumentException( "Given RepeatedSubSamlpingAssessParameterSet contains " + "a invalid parameter trainDataNumbers. trainDataNumbers (int[], number of train-data items that should be subsampled "
												+ "for each class) must contain as much entries "
												+ "as classes the local classifers are able to distinguish." );
		}

		if( testN.length != this.myAbstractClassifier[0].getNumberOfClasses() ) {
			throw new IllegalArgumentException( "Given RepeatedSubSamlpingAssessParameterSet contains " + "a invalid parameter testDataNumbers. testDataNumbers (int[], number of test-data items that should be subsampled "
												+ "for each class) must contain as much entries "
												+ "as classes the local classifers are able to distinguish." );
		}

		DataSet[][] sTrainTestClassWise = new DataSet[2][s.length];
		double[][][] weightsTrainTestClassWise = new double[2][s.length][];

		pU.setMax( repeats );

		Pair<DataSet, double[]> p;
		for( int iteration = 0; iteration < repeats; iteration++ ) {
			for( int classes = 0; classes < s.length; classes++ ) {
				//train
				p = s[classes].subSampling(trainN[classes], weights[classes]);
				sTrainTestClassWise[0][classes] = p.getFirstElement();
				weightsTrainTestClassWise[0][classes] = p.getSecondElement();

				//test
				p = s[classes].subSampling(testN[classes], weights[classes]);
				p = p.getFirstElement().resize( p.getSecondElement(), subSeqL );
				sTrainTestClassWise[1][classes] = p.getFirstElement();
				weightsTrainTestClassWise[1][classes] = p.getSecondElement();				
			}

			train( sTrainTestClassWise[0], weightsTrainTestClassWise[0] );
			test( mp, exceptionIfMPNotComputable, sTrainTestClassWise[1], weightsTrainTestClassWise[1] );

			pU.setValue( iteration + 1 );
			if( pU.isCancelled() ) {
				break;
			}
		}
	}

	public RepeatedSubSamplingAssessParameterSet getAssessParameterSet() throws Exception {
		return new RepeatedSubSamplingAssessParameterSet();
	}
}