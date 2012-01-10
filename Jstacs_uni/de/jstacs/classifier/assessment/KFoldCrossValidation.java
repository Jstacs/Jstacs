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

package de.jstacs.classifier.assessment;

import java.util.LinkedList;

import de.jstacs.WrongAlphabetException;
import de.jstacs.classifier.AbstractClassifier;
import de.jstacs.classifier.ClassDimensionException;
import de.jstacs.classifier.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.MeanResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.trainableStatisticalModels.TrainableStatisticalModel;
import de.jstacs.utils.ProgressUpdater;

/**
 * This class implements a k-fold crossvalidation. A k-fold crossvalidation
 * assesses classifiers using the following methodology. The user supplies
 * datasets (one for each class the classifiers are capable to distinguish). The
 * data is randomly, mutually exclusive partitioned into k parts. Each of those
 * parts is used once as a test dataset while the remaining k-1 parts are used
 * as train dataset. In each of the k iterations, the train datasets are used to
 * train the classifier and the test datasets are used to assess the
 * classifier's performance to predict the elements therein. Additional the user
 * has to define which assessment measures should be used.
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com), Jens Keilwagen
 */
public class KFoldCrossValidation extends ClassifierAssessment {

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
	 * Creates a new {@link KFoldCrossValidation} from an array of
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
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param checkAlphabetConsistencyAndLength
	 *            indicates if alphabets and lengths shall be checked for
	 *            consistency
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
	protected KFoldCrossValidation( AbstractClassifier[] aCs, TrainableStatisticalModel[][] aMs, boolean buildClassifiersByCrossProduct,
									boolean checkAlphabetConsistencyAndLength ) throws IllegalArgumentException, WrongAlphabetException,
																				CloneNotSupportedException, ClassDimensionException {
		super( aCs, aMs, buildClassifiersByCrossProduct, checkAlphabetConsistencyAndLength );
	}

	/**
	 * Creates a new {@link KFoldCrossValidation} from a set of
	 * {@link AbstractClassifier}s.
	 * 
	 * @param aCs
	 *            contains the classifiers to be assessed.<br>
	 *            If model based classifiers are trained, the order of models in
	 *            classifiers determines, which model will be trained using
	 *            which sample in method <code>assess( ... )</code>.<br>
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
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier...)
	 */
	public KFoldCrossValidation( AbstractClassifier... aCs ) throws IllegalArgumentException, WrongAlphabetException,
															CloneNotSupportedException, ClassDimensionException {
		super( aCs );
	}

	/**
	 * Creates a new {@link KFoldCrossValidation} from a set of {@link TrainableStatisticalModel}s.
	 * The argument <code>buildClassifiersByCrossProduct</code> determines how
	 * these {@link TrainableStatisticalModel}s are combined to classifiers.
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
	 *            If models are trained directly (during assessment), the order
	 *            of given models during initiation of this assessment object
	 *            determines, which sample will be used for training which
	 *            model. In general the first model will be trained using the
	 *            first sample in <code>s</code>... . <br>
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
	public KFoldCrossValidation( boolean buildClassifiersByCrossProduct, TrainableStatisticalModel[]... aMs ) throws IllegalArgumentException,
																							WrongAlphabetException,
																							CloneNotSupportedException,
																							ClassDimensionException {
		super( buildClassifiersByCrossProduct, aMs );
	}

	/**
	 * This constructor allows to assess a collection of given
	 * {@link AbstractClassifier}s and those constructed using the given
	 * {@link de.jstacs.trainableStatisticalModels.TrainableStatisticalModel}s by a {@link KFoldCrossValidation}
	 * . <br>
	 * 
	 * @param aCs
	 *            contains some {@link AbstractClassifier}s that should be
	 *            assessed in addition to the {@link AbstractClassifier}
	 *            constructed using the given
	 *            {@link de.jstacs.trainableStatisticalModels.TrainableStatisticalModel}s
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
	 *            determines, which sample will be used for training which
	 *            model. In general the first model will be trained using the
	 *            first sample in <code>s</code>... . <br>
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
	public KFoldCrossValidation( AbstractClassifier[] aCs, boolean buildClassifiersByCrossProduct, TrainableStatisticalModel[]... aMs )
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
	 * Evaluates a classifier.
	 * 
	 * @param mp
	 *            defines which performance measures are used to assess
	 *            classifiers
	 * @param pU
	 *            the progress updater which shows the progress of the k-fold
	 *            crossvalidation
	 * @param s
	 *            contains the data to be used for assessment. The order of
	 *            samples is important. <br>
	 *            If model based classifiers are trained, the order of the
	 *            models in the classifiers determines, which model will be
	 *            trained using which sample. The first model in the classifier
	 *            will be trained using the first sample in <code>s</code>. If
	 *            the models are trained directly, the order of the given models
	 *            during initiation of this assessment object determines, which
	 *            sample will be used for training which model. In general the
	 *            first model will be trained using the first sample in
	 *            <code>s</code>... . <br>
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
	 * @param assessPS
	 *            contains parameters for a run of this
	 *            {@link KFoldCrossValidation}. Must be of type
	 *            {@link KFoldCrossValidationAssessParameterSet}.
	 * @throws IllegalArgumentException
	 *             if the given <code>assessPS</code> is not of type
	 *             {@link KFoldCrossValidationAssessParameterSet}
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see ClassifierAssessment#evaluateClassifier(NumericalPerformanceMeasureParameterSet, ClassifierAssessmentAssessParameterSet, DataSet[], ProgressUpdater)
	 * 
	 */
	@Override
	protected void evaluateClassifier( NumericalPerformanceMeasureParameterSet mp, ClassifierAssessmentAssessParameterSet assessPS, DataSet[] s, ProgressUpdater pU ) throws IllegalArgumentException,
			Exception {

		KFoldCrossValidationAssessParameterSet tempAssessPS = null;

		try {
			tempAssessPS = (KFoldCrossValidationAssessParameterSet)assessPS;
		} catch ( ClassCastException e ) {
			throw new IllegalArgumentException( "Given AssessParameterSet is not of type " + "KFoldCVAssessParameterSet." );
		}

		PartitionMethod splitMethod = tempAssessPS.getDataSplitMethod();
		int k = tempAssessPS.getK();

		DataSet[][] sInParts = new DataSet[s.length][];
		try {
			for( int i = 0; i < sInParts.length;
			//im bg-Fall sollte das Sample aus den unterschiedlichen grossen Sequenzen bestehen
			sInParts[i] = s[i++].partition( k, splitMethod ) );
		} catch ( EmptyDataSetException e ) {
			throw new IllegalArgumentException( "Given Sample s seems to contain to few elements " + "for a "
												+ k
												+ "-fold crossvalidation since at least one empty subset occured "
												+ "during splitting given data into "
												+ k
												+ " non-overlapping parts." );
		}
		evaluate( mp, tempAssessPS, pU, sInParts );
	}

	/**
	 * This method implements the core of the k-fold crossvalidation.
	 * 
	 * @param mp
	 *            defines which performance measures are used to assess
	 *            classifiers
	 * @param caaps
	 *            contains the defined element length and choice whether to
	 *            throw an {@link Exception} if a measure could not be computed
	 * @param pU
	 *            the progress updater which shows the progress of the k-fold
	 *            crossvalidation
	 * @param splitData
	 *            the previously split data; <code>splitData[i]</code> contains
	 *            the splits for class <code>i</code>; therefore the length of
	 *            each subarray <code>splitData[i]</code> has to be identical
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	private void evaluate( NumericalPerformanceMeasureParameterSet mp, ClassifierAssessmentAssessParameterSet caaps, ProgressUpdater pU, DataSet[]... splitData ) throws Exception {
		int subSeqL = caaps.getElementLength();
		boolean exceptionIfMPNotComputable = caaps.getExceptionIfMPNotComputable();

		int clazz = splitData.length, k = splitData[0].length, j = 1;
		DataSet[][] sTrainTestClassWise = new DataSet[2][clazz];
		boolean[] tempBool = new boolean[k];
		java.util.Arrays.fill( tempBool, true );
		while( j < splitData.length && splitData[j].length == k ) {
			j++;
		}
		if( j != splitData.length ) {
			throw new IllegalArgumentException( "Please check the number of predefined splits per class. Compare class 0 with class " + j );
		}
		pU.setMax( k );
		for( int i = 0; i < k; i++ ) {
			//prepare data
			tempBool[i] = false;
			for( j = 0; j < clazz; j++ ) {

				//train-data
				//im bg-Fall beinhalten die Teil-Mengen die gesamten unterschiedlichen Langen Sequenzen
				//das bedeutet, dass nach der Anzahl der Symbole zu urteilen, die Partitionen sehr
				//wahrscheinlich nicht gleich gross sind
				sTrainTestClassWise[0][j] = DataSet.union( splitData[j], tempBool );

				//test-data
				//muss auf richtige subSeqL gesetzt werden, damit auch die fg-Modelle
				//auf diese angewendet werden koennen.
				sTrainTestClassWise[1][j] = new DataSet( splitData[j][i], subSeqL );
			}
			tempBool[i] = true;

			train( sTrainTestClassWise[0] );

			test( mp, exceptionIfMPNotComputable, sTrainTestClassWise[1] );
			pU.setValue( i + 1 );
		}
	}

	/**
	 * This method implements a k-fold crossvalidation on previously split data.
	 * This might be useful if you like to compare the results of your
	 * classifier(s) with those of a previous study in a paper or manuscript.
	 * 
	 * @param mp
	 *            defines which performance measures are used to assess
	 *            classifiers
	 * @param caaps
	 *            contains the defined element length and choice whether an
	 *            exception should be thrown if a measure could not be computed
	 * @param pU
	 *            the progress updater which shows the progress of the k-fold
	 *            crossvalidation
	 * @param splitData
	 *            the previously split data; <code>splitData[i]</code> contains
	 *            the splits for class <code>i</code>; therefore the length of
	 *            each subarray <code>splitData[i]</code> has to to be identical
	 * 
	 * @return the result packed in a {@link ListResult}
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */

	public ListResult assessWithPredefinedSplits( NumericalPerformanceMeasureParameterSet mp, ClassifierAssessmentAssessParameterSet caaps, ProgressUpdater pU,
			DataSet[]... splitData ) throws Exception {
		int clazz = myAbstractClassifier[0].getNumberOfClasses();
		if( splitData.length != clazz ) {
			throw new IllegalArgumentException( "The number of classes in the data array and the classifier differs." );
		}

		this.myTempMeanResultSets = new MeanResultSet[this.myAbstractClassifier.length];
		for( int i = 0; i < this.myAbstractClassifier.length; i++ ) {
			myTempMeanResultSets[i] = new MeanResultSet( myAbstractClassifier[i].getClassifierAnnotation() );
		}

		evaluate( mp, caaps, pU, splitData );

		LinkedList<Result> annotation = new LinkedList<Result>();
		annotation.add( new CategoricalResult( "kind of assessment", "a description or name of the assessment", getNameOfAssessment() ) );
		annotation.addAll( caaps.getAnnotation() );
		StringBuffer sb = new StringBuffer( 1000 );
		sb.append( "[" + DataSet.getAnnotation( splitData[0] ) );
		for( int i = 1; i < splitData.length; i++ ) {
			sb.append( ", " + DataSet.getAnnotation( splitData[i] ) );
		}
		sb.append( "]" );
		annotation.add( new CategoricalResult( "samples", "annotation of used samples", "predefined splits: " + sb ) );

		return new ListResult( this.getNameOfAssessment(),
				"the results of a " + this.getNameOfAssessment() + " of predefined splits",
				new ResultSet( annotation ),
				this.myTempMeanResultSets );
	}
}
