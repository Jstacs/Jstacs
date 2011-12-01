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

import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.classifier.AbstractClassifier;
import de.jstacs.classifier.ClassDimensionException;
import de.jstacs.classifier.measures.MeasureParameters;
import de.jstacs.classifier.modelBased.ModelBasedClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.io.ArrayHandler;
import de.jstacs.models.Model;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.MeanResultSet;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.MeanResultSet.AdditionImpossibleException;
import de.jstacs.results.MeanResultSet.InconsistentResultNumberException;
import de.jstacs.utils.NullProgressUpdater;
import de.jstacs.utils.ProgressUpdater;

/**
 * Class defining an assessment of classifiers. <br>
 * It should be used as a superclass for specialized classifier assessments like
 * k-fold-crossvalidation or subsampling-crossvalidation.<br>
 * Several standard tasks like classifier or model management (testing,
 * training) are implemented. The method <code>assess( ... )</code> should be
 * used as standard method to start a classifier assessment. Subclasses have to
 * implement the method <code>evaluateClassifier( ... )</code>. This method
 * mainly has to execute the construction of test and training subsets of the
 * given data. These test and training subsets may be used by the methods
 * <code>test( ... )</code> and <code>train( ... )</code> which already are
 * implemented in a standard way.
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com), Jens Keilwagen
 */
public abstract class ClassifierAssessment {

	// ***********************
	// static variables
	// ***********************

	AbstractClassifier[] NULL_AC_ARRAY = new AbstractClassifier[0];

	Model[] NULL_AM_ARRAY = new Model[0];

	Model[][] NULL_AM_2DARRAY = new Model[0][0];

	// ***********************
	// static methods
	// ***********************

	/**
	 * Constructs an array of classifiers using the given models.
	 * 
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
	 *            constructed . In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param aMs
	 *            Contains the models in the following way (suppose a k-class
	 *            problem): the first dimension encodes the class (here it is
	 *            k), the second dimension (<code>aMs[i]</code>) contains the
	 *            models according to class <code>i</code>.
	 * 
	 * @return contains the constructed classifiers
	 */
	private static AbstractClassifier[] turnIntoClassifierArray( boolean buildClassifiersByCrossProduct, Model[][] aMs ) throws IllegalArgumentException,
			WrongAlphabetException,
			CloneNotSupportedException,
			ClassDimensionException {
		// evl. Ueberpruefung, ob alle Modelle auf gleichem AC arbeiten
		// Bestimmung der Laenge der Seq, die klassifiziert werden koennen, ect. pp. ....
		// erzeugt aus den aMs Klassifikatoren als Kreuzprodukt der aMs oder einfach spaltenweise

		if( aMs.length < 2 ) {
			throw new IllegalArgumentException( "Given classes of models is less than 2 but a classifier consists of at least 2 models." );
		}

		boolean allFirstDimsEqual = true;
		int tempDim = aMs[0].length;

		// determine length
		int l = 0;
		for( int counter1 = 0; counter1 < aMs.length; counter1++ ) {
			if( aMs[counter1].length == 0 ) {
				return new AbstractClassifier[0];
			}
			if( tempDim != aMs[counter1].length ) {
				allFirstDimsEqual = false;
			}

			for( int counter2 = 0; counter2 < aMs[counter1].length; counter2++ ) {
				if( aMs[counter1][counter2].getLength() != 0 ) {
					if( l == 0 ) {
						l = aMs[counter1][counter2].getLength();
					} else {
						if( aMs[counter1][counter2].getLength() != l ) {
							throw new IllegalArgumentException( "All models have to be able to use sequences of the same length " + "but at least 2 given models have different length (!=0)." );
						}
					}// else
				}// else
			}// for-counter2
		}// for-counter1

		// create classifiers
		Model[] m = new Model[aMs.length];
		ModelBasedAssessmentClassifier[] erg = null;

		if( buildClassifiersByCrossProduct ) {

			int anz = 1, counter2;
			int[] current = new int[aMs.length], max = new int[aMs.length];
			for( int counter1 = 0; counter1 < aMs.length; counter1++ ) {
				anz *= aMs[counter1].length;
				max[counter1] = aMs[counter1].length - 1;
			}

			erg = new ModelBasedAssessmentClassifier[anz];
			anz = m.length - 1;

			for( int counter1 = 0; counter1 < erg.length; counter1++ ) {

				for( counter2 = 0; counter2 < m.length; counter2++ ) {
					m[counter2] = aMs[counter2][current[counter2]];
				}

				erg[counter1] = new ModelBasedAssessmentClassifier( m );

				counter2 = 0;
				while( counter2 < anz && current[counter2] == max[counter2] ) {
					current[counter2++] = 0;
				}
				current[counter2]++;
			}
		}// if(buildClassifierByCrossProduct)
		else {
			if( !allFirstDimsEqual ) {
				throw new IllegalArgumentException( "Not all first dimensions of given aMs are equal but this is necessary if " + "classfiers should not be constructed by cross-product of given models." );
			}

			// check AlphabetConsistency in first dim
			AlphabetContainer abc = aMs[0][0].getAlphabetContainer();
			for( int counter1 = 1; counter1 < aMs[0].length; counter1++ ) {
				if( !abc.checkConsistency( aMs[0][counter1].getAlphabetContainer() ) ) {
					throw new WrongAlphabetException( "At least 2 models use different alphabets." );
				}
			}// for

			erg = new ModelBasedAssessmentClassifier[tempDim];

			for( int i = 0; i < erg.length; i++ ) {

				for( int j = 0; j < aMs.length; m[j] = aMs[j++][i] );

				erg[i] = new ModelBasedAssessmentClassifier( m );
			}

		}// else

		return erg;
	}

	// ************************
	// member variables
	// ************************

	/**
	 * This array contains the internal used classifiers.
	 */
	protected AbstractClassifier[] myAbstractClassifier;

	/**
	 * This array contains for each class the internal used models. This is
	 * helpful to allow to train some models only once while evaluating them in
	 * combination with other models.
	 */
	protected Model[][] myModel;

	/**
	 * The temporary result set.
	 */
	protected MeanResultSet[] myTempMeanResultSets;

	/**
	 * Skip last classifier.
	 */
	protected int skipLastClassifiersDuringClassifierTraining;

	//*************************
	//constructors
	//*************************

	/**
	 * Creates a new {@link ClassifierAssessment} from an array of
	 * {@link AbstractClassifier}s and a two-dimensional array of {@link Model}
	 * s, which are combined to additional classifiers. If
	 * <code>buildClassifiersByCrossProduct</code> is <code>true</code>, the
	 * cross product of all {@link Model}s in <code>aMs</code> is built to
	 * obtain these classifiers.
	 * 
	 * @param aCs
	 *            the predefined classifiers
	 * @param aMs
	 *            the {@link Model}s that are used to build additional
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
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 */
	protected ClassifierAssessment( AbstractClassifier[] aCs, Model[][] aMs, boolean buildClassifiersByCrossProduct,
									boolean checkAlphabetConsistencyAndLength ) throws IllegalArgumentException, WrongAlphabetException,
																				CloneNotSupportedException, ClassDimensionException {

		// construct AbstractClassifiers using given Models

		AbstractClassifier[] tempACs;
		if( aMs == null || aMs.length == 0 ) {
			this.myModel = NULL_AM_2DARRAY;
			tempACs = NULL_AC_ARRAY;
		} else {
			// clone given Models
			this.myModel = new Model[aMs.length][];
			for( int i = 0; i < this.myModel.length; i++ ) {
				this.myModel[i] = ArrayHandler.clone( aMs[i] );
			}
			tempACs = turnIntoClassifierArray( buildClassifiersByCrossProduct, this.myModel );
		}

		this.skipLastClassifiersDuringClassifierTraining = tempACs.length;

		//check and clone given AbstractClassifiers

		if( aCs == null ) {
			aCs = NULL_AC_ARRAY;
		}

		if( aCs.length > 0 && checkAlphabetConsistencyAndLength ) {

			int tempL = aCs[0].getLength();
			AlphabetContainer abc = aCs[0].getAlphabetContainer();
			for( int counter1 = 1; counter1 < aCs.length; counter1++ ) {
				if( !abc.checkConsistency( aCs[counter1].getAlphabetContainer() ) ) {
					throw new WrongAlphabetException( "At least 2 classifiers use different alphabets." );
				}

				if( aCs[counter1].getLength() != 0 ) {
					if( tempL == 0 ) {
						aCs[counter1].getLength();
					} else {
						if( aCs[counter1].getLength() != tempL ) {
							throw new IllegalArgumentException( "At least 2 classifiers have different length (!=0)." );
						}
					}// else

				}// if

			}// for

		}// if-checkAlphabetConsistencyAndLength

		this.myAbstractClassifier = new AbstractClassifier[aCs.length + tempACs.length];

		for( int i = 0; i < aCs.length; this.myAbstractClassifier[i] = aCs[i++].clone() );

		for( int i = 0; i < tempACs.length; this.myAbstractClassifier[i + aCs.length] = tempACs[i++] );

		//this.myParameterSet=null;
	}

	/**
	 * Creates a new {@link ClassifierAssessment} from a set of
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
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      Model[][], boolean, boolean)
	 */
	public ClassifierAssessment( AbstractClassifier... aCs ) throws IllegalArgumentException, WrongAlphabetException,
															CloneNotSupportedException, ClassDimensionException {
		this( aCs, null, false, true );
	}

	/**
	 * Creates a new {@link ClassifierAssessment} from a set of {@link Model}s.
	 * The argument <code>buildClassifiersByCrossProduct</code> determines how
	 * these {@link Model}s are combined to classifiers.
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
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      Model[][], boolean, boolean)
	 */
	public ClassifierAssessment( boolean buildClassifiersByCrossProduct, Model[]... aMs ) throws IllegalArgumentException,
																							WrongAlphabetException,
																							CloneNotSupportedException,
																							ClassDimensionException {
		this( null, aMs, buildClassifiersByCrossProduct, false );
	}

	/**
	 * This constructor allows to assess a collection of given
	 * {@link AbstractClassifier}s and, in addition, classifiers that will be
	 * constructed using the given {@link Model}s.
	 * 
	 * @param aCs
	 *            contains some {@link AbstractClassifier}s that should be
	 *            assessed in addition to the {@link AbstractClassifier}s
	 *            constructed using the given {@link Model}s
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
	 * @param aMs
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
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      Model[][], boolean, boolean)
	 */
	public ClassifierAssessment( AbstractClassifier[] aCs, boolean buildClassifiersByCrossProduct, Model[]... aMs )
																													throws IllegalArgumentException,
																													WrongAlphabetException,
																													CloneNotSupportedException,
																													ClassDimensionException {
		this( aCs, aMs, buildClassifiersByCrossProduct, true );
	}

	//***************************
	//member methods
	//***************************

	/**
	 * Assesses the contained classifiers.
	 * 
	 * @param s
	 *            contains the data to be used for assessment. The order of the
	 *            samples is important. <br>
	 *            If model based classifiers are trained, the order of the
	 *            models in the classifiers determines, which model will be
	 *            trained using which sample. The first model in the classifier
	 *            will be trained using the first sample in <code>s</code>. If
	 *            the models are trained directly, the order of given models
	 *            during initiation of this assessment object determines, which
	 *            sample will be used for training which model. In general the
	 *            first model will be trained using the first sample in
	 *            <code>s</code> ... . <br>
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
	 * @param mp
	 *            defines which performance measure should be used to assess
	 *            classifiers
	 * @param assessPS
	 *            contains some parameters necessary for assessment (depends on
	 *            the kind of assessment!)
	 * 
	 * @return a {@link ListResult} that contains the results (mean and standard
	 *         errors) of user specified performance measures. These performance
	 *         measures are user specified via the given
	 *         {@link MeasureParameters}.
	 * 
	 * @throws IllegalArgumentException
	 *             if the given <code>assessPS</code> is not of the right type
	 *             (see method <code>evaluateClassifier( ... )</code>)
	 * @throws WrongAlphabetException
	 *             if the given samples <code>s</code> do not use the same
	 *             {@link AlphabetContainer} as contained classifiers/models
	 * @throws Exception
	 *             forwarded from training/testing of classifiers/models
	 * 
	 * @see ClassifierAssessment#assess(MeasureParameters,
	 *      ClassifierAssessmentAssessParameterSet, ProgressUpdater, Sample...)
	 */
	public ListResult assess( MeasureParameters mp, ClassifierAssessmentAssessParameterSet assessPS, Sample... s ) throws IllegalArgumentException,
			WrongAlphabetException,
			Exception {
		return assess( mp, assessPS, null, s );
	}

	/**
	 * Assesses the contained classifiers.
	 * 
	 * @param s
	 *            contains the data to be used for assessment. The order of the
	 *            samples is important. <br>
	 *            If model based classifiers are trained, the order of the
	 *            models in the classifiers determines which model will be
	 *            trained using which sample. The first model in the classifier
	 *            will be trained using the first sample in <code>s</code>. If
	 *            the models are trained directly, the order of given models
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
	 * @param mp
	 *            defines which performance measure should be used to assess
	 *            classifiers
	 * @param assessPS
	 *            contains some parameters necessary for assessment (depends on
	 *            the kind of assessment!)
	 * @param pU
	 *            this {@link ProgressUpdater} may be used to cancel this method
	 *            <code>assess()</code> by setting
	 *            <code>pU.isCancelled()=true</code>. In that case,
	 *            <code>assess()</code> will abort but return results already
	 *            computed.<br>
	 *            In certain cases aborting a classifier assessment will not be
	 *            allowed for example in case of {@link KFoldCrossValidation}.
	 *            In this case it might be wise to override this method such
	 *            that it just returns an error message. <br>
	 *            <code>pU</code> is allowed to be <code>null</code> although in
	 *            this case it may be more convenient to use the second method
	 *            <code>code()</code> not requiring a {@link ProgressUpdater} .
	 * @return a {@link ListResult} that contains the results (mean and standard
	 *         errors) of user specified performance measures. These performance
	 *         measures are user specified via the given
	 *         {@link MeasureParameters}.
	 * @throws IllegalArgumentException
	 *             if the given <code>assessPS</code> is not of the right type
	 *             (see method <code>evaluateClassifier( ... )</code>)
	 * @throws WrongAlphabetException
	 *             if the given samples <code>s</code> do not use the same
	 *             {@link AlphabetContainer} as contained classifiers/models
	 * @throws Exception
	 *             forwarded from training/testing of classifiers/models
	 */
	public ListResult assess( MeasureParameters mp, ClassifierAssessmentAssessParameterSet assessPS, ProgressUpdater pU, Sample... s ) throws IllegalArgumentException,
			WrongAlphabetException,
			Exception {
		if( !mp.isNumeric() ) {
			throw new IllegalArgumentException();//TODO
		}
		if( pU == null ) {
			pU = NullProgressUpdater.getImmutableInstance();
		}

		prepareAssessment( s );

		LinkedList<Result> annotation = new LinkedList<Result>();
		annotation.add( new CategoricalResult( "kind of assessment", "a description or name of the assessment", getNameOfAssessment() ) );
		annotation.addAll( assessPS.getAnnotation() );
		annotation.add( new CategoricalResult( "samples", "annotation of used samples", Sample.getAnnotation( s ) ) );

		this.evaluateClassifier( mp, assessPS, s, pU );

		return new ListResult( this.getNameOfAssessment(),
				"the results of a " + this.getNameOfAssessment(),
				new ResultSet( annotation ),
				this.myTempMeanResultSets );
	}

	/**
	 * Assesses the contained classifiers. In contrast to the other
	 * <code>assess()</code>-methods this one allows the user to predefine all
	 * train and test data beforehand which should be used for assessment.
	 * 
	 * @param s
	 *            contains the data to be used for assessment.<br>
	 *            The order of the samples in <code>s</code> are important.<br>
	 *            <code>s[iter][train/test][]</code> -> the first dimension
	 *            codes for which samples (train, test) are used in iteration
	 *            <code>iter</code>. <br>
	 *            The second dimension codes for training:
	 *            <code>s[iter][0]</code> or test: <code>s[iter][1]</code>.
	 *            <code>s[iter][0]</code> contains for each class a training
	 *            sample. Analog <code>s[iter][1]</code> contains the test
	 *            samples. The order of the samples is important. For further
	 *            details see comment of method
	 *            {@link #assess(MeasureParameters, ClassifierAssessmentAssessParameterSet, Sample...)}
	 *            .<br>
	 *            The user is responsible to take care or not to take care of
	 *            the given test and training dataset to be not overlapping.
	 * @param mp
	 *            defines which performance measure should be used to assess
	 *            classifiers
	 * @param assessPS
	 *            contains some parameters necessary for assessment. Must be of
	 *            type {@link ClassifierAssessmentAssessParameterSet}
	 * @param pU
	 *            this {@link ProgressUpdater} allows to abort this classifier
	 *            assessment. If <code>pU.isCancelled()=true</code>, all results
	 *            already computed will be returned. It is allowed to give a
	 *            null reference.
	 * 
	 * @return a {@link ListResult} that contains the results (mean and standard
	 *         errors) of user specified performance measures. These performance
	 *         measures are user specified via the given
	 *         {@link MeasureParameters}.
	 * @throws IllegalArgumentException
	 *             if the given <code>assessPS</code> is not of the right type
	 *             (see method <code>evaluateClassifier( ... )</code>)
	 * @throws WrongAlphabetException
	 *             if the given samples <code>s</code> do not use the same
	 *             {@link AlphabetContainer} as contained classifiers/models
	 * @throws Exception
	 *             forwarded from training/testing of classifiers/models
	 * 
	 */
	public ListResult assess( MeasureParameters mp, ClassifierAssessmentAssessParameterSet assessPS, ProgressUpdater pU, Sample[][]... s ) throws IllegalArgumentException,
			WrongAlphabetException,
			Exception {
		if( !mp.isNumeric() ) {
			throw new IllegalArgumentException();//TODO
		}
		if( pU == null ) {
			pU = NullProgressUpdater.getImmutableInstance();
		}

		if( s == null ) {
			throw new IllegalArgumentException( "Please check the input: No samples are given!" );
		}

		int i, j;
		int numOfClasses = this.myAbstractClassifier[0].getNumberOfClasses();
		int eL = assessPS.getElementLength();
		boolean exceptionIfMPNotComputable = assessPS.getExceptionIfMPNotComputable();

		Sample[][][] correctedS = new Sample[s.length][2][s[0][0].length];

		AlphabetContainer abc = this.myAbstractClassifier[0].getAlphabetContainer();
		for( i = 0; i < s.length; i++ ) {
			if( s[i].length != 2 ) {
				throw new IllegalArgumentException( "For each iteration the given sample-3d-array has to contain " + "exactly one array of training-samples and one array of test-samples but it doesn't." );
			}

			if( s[i][0].length != numOfClasses ) {
				throw new IllegalArgumentException( "For each iteration the given sample-3d-array has to contain " + "a array of training-samples containing one training-sample for each class but it doesn't." );
			}

			if( s[i][1].length != numOfClasses ) {
				throw new IllegalArgumentException( "For each iteration the given sample-3d-array has to contain " + "a array of test-samples containing one test-sample for each class but it doesn't." );
			}

			for( j = 0; j < s[i][0].length; j++ ) {
				for( int k=0; k < 2; k++ ) {
					abc.checkConsistency( s[i][k][j].getAlphabetContainer() );
					if( s[i][k][j].getElementLength() != eL ) {
						correctedS[i][k][j] = new Sample( s[i][k][j], eL );
					} else {
						correctedS[i][k][j] = s[i][k][j];
					}
				}
			}
		}

		this.myTempMeanResultSets = new MeanResultSet[this.myAbstractClassifier.length];
		for( i = 0; i < this.myAbstractClassifier.length; i++ ) {
			this.myTempMeanResultSets[i] = new MeanResultSet( this.myAbstractClassifier[i].getClassifierAnnotation() );
		}

		// evaluate

		LinkedList<Result> annotation = new LinkedList<Result>();
		annotation.add( new CategoricalResult( "kind of assessment",
				"a description or name of the assessment",
				"assessment using a series of user-given pairs of test- and train-datasets" ) );

		pU.setMax( s.length );

		for( i = 0; i < s.length; i++ ) {
			this.train( correctedS[i][0] );
			this.test( mp, exceptionIfMPNotComputable, correctedS[i][1] );

			pU.setValue( i + 1 );
			if( pU.isCancelled() ) {
				break;
			}
		}

		return new ListResult( this.getNameOfAssessment(),
				"the results of an assessment using a series of user-given pairs of test- and train-datasets",
				new ResultSet( annotation ),
				this.myTempMeanResultSets );
	}

	/**
	 * Returns a deep copy of all classifiers that have been or will be used in
	 * this assessment.
	 * 
	 * @return a deep copy of all used classifiers in this assessment
	 * 
	 * @throws CloneNotSupportedException
	 *             if it is impossible to get a deep copy for at least one
	 *             classifier (if the classifier could not be cloned)
	 */
	public AbstractClassifier[] getClassifier() throws CloneNotSupportedException {
		AbstractClassifier[] res = new AbstractClassifier[myAbstractClassifier.length];
		for( int i = 0; i < res.length; i++ ) {
			res[i] = myAbstractClassifier[i].clone();
		}
		return res;
	}

	/**
	 * Returns the name of this class.
	 * 
	 * @return name of this class
	 */
	public String getNameOfAssessment() {
		return this.getClass().getSimpleName();
	}

	/**
	 * This method must be implemented in all subclasses. It should perform the
	 * following tasks: <br>
	 * 1.) create test and train datasets <br>
	 * 2.) call method <code>train()</code> to train classifiers/models using
	 * train data <br>
	 * 3.) call method <code>test()</code> to cause evaluation (test) of trained
	 * classifiers
	 * 
	 * @param mp
	 *            defines which performance measures are used to assess
	 *            classifiers
	 * @param assessPS
	 *            contains assessment specific parameters (like: number of
	 *            iterations of a k-fold-crossvalidation)
	 * @param s
	 *            data to be used for assessment (both: test and train data)
	 * @param pU
	 *            a {@link ProgressUpdater} that mainly has to be used to allow
	 *            the user to cancel a current running classifier assessment.
	 *            This {@link ProgressUpdater} is guaranteed to be not
	 *            <code>null</code>. In certain cases aborting a classifier
	 *            assessment will not be allowed for example in case of
	 *            {@link KFoldCrossValidation}. In this case the given
	 *            {@link ProgressUpdater} should be ignored. <br>
	 * <br>
	 *            Usage:<br>
	 *            <ul>
	 *            <li>
	 *            <code>pU.setMax()</code>= number of iterations of the assessment loop
	 *            <li>iteration=0;
	 *            <li>assessment loop
	 *            <ul>
	 *            <li> <code>pU.setValue()</code>=iteration+1;
	 *            <li> sample treatment
	 *            <li> <code>train()</code>;
	 *            <li> <code>test()</code>;
	 *            <li>iteration++;
	 *            </ul>
	 *            <li>repeat unless(ready or not(<code>pU.isCancelled()</code>))
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if the given {@link ClassifierAssessmentAssessParameterSet}
	 *             is of wrong type
	 * @throws Exception
	 *             that occurred during training or using classifiers/models
	 */
	protected abstract void evaluateClassifier( MeasureParameters mp, ClassifierAssessmentAssessParameterSet assessPS, Sample[] s,
			ProgressUpdater pU ) throws IllegalArgumentException, Exception;

	/**
	 * Prepares an assessment. If the given {@link Sample} may not be used for
	 * this assessment, this method throws an {@link Exception}. <br>
	 * Further {@link MeanResultSet}s are initiated for this assessment (one for
	 * each contained classifier).
	 * 
	 * @param s
	 *            the {@link Sample} to be checked
	 *            
	 * @throws WrongAlphabetException
	 *             if <br>
	 *             <ol>
	 *             <li> <code>s</code> is <code>null</code> or not of required
	 *             length (number of classes)</li>
	 *             <li> {@link AlphabetContainer}s of <code>s</code> are not
	 *             consistent with {@link AlphabetContainer} of local models or
	 *             classifiers</li>
	 *             </ol>
	 * @throws IllegalArgumentException
	 * 				if the given samples are not suitable
	 */
	protected void prepareAssessment( Sample... s ) throws IllegalArgumentException, WrongAlphabetException {

		if( s == null || s.length != this.myAbstractClassifier[0].getNumberOfClasses() ) {
			throw new IllegalArgumentException( "Either no samples are given or the number of samples " + "is not equal to the number of different classes "
												+ "the local classifiers are able to distinguish" );
		}

		AlphabetContainer abc = this.myAbstractClassifier[0].getAlphabetContainer();
		for( int i = 0; i < s.length; i++ ) {
			if( !abc.checkConsistency( s[i].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException( "The AlphabetContainer of Sample " + i + " does not match with tha AlphabetContainer of the first classifier." );
			}
		}

		this.myTempMeanResultSets = new MeanResultSet[this.myAbstractClassifier.length];
		for( int i = 0; i < this.myAbstractClassifier.length; this.myTempMeanResultSets[i] = new MeanResultSet( this.myAbstractClassifier[i++].getClassifierAnnotation() ) );
	}

	/**
	 * Uses the given test samples to call the <code>evaluate( ... )</code>
	 * -methods of the local {@link AbstractClassifier}s. The returned
	 * {@link de.jstacs.results.NumericalResult}s as well as the numerical
	 * characteristics are added to each classifiers {@link MeanResultSet}. <br>
	 * It should not be necessary to override this method in subclasses.
	 * 
	 * @param mp
	 *            determines which performance measures are used to assess the
	 *            classifiers
	 * @param exception
	 *            whether an {@link Exception} should be thrown if some
	 *            {@link de.jstacs.classifier.MeasureParameters.Measure} could
	 *            not be evaluated
	 * @param testS
	 *            samples used as test sets (has to contain one {@link Sample}
	 *            for each class)
	 * 
	 * @throws IllegalValueException
	 *             if a parameter is not valid
	 * @throws InconsistentResultNumberException
	 *             if the number of results between the different result sets
	 *             differ
	 * @throws AdditionImpossibleException
	 *             if added result sets do not match
	 * @throws Exception
	 *             if necessary
	 * @throws IllegalArgumentException
	 *             if the length of <code>testS</code> is not equal to the
	 *             dimension of the classification problem (
	 *             <code>testS.length!=this.myAbstractClassifier
	 *             [0].getNumberOfClasses()</code>)
	 * 
	 * @see AbstractClassifier#evaluate(MeasureParameters, boolean, Sample...)
	 */
	protected void test( MeasureParameters mp, boolean exception, Sample... testS ) throws IllegalValueException,
			InconsistentResultNumberException,
			AdditionImpossibleException,
			Exception {

		if( testS.length != this.myAbstractClassifier[0].getNumberOfClasses() ) {
			throw new IllegalArgumentException( "Dimension of given testSample-array is not " + "equal to problem-dimension (classifier.getNumberOfClasses)." );
		}

		for( int i = 0; i < this.myAbstractClassifier.length; i++ ) {
			this.myTempMeanResultSets[i].addResults( (NumericalResultSet) this.myAbstractClassifier[i].evaluate( mp, exception, testS ),
					this.myAbstractClassifier[i].getNumericalCharacteristics() );

		}
	}

	/**
	 * Trains the local classifiers using the given training samples. <br>
	 * The classifiers are either directly trained or via training of the local
	 * models. The second option always is used if the
	 * {@link ClassifierAssessment}-object was constructed using
	 * {@link Model}s. <br>
	 * <br>
	 * It should not be necessary to override this method in subclasses.
	 * 
	 * @param trainS
	 *            samples used as training sets (has to contain one
	 *            {@link Sample} for each class)
	 * 
	 * @throws IllegalArgumentException
	 *             if the length of <code>trainS</code> is not equal to the
	 *             dimension of the classification problem (
	 *             <code>trainS.length!=this.myAbstractClassifier
	 *             [0].getNumberOfClasses()</code>)
	 * @throws Exception
	 *             if necessary
	 */
	protected void train( Sample... trainS ) throws IllegalArgumentException, Exception {

		if( trainS.length != this.myAbstractClassifier[0].getNumberOfClasses() ) {
			throw new IllegalArgumentException( "Dimension of given trainSample-array is not " + "equal to problem-dimension (classifier.getNumberOfClasses)." );
		}

		for( int i = 0; i < this.myAbstractClassifier.length - this.skipLastClassifiersDuringClassifierTraining; i++ ) {
			this.myAbstractClassifier[i].train( trainS );
		}

		//FIXME 
		//hier werden Klassifikatoren gelernt indem die enthaltenen Modelle gelernt werden
		//Problem: die entsprechenden Klassifikatoren haben dann noch keine KlassenWahrscheinlichkeiten
		//gelernt denn dies wuerden diese nur machen, wenn sie direkt gelernt werden

		// train classifier via training of models
		for( int classes = 0; classes < this.myModel.length; classes++ ) {
			for( int models = 0; models < this.myModel[classes].length; this.myModel[classes][models++].train( trainS[classes] ) );
		}
	}

	/**
	 * Used only in {@link ClassifierAssessment} to allow the construction of
	 * {@link ModelBasedClassifier} without cloning of the {@link Model}
	 * s.
	 * 
	 * @author Jens Keilwagen
	 */
	private static class ModelBasedAssessmentClassifier extends ModelBasedClassifier {

		private ModelBasedAssessmentClassifier( Model... models ) throws IllegalArgumentException, CloneNotSupportedException,
																	ClassDimensionException {
			super( false, models );
		}

		private ModelBasedAssessmentClassifier( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}

		/**
		 * Returns an instance of {@link ModelBasedClassifier}, since cloning
		 * for this class does not make sense and complicates saving and loading
		 * (the clone).
		 * 
		 * @return an instance of {@link ModelBasedClassifier}
		 * 
		 * @throws CloneNotSupportedException
		 *             if something went wrong while cloning
		 */
		@Override
		public ModelBasedClassifier clone() throws CloneNotSupportedException {
			try {
				ModelBasedClassifier mbc = new ModelBasedClassifier( models );
				mbc.setClassWeights( false, this.getClassWeights() );
				return mbc;
			} catch ( ClassDimensionException cde ) {
				// does not happen
				throw new CloneNotSupportedException( cde.getMessage() );
			}
		}
	}
}