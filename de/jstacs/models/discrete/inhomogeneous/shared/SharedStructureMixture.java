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

package de.jstacs.models.discrete.inhomogeneous.shared;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.io.XMLParser;
import de.jstacs.models.discrete.DiscreteGraphicalModel;
import de.jstacs.models.discrete.inhomogeneous.FSDAGModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.models.discrete.inhomogeneous.parameters.IDGMParameterSet;
import de.jstacs.models.mixture.MixtureModel;

/**
 * This class handles a mixture of models with the same structure that is
 * learned via EM. One well known example is a <a
 * href="http://people.csail.mit.edu/mmp/thesis.html">mixture of trees</a>.
 * 
 * @author Jens Keilwagen
 */
public class SharedStructureMixture extends MixtureModel {

	private static final long serialVersionUID = -9019277262448497915L;

	private StructureLearner sl;

	private ModelType modelType;

	private byte order;

	private LearningType method;

	private static double[] getClassHyperParams( DiscreteGraphicalModel[] m ) {
		double[] erg = new double[m.length];
		for( int i = 0; i < erg.length; i++ ) {
			erg[i] = m[i].getESS();
		}
		return erg;
	}

	/**
	 * Creates a new {@link SharedStructureMixture} instance which estimates the
	 * component probabilities/weights.
	 * 
	 * @param m
	 *            the single models building the mixture model
	 * @param model
	 *            the type of the model
	 * @param order
	 *            the order of the model
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param alpha
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas, it is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex)
	 * @param tc
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequences of same
	 *             length
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null && weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> (hyperparameters for
	 *             the component assignment prior) are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all models work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the models can not be cloned
	 * 
	 * @see ModelType
	 * @see SharedStructureMixture#SharedStructureMixture(FSDAGModel[],
	 *      de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType, byte, int, boolean, double[], double, TerminationCondition)
	 */
	public SharedStructureMixture( FSDAGModel[] m, ModelType model, byte order, int starts, double alpha, TerminationCondition tc )
																															throws IllegalArgumentException,
																															WrongAlphabetException,
																															CloneNotSupportedException {
		this( m, model, order, starts, true, null, alpha, tc );
	}

	/**
	 * Creates a new {@link SharedStructureMixture} instance with fixed
	 * component weights.
	 * 
	 * @param m
	 *            the single models building the mixture model
	 * @param model
	 *            the type of the model
	 * @param order
	 *            the order of the model
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == models.length</code>)
	 * @param alpha
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas, it is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex)
	 * @param tc
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequences of same
	 *             length
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null && weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> (hyperparameters for
	 *             the component assignment prior) are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all models work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the models can not be cloned
	 * 
	 * @see ModelType
	 * @see SharedStructureMixture#SharedStructureMixture(FSDAGModel[],
	 *      de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType, byte, int, boolean, double[], double, TerminationCondition)
	 */
	public SharedStructureMixture( FSDAGModel[] m, ModelType model, byte order, int starts, double[] weights, double alpha, TerminationCondition tc )
																																			throws IllegalArgumentException,
																																			WrongAlphabetException,
																																			CloneNotSupportedException {
		this( m, model, order, starts, false, weights, alpha, tc );
	}

	/**
	 * Creates a new {@link SharedStructureMixture} instance with all relevant
	 * values. This constructor is used by the other main constructors.
	 * 
	 * @param m
	 *            the single models building the mixture model
	 * @param model
	 *            the type of the model
	 * @param order
	 *            the order of the model
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param estimateComponentProbs
	 *            the switch for estimating the component probabilities in the
	 *            algorithm or to hold them fixed; if the component parameters
	 *            are fixed, the values of <code>weights</code> will be used,
	 *            otherwise the <code>componentHyperParams</code>
	 *            (hyperparameters for the component assignment prior) will be
	 *            incorporated in the adjustment
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == models.length</code>)
	 * @param alpha
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas, it is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex)
	 * @param tc
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequences of same
	 *             length
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null && weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> (hyperparameters for
	 *             the component assignment prior) are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all models work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the models can not be cloned
	 * 
	 * @see ModelType
	 * @see MixtureModel#MixtureModel(int, de.jstacs.models.Model[], int,
	 *      boolean, double[], double[], de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double, TerminationCondition,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int, int,
	 *      de.jstacs.sampling.BurnInTest)
	 */
	protected SharedStructureMixture( FSDAGModel[] m, ModelType model, byte order, int starts, boolean estimateComponentProbs,
										double[] weights, double alpha, TerminationCondition tc ) throws IllegalArgumentException,
																						WrongAlphabetException, CloneNotSupportedException {
		super( m[0].getLength(),
				m,
				starts,
				estimateComponentProbs,
				getClassHyperParams( m ),
				weights,
				Algorithm.EM,
				alpha,
				tc,
				Parameterization.LAMBDA,
				0,
				0,
				null );
		sl = new StructureLearner( m[0].getAlphabetContainer(), getLength() );
		this.modelType = model;
		if( order < 0 ) {
			throw new IllegalArgumentException( "The value of order has to be non-negative." );
		}
		this.order = order;
		this.method = LearningType.ML_OR_MAP;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SharedStructureMixture} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SharedStructureMixture} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see MixtureModel#MixtureModel(StringBuffer)
	 */
	public SharedStructureMixture( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#clone()
	 */
	@Override
	public SharedStructureMixture clone() throws CloneNotSupportedException {
		SharedStructureMixture clone = (SharedStructureMixture)super.clone();
		clone.sl = new StructureLearner( alphabets, length );
		return clone;
	}

	/**
	 * Returns a {@link String} representation of the structure of the used
	 * models.
	 * 
	 * @return a {@link String} representation of the structure of the used
	 *         models
	 * 
	 * @throws NotTrainedException
	 *             if the classifier is not trained yet
	 * 
	 * @see FSDAGModel#getStructure()
	 */
	public String getStructure() throws NotTrainedException {
		return ( (FSDAGModel)model[0] ).getStructure();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "ssMixModel(" + dimension
				+ " "
				+ IDGMParameterSet.getModelInstanceName( modelType, order, method, ( (FSDAGModel)model[0] ).getESS() )
				+ ")";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( xml, modelType, "model" );
		XMLParser.appendObjectWithTags( xml, order, "order" );
		XMLParser.appendObjectWithTags( xml, method, "method" );
		xml.append( super.toXML() );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	private static final String XML_TAG = "SharedStructureMixture";

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( representation, XML_TAG );
		modelType = XMLParser.extractObjectForTags( xml, "model", ModelType.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		order = XMLParser.extractObjectForTags( xml, "order", byte.class );
		method = XMLParser.extractObjectForTags( xml, "method", LearningType.class );
		sl = new StructureLearner( getAlphabetContainer(), getLength() );
		super.fromXML( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#getNewParameters(int, double[][], double[])
	 */
	@Override
	protected void getNewParameters( int iteration, double[][] seqWeights, double[] w ) throws Exception {
		SymmetricTensor[] parts = new SymmetricTensor[dimension];
		double[] x = new double[dimension];
		for( int i = 0; i < dimension; i++ ) {
			sl.setESS( ( (FSDAGModel)model[i] ).getESS() );
			parts[i] = sl.getTensor( sample[0], seqWeights[i], order, method );
			x[i] = 1d;
		}
		FSDAGModel.train( model, StructureLearner.getStructure( new SymmetricTensor( parts, x ), modelType, order ), seqWeights, sample[0] );

		getNewComponentProbs( w );
	}
}
