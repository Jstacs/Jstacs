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

package de.jstacs.scoringFunctions.mix;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction;
import de.jstacs.scoringFunctions.NormalizableScoringFunction;
import de.jstacs.scoringFunctions.NormalizedScoringFunction;
import de.jstacs.scoringFunctions.SamplingScoringFunction;
import de.jstacs.scoringFunctions.mix.motifSearch.DurationScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This main abstract class for any mixture scoring function (e.g.
 * &quot;real&quot; mixture, strand mixture, hidden motif mixture, ...). The
 * potential for the hidden variables is parameterized depending on the
 * parameterization of the given {@link NormalizableScoringFunction}s. If these
 * are already normalized (see
 * {@link NormalizableScoringFunction#isNormalized()} ) the potential is
 * parameterized using the Meila-parameterization, otherwise it is parameterized
 * using the unnormalized MRF(Markov Random Fields)-parameterization.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractMixtureScoringFunction extends AbstractNormalizableScoringFunction implements SamplingScoringFunction {

	private int starts;

	/**
	 * This array contains the references/indices for the parameters. Only the
	 * start index for each new function is stored.
	 */
	protected int[] paramRef;

	/**
	 * This <code>boolean</code> indicates whether to optimize the hidden
	 * variables of this instance. (It is not used recursive.)
	 */
	protected boolean optimizeHidden;

	/**
	 * This <code>boolean</code> indicates whether free parameterization or all
	 * parameters are used.
	 */
	protected boolean freeParams;

	/**
	 * This <code>boolean</code> indicates whether to use a plug-in strategy to
	 * initialize the instance.
	 */
	private boolean plugIn;

	/**
	 * This array contains the internal
	 * {@link de.jstacs.scoringFunctions.ScoringFunction}s that are used to
	 * determine the score.
	 */
	protected NormalizableScoringFunction[] function;

	/**
	 * This array contains the hidden parameters of the instance.
	 */
	protected double[] hiddenParameter;

	/**
	 * This array contains the logarithm of the hidden potentials of the
	 * instance.
	 */
	protected double[] logHiddenPotential;

	/**
	 * This array contains the hidden potentials of the instance.
	 */
	protected double[] hiddenPotential;

	/**
	 * This array is used while computing the score. It stores the scores of the
	 * components and is used to avoid creating a new array every time.
	 */
	protected double[] componentScore;

	/**
	 * This array contains the partial normalization constants, i.e. the
	 * normalization constant for each component.
	 */
	protected double[] partNorm;

	/**
	 * This <code>double</code> contains the normalization constant of the
	 * instance.
	 */
	protected double norm;

	/**
	 * This <code>double</code> contains the logarithm of the normalization
	 * constant of hidden parameters of the instance.
	 */
	protected double logHiddenNorm;

	/**
	 * This <code>double</code> contains the sum of the logarithm of the gamma
	 * functions used in the prior.
	 * 
	 * @see AbstractMixtureScoringFunction#computeLogGammaSum()
	 */
	protected double logGammaSum;

	/**
	 * This array contains some {@link DoubleList}s that are used while
	 * computing the partial derivation.
	 */
	protected DoubleList[] dList;

	/**
	 * This array contains some {@link IntList}s that are used while computing
	 * the partial derivation.
	 */
	protected IntList[] iList;

	/**
	 * This boolean indicates whether this instance is a normalized one or not.
	 */
	private boolean isNormalized;

	/**
	 * This constructor creates a new {@link AbstractMixtureScoringFunction}.
	 * 
	 * @param length
	 *            the sequence length that should be modeled
	 * @param starts
	 *            the number of starts that should be done in an optimization
	 * @param dimension
	 *            the number of different mixture components
	 * @param optimizeHidden
	 *            indicates whether the parameters for the hidden variables
	 *            should be optimized or not
	 * @param plugIn
	 *            indicates whether the initial parameters for an optimization
	 *            should be related to the data or randomly drawn
	 * @param function
	 *            the {@link de.jstacs.scoringFunctions.ScoringFunction}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element of <code>function</code> could not be cloned
	 */
	protected AbstractMixtureScoringFunction( int length, int starts, int dimension, boolean optimizeHidden, boolean plugIn,
											NormalizableScoringFunction... function ) throws CloneNotSupportedException {
		super( function[0].getAlphabetContainer(), length );
		this.function = ArrayHandler.clone( function );
		if( starts < 1 ) {
			throw new IllegalArgumentException( "The number of recommended starts has to be positive." );
		} else {
			this.starts = starts;
		}
		if( dimension == 0 ) {
			throw new IllegalArgumentException( "The number of components has to be positive." );
		}

		isNormalized = determineIsNormalized();
		hiddenParameter = new double[dimension];
		logHiddenPotential = new double[dimension];
		hiddenPotential = new double[dimension];
		partNorm = new double[dimension];
		setHiddenParameters( hiddenParameter, 0 );

		componentScore = new double[dimension];
		this.optimizeHidden = optimizeHidden && dimension > 1;
		this.plugIn = plugIn;
		paramRef = null;
		init( freeParams );
		norm = Double.NaN;
	
	}

	/**
	 * This method is used to pre-compute the sum of the logarithm of the gamma
	 * functions that is used in the prior.
	 */
	protected void computeLogGammaSum() {
		logGammaSum = 0;
		int i = 0, n = getNumberOfComponents();
		if( n > 1 && getEss() > 0 ) {
			double sum = 0, h;
			for( ; i < n; i++ ) {
				h = getHyperparameterForHiddenParameter( i );
				sum += h;
				logGammaSum -= Gamma.logOfGamma( h );
			}
			logGammaSum += Gamma.logOfGamma( sum );
		}
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AbstractMixtureScoringFunction} out of a
	 * {@link StringBuffer} as returned by {@link #toXML()}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the representation could not be parsed
	 */
	protected AbstractMixtureScoringFunction( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#clone()
	 */
	@Override
	public AbstractMixtureScoringFunction clone() throws CloneNotSupportedException {
		AbstractMixtureScoringFunction clone = (AbstractMixtureScoringFunction)super.clone();
		clone.cloneFunctions( function );
		clone.hiddenParameter = hiddenParameter.clone();
		clone.logHiddenPotential = logHiddenPotential.clone();
		clone.hiddenPotential = hiddenPotential.clone();
		clone.componentScore = componentScore.clone();
		clone.partNorm = partNorm.clone();
		clone.iList = null;
		clone.paramRef = null;
		clone.init( freeParams );
		return clone;
	}

	/**
	 * This method clones the given array of functions and enables the user to
	 * do some post-processing. This method is only used in {@link #clone()}.
	 * 
	 * @param originalFunctions
	 *            the array of functions to be cloned
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element of <code>originalFunctions</code> could not be
	 *             cloned
	 */
	protected void cloneFunctions( NormalizableScoringFunction[] originalFunctions ) throws CloneNotSupportedException {
		function = ArrayHandler.clone( originalFunctions );
	}

	/**
	 * This method returns the hyperparameter for the hidden parameter with
	 * index <code>index</code>.
	 * 
	 * @param index
	 *            the index of the hidden parameter
	 * 
	 * @return the hyperparameter for the hidden parameter
	 */
	public abstract double getHyperparameterForHiddenParameter( int index );

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.NormalizableScoringFunction#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		double val = 0, h, sum = 0;
		for( int i = 0; i < hiddenParameter.length; i++ ) {
			h = getHyperparameterForHiddenParameter( i );
			sum += h;
			val += hiddenParameter[i] * h;
		}
		if( isNormalized() ) {
			val -= sum * logHiddenNorm;
		}
		for( int i = 0; i < function.length; i++ ) {
			if( function[i] != null ) {
				val += function[i].getLogPriorTerm();
			}
		}
		return val + logGammaSum;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		int i = 0, j;
		for( ; i < function.length; i++ ) {
			if( function[i] != null ) {
				function[i].addGradientOfLogPriorTerm( grad, start + paramRef[i] );
			}
		}
		j = start + paramRef[function.length + 1];
		start += paramRef[function.length];
		double e = getEss();
		for( i = 0; start < j; i++, start++ ) {
			grad[start] += getHyperparameterForHiddenParameter( i ) - ( isNormalized() ? e * hiddenPotential[i] : 0 );
		}
	}

	/**
	 * Returns the index of the component that has the greatest impact on the
	 * complete score for a {@link Sequence}.
	 * 
	 * @param seq
	 *            the sequence
	 * @param start
	 *            the start position
	 * 
	 * @return the index of the component
	 */
	public int getIndexOfMaximalComponentFor( Sequence seq, int start ) {
		fillComponentScores( seq, start );
		return ToolBox.getMaxIndex( componentScore );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		int numPars = this.getNumberOfParameters();
		if( numPars == UNKNOWN ) {
			throw new Exception( "No parameters exists, yet." );
		} else {
			double[] part, current = new double[numPars];
			int i = 0, j = function.length;
			while( i < j ) {
				if( function[i] != null ) {
					part = function[i].getCurrentParameterValues();
					System.arraycopy( part, 0, current, paramRef[i], part.length );
				}
				i++;
			}
				
			System.arraycopy( hiddenParameter, 0, current, paramRef[j], paramRef[j + 1] - paramRef[j] );
			return current;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data.Sequence, int)
	 */
	public double getLogScore( Sequence seq, int start ) {
		fillComponentScores( seq, start );
		return Normalisation.getLogSum( componentScore );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getNormalizationConstant()
	 */
	public final double getLogNormalizationConstant() {
		if( isNormalized() ) {
			return 0;
		} else {
			if( Double.isNaN( norm ) ) {
				precomputeNorm();
			}
			return norm;
		}
	}

	/**
	 * Returns the number of different components of this
	 * {@link AbstractMixtureScoringFunction}.
	 * 
	 * @return the number of different components
	 */
	public final int getNumberOfComponents() {
		return componentScore.length;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	public final int getNumberOfParameters() {
		if( paramRef == null ) {
			return UNKNOWN;
		} else {
			return paramRef[paramRef.length - 1];
		}
	}

	public final int getNumberOfRecommendedStarts() {
		return starts;
	}

	/**
	 * Returns the probabilities for each component given a {@link Sequence}.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * 
	 * @return an array containing the probability of component <code>i</code>
	 *         (<code>=p(i|class,seq)</code>) in entry <code>i</code>
	 */
	public double[] getProbsForComponent( Sequence seq ) {
		fillComponentScores( seq, 0 );
		double[] p = new double[componentScore.length];
		Normalisation.logSumNormalisation( componentScore, 0, p.length, p, 0 );
		return p;
	}

	/**
	 * Returns a deep copy of all internal used
	 * {@link de.jstacs.scoringFunctions.ScoringFunction}s.
	 * 
	 * @return a deep copy of all internal used
	 *         {@link de.jstacs.scoringFunctions.ScoringFunction}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if one element of the internal used
	 *             {@link de.jstacs.scoringFunctions.ScoringFunction}s could not
	 *             be cloned
	 */
	public NormalizableScoringFunction[] getScoringFunctions() throws CloneNotSupportedException {
		return ArrayHandler.clone( function );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		int[] ind = getIndices( index );
		if( ind[0] == function.length ) {
			return hiddenParameter.length;
		} else {
			return function[ind[0]].getSizeOfEventSpaceForRandomVariablesOfParameter( ind[1] );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int, boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction( int index, boolean freeParams, Sample[] data, double[][] weights ) throws Exception {
		if( plugIn ) {
			initializeUsingPlugIn( index, freeParams, data, weights );
			init( freeParams );
		} else {
			initializeFunctionRandomly( freeParams );
		}
	}

	/**
	 * This method initializes the functions using the data in some way.
	 * 
	 * @param index
	 *            the class index
	 * @param freeParams
	 *            if <code>true<code>, the (reduced) parameterization is used
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights for the data
	 * 
	 * @throws Exception
	 *             if the initialization could not be done
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 *      boolean, Sample[], double[][])
	 */
	protected abstract void initializeUsingPlugIn( int index, boolean freeParams, Sample[] data, double[][] weights ) throws Exception;

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly(boolean)
	 */
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		for( int i = 0; i < function.length; i++ ) {
			if( function[i] != null ) {
				function[i].initializeFunctionRandomly( freeParams );
			}
		}
		if( optimizeHidden ) {
			initializeHiddenPotentialRandomly();
		}
		init( freeParams );
	}

	/**
	 * This method initializes the hidden potential (and the corresponding
	 * parameters) randomly.
	 */
	protected void initializeHiddenPotentialRandomly() {
		double[] h = new double[this.getNumberOfComponents()];
		if( getEss() == 0 ) {
			Arrays.fill( h, 1 );
		} else {
			for( int j = 0; j < h.length; j++ ) {
				h[j] = getHyperparameterForHiddenParameter( j );
			}
		}
		DirichletMRGParams param = new DirichletMRGParams( h );
		DirichletMRG.DEFAULT_INSTANCE.generate( hiddenPotential, 0, hiddenPotential.length, param );

		computeHiddenParameter( hiddenPotential, false );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	public boolean isInitialized() {
		int i = 0;
		while( i < function.length && (function[i] == null || function[i].isInitialized()) ) {
			i++;
		}
		return paramRef != null && i == function.length;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[], int)
	 */
	public void setParameters( double[] params, int start ) {
		int i = 0;
		for( ; i < function.length; i++ ) {
			if( function[i] != null ) {
				setParametersForFunction( i, params, start + paramRef[i] );
			}
		}
		isNormalized = determineIsNormalized();
		if( optimizeHidden ) {
			setHiddenParameters( params, start + paramRef[i] );
		} else {
			if( isNormalized ) {
				norm = 0;
			} else {
				norm = Double.NaN;
			}
		}
	}
	
	/**
	 * This method is used to determine the value that is returned by the method {@link #isNormalized()}.
	 * 
	 * @return the value that is set to an internal variable and that will be returned by {@link #isNormalized()}
	 */
	protected boolean determineIsNormalized(){
		return isNormalized( function );
	}

	/**
	 * This method initializes the hidden parameters of the instance uniformly.
	 */
	public void initializeHiddenUniformly() {
		int i = 0, c = getNumberOfComponents();
		for( ; i < function.length; i++ ) {
			if( function[i] != null ) {
				if( function[i] instanceof AbstractMixtureScoringFunction ) {
					( (AbstractMixtureScoringFunction)function[i] ).initializeHiddenUniformly();
				} else if( function[i] instanceof NormalizedScoringFunction ) {
					( (NormalizedScoringFunction)function[i] ).initializeHiddenUniformly();
				} else if( function[i] instanceof DurationScoringFunction ) {
					( (DurationScoringFunction)function[i] ).initializeUniformly();
				}
			}
		}
		if( optimizeHidden ) {
			double[] pars = new double[c];
			double d;
			if( freeParams ) {
				d = getLogNormalizationConstantForComponent( c );
			} else {
				d = 0;
			}
			for( i = 0; i < c; i++ ) {
				pars[i] = d - getLogNormalizationConstantForComponent( i );
			}
			setHiddenParameters( pars, 0 );
		}
		init( freeParams );
	}

	/**
	 * This method sets the hidden parameters of the model.
	 * 
	 * @param params
	 *            the parameter vector
	 * @param start
	 *            the start index in <code>params</code>
	 */
	protected void setHiddenParameters( double[] params, int start ) {
		int i, len = hiddenParameter.length - ( freeParams ? 1 : 0 );

		for( i = 0; i < len; i++, start++ ) {
			hiddenParameter[i] = params[start];
		}
		if( freeParams ) {
			hiddenParameter[i] = 0;
		}
		if( isNormalized() ) {
			logHiddenNorm = Normalisation.getLogSum( hiddenParameter );
			norm = 0;
		} else {
			logHiddenNorm = 0;
			norm = Double.NaN;
		}

		for( i = 0; i < logHiddenPotential.length; i++ ) {
			logHiddenPotential[i] = hiddenParameter[i] - logHiddenNorm;
			hiddenPotential[i] = Math.exp( logHiddenPotential[i] );
			partNorm[i] = logHiddenPotential[i];
		}		
	}

	/**
	 * This method allows to set the parameters for specific functions.
	 * 
	 * @param index
	 *            the function index
	 * @param params
	 *            the parameter vector
	 * @param start
	 *            the start index in <code>params</code>
	 */
	protected void setParametersForFunction( int index, double[] params, int start ) {
		function[index].setParameters( params, start );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public final StringBuffer toXML() {
		StringBuffer b = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( b, length, "length" );
		XMLParser.appendObjectWithTags( b, starts, "starts" );
		XMLParser.appendObjectWithTags( b, freeParams, "freeParams" );
		XMLParser.appendObjectWithTags( b, function, "function" );
		XMLParser.appendObjectWithTags( b, optimizeHidden, "optimizeHidden" );
		XMLParser.appendObjectWithTags( b, plugIn, "plugIn" );
		XMLParser.appendObjectWithTags( b, hiddenParameter, "hiddenParameter" );
		b.append( getFurtherInformation() );
		XMLParser.addTags( b, getXMLTag() );
		return b;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#fromXML(java.lang.StringBuffer)
	 */
	protected final void fromXML( StringBuffer b ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( b, getXMLTag() );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		starts = XMLParser.extractObjectForTags( xml, "starts", int.class );
		freeParams = XMLParser.extractObjectForTags( xml, "freeParams", boolean.class );
		function = XMLParser.extractObjectForTags( xml, "function", NormalizableScoringFunction[].class );
		alphabets = function[0].getAlphabetContainer();
		optimizeHidden = XMLParser.extractObjectForTags( xml, "optimizeHidden", boolean.class );
		plugIn = XMLParser.extractObjectForTags( xml, "plugIn", boolean.class );
		hiddenParameter = XMLParser.extractObjectForTags( xml, "hiddenParameter", double[].class );
		hiddenPotential = new double[hiddenParameter.length];
		logHiddenPotential = new double[hiddenParameter.length];
		partNorm = new double[hiddenParameter.length];
		componentScore = new double[hiddenParameter.length];
		extractFurtherInformation( xml );
		isNormalized = determineIsNormalized();
		setHiddenParameters( hiddenParameter, 0 );
		norm = Double.NaN;
		init( freeParams );
		computeLogGammaSum();
	}

	/**
	 * This method is used to append further information of the instance to the
	 * XML representation. This method is designed to allow subclasses to add
	 * information to the XML representation.
	 * 
	 * @return the further information as XML code in a {@link StringBuffer}
	 * 
	 * @see AbstractMixtureScoringFunction#extractFurtherInformation(StringBuffer)
	 */
	protected StringBuffer getFurtherInformation() {
		return new StringBuffer( 1 );
	}

	/**
	 * This method is the opposite of {@link #getFurtherInformation()}. It
	 * extracts further information of the instance from a XML representation.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} containing the information to be
	 *            extracted as XML code
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 * 
	 * @see AbstractMixtureScoringFunction#getFurtherInformation()
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {}

	/**
	 * This array is used to compute the relative indices of a parameter index.
	 * 
	 * @param index
	 *            the parameter index
	 * 
	 * @return the relative indices of the parameter index
	 * 
	 * @see AbstractMixtureScoringFunction#paramRef
	 */
	protected int[] getIndices( int index ) {
		int[] erg = { 0, -1 };
		while( index >= paramRef[erg[0]] ) {
			erg[0]++;
		}
		erg[0]--;
		erg[1] = index - paramRef[erg[0]];
		return erg;
	}

	/**
	 * This method returns the XML tag of the instance that is used to build a
	 * XML representation.
	 * 
	 * @return the XML tag of the instance
	 * 
	 * @see Class#getSimpleName()
	 */
	protected String getXMLTag() {
		return getClass().getSimpleName();
	}

	/**
	 * This method creates the underlying structure for the parameters.
	 * 
	 * @param freeParams
	 *            indicates whether to use only free parameters or all
	 *            parameters
	 */
	protected void init( boolean freeParams ) {
		initWithLength( freeParams, function.length + 2 );
	}

	/**
	 * This method is used to create the underlying structure, e.g.
	 * {@link #paramRef}.
	 * 
	 * @param freeParams
	 *            indicates whether to use free parameters or all
	 * @param len
	 *            the length of the array {@link #paramRef}
	 */
	protected final void initWithLength( boolean freeParams, int len ) {
		if( paramRef == null || paramRef.length != len ) {
			paramRef = new int[len];
		}
		int h, i = 0;
		if( iList == null ) {
			iList = new IntList[Math.max( function.length, hiddenParameter.length )];
			dList = new DoubleList[iList.length];
			for( ; i < iList.length; i++ ) {
				iList[i] = new IntList();
				dList[i] = new DoubleList();
			}
		}
		for( i = 0; i < function.length; i++ ) {
			h = function[i] == null ? 0 : function[i].getNumberOfParameters();
			if( h != UNKNOWN ) {
				paramRef[i + 1] = paramRef[i] + h;
			} else {
				paramRef = null;
				return;
			}
		}
		if( optimizeHidden ) {
			paramRef[i + 1] = paramRef[i] + hiddenParameter.length - ( freeParams ? 1 : 0 );
		} else {
			paramRef[i + 1] = paramRef[i];
		}
		this.freeParams = freeParams;
	}

	/**
	 * This method has to be invoked during an initialization.
	 * 
	 * @param statistic
	 *            a statistic for the initialization of the hidden parameters
	 * @param add
	 *            a switch for adding hyperparameters to the statistic
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 *      boolean, Sample[], double[][])
	 */
	protected void computeHiddenParameter( double[] statistic, boolean add ) {
		int i, j = hiddenParameter.length;
		double sum = 0, offset = 0;
		for( i = 0; i < j; i++ ) {
			if( add ) {
				statistic[i] += getHyperparameterForHiddenParameter( i );
			}
			sum += statistic[i];
		}
		sum = Math.log( sum );
		if( freeParams ) {
			j--;
			offset = Math.log( statistic[j] ) - sum - getLogNormalizationConstantForComponent( j );
			hiddenParameter[j] = 0;
		}
		for( i = 0; i < j; i++ ) {
			hiddenParameter[i] = Math.log( statistic[i] ) - sum -getLogNormalizationConstantForComponent( i ) - offset;
		}
		setHiddenParameters( hiddenParameter, 0 );
	}

	/**
	 * Pre-computes the normalization constant.
	 */
	protected void precomputeNorm() {
		for( int i = 0; i < hiddenPotential.length; i++ ) {
			partNorm[i] = logHiddenPotential[i] + getLogNormalizationConstantForComponent( i );
		}
		norm = Normalisation.getLogSum( partNorm );
	}

	/**
	 * Computes the logarithm of the normalization constant for the component <code>i</code>.
	 * 
	 * @param i
	 *            the index of the component
	 * 
	 * @return the logarithm of the normalization constant of the component
	 */
	protected abstract double getLogNormalizationConstantForComponent( int i );

	/**
	 * Fills the internal array {@link #componentScore} with the logarithmic
	 * scores of the components given a {@link Sequence}.
	 * 
	 * @param seq
	 *            the sequence
	 * @param start
	 *            the start position in <code>seq</code>
	 */
	protected abstract void fillComponentScores( Sequence seq, int start );

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#isNormalized
	 * ()
	 */
	@Override
	public final boolean isNormalized() {
		return isNormalized;
	}

	/**
	 * This method returns a specific internal function.
	 * 
	 * @param index
	 *            the index of the specific function
	 * 
	 * @return a clone of the specific function
	 * 
	 * @throws CloneNotSupportedException
	 *             if the function could not be cloned
	 */
	public NormalizableScoringFunction getFunction( int index ) throws CloneNotSupportedException {
		return (function[index]!= null ? (NormalizableScoringFunction)function[index].clone() : null);
	}

	/**
	 * This method returns an array of clones of the internal used functions.
	 * 
	 * @return an array of clones of the internal used functions
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one function could not be cloned
	 */
	public NormalizableScoringFunction[] getFunctions() throws CloneNotSupportedException {
		return ArrayHandler.clone( function );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.SamplingScoringFunction#getSamplingGroups(int)
	 */
	@Override
	public int[][] getSamplingGroups( int parameterOffset ) {
		LinkedList<int[]> list = new LinkedList<int[]>();
		int i, j, o = 0;
		for( i = 0; i < function.length; i++, o++ ){
			if( function[i] instanceof SamplingScoringFunction ) {
				int[][] current = ((SamplingScoringFunction)function[i]).getSamplingGroups( parameterOffset+paramRef[o] );
				for( j = 0; j < current.length; j++ ) {
					list.add( current[j] );
				}
			} else {
				int[] current = new int[function[i].getNumberOfParameters()];
				for( j = 0; j < current.length; j++ ) {
					current[j] = parameterOffset + paramRef[i]+j;
				}
				list.add( current );
			}
		}
		int[] current = new int[hiddenParameter.length];
		for( j = 0; j < current.length; j++ ) {
			current[j] = parameterOffset + paramRef[function.length]+j;
		}
		list.add( current );
		return list.toArray( new int[0][] );
	}
}
