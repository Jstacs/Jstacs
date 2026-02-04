package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;

/**
 * A higher order discrete emission for a HMM.
 * 
 * @author Jens Keilwagen
 * 
 * @see AbstractHMM
 */
public class HigherOrderDiscreteEmission extends AbstractConditionalDiscreteEmission {

	/**
	 * The powers of the alphabet size.
	 */
	protected int[] power; 
	
	/**
	 * The Markov order.
	 */
	protected int order;
	
	/**
	 * Creates a higher order discrete {@link Emission} with given order.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param ess the ess that is equally distributed to all parameters
	 * @param order the non negative Markov order
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public HigherOrderDiscreteEmission(AlphabetContainer con, double ess, int order ) throws WrongAlphabetException {
		this( con, ess, ess, order );
	}
	
	/**
	 * Creates a higher order discrete {@link Emission} with given order.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param ess the ess that is equally distributed to all parameters
	 * @param initESS the ess for the initialization
	 * @param order the non negative Markov order
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public HigherOrderDiscreteEmission(AlphabetContainer con, double ess, double initESS, int order ) throws WrongAlphabetException {
		this(con,ess,initESS,order,true);
	}
	
	/**
	 * Creates a higher order discrete {@link Emission} with given order.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param ess the ess that is equally distributed to all parameters
	 * @param initESS the ess for the initialization
	 * @param order the non negative Markov order
	 * @param norm whether a normalized or unnormalized variant should be created
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public HigherOrderDiscreteEmission(AlphabetContainer con, double ess, double initESS, int order, boolean norm ) throws WrongAlphabetException {
		this(con,ess,initESS,order, null,norm);
	}

	/**
	 * Creates a higher order discrete {@link Emission} with given order.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param ess the ess that is equally distributed to all parameters
	 * @param initESS the ess for the initialization
	 * @param order the non negative Markov order
	 * @param frac the fraction of the <code>ess</code> that is used for each parameter as hyper-parameter
	 * @param norm whether a normalized or unnormalized variant should be created
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public HigherOrderDiscreteEmission(AlphabetContainer con, double ess, double initESS, int order, double[] frac, boolean norm ) throws WrongAlphabetException {
		super(con, (int) Math.pow(con.getAlphabetLengthAt(0), order), ess, initESS, frac, norm);
		if( order < 0 ) {
			throw new IllegalArgumentException("The order must be non-negative");
		}	
		this.order = order;
		createPower();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link HigherOrderDiscreteEmission} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HigherOrderDiscreteEmission} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public HigherOrderDiscreteEmission(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	public void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		order = (Integer) XMLParser.extractObjectForTags(xml, "order");
		createPower();
	}
	
	/**
	 * Creates and fills {@link #power}.
	 */
	protected void createPower() {
		power = new int[order];
		if( order > 0 ) {
			power[0]=1;
			int a = (int) con.getAlphabetLengthAt(0);
			for( int i = 1; i < order; i++ ) {
				power[i] = power[i-1]*a;
			}
		}
	}
	
	protected void appendFurtherInformation( StringBuffer xml ) {
		XMLParser.appendObjectWithTags(xml, order, "order");
	}
	
	public HigherOrderDiscreteEmission clone() throws CloneNotSupportedException {
		HigherOrderDiscreteEmission clone = (HigherOrderDiscreteEmission) super.clone();
		clone.power = power.clone();
		return clone;
	}

	@Override
	protected int getConditionIndex(boolean forward, int seqPos, Sequence seq) {
		int res = 0;
		for( int i = 0; i < order; i++ ) {
			int p = seqPos-(i+1);
			if( p < 0 ) return -1;
			res += power[i]*seq.discreteVal(p);
		}
		return res;
	}

	@Override
	protected String getCondition(int i) {
		StringBuffer sb = new StringBuffer();
		sb.append(" ");
		int a = (int) con.getAlphabetLengthAt(0);
		for( int j = 0; j < order; j++ ) {
			int x = i % a;
			sb.append( (j==0?"|":",")+" X_{-"+(j+1)+"}=" + con.getSymbol(0, x) );
			i /= a;
		}
		return sb.toString();
	}
}