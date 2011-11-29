package de.jstacs.models.hmm.states.emissions.discrete;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DiscreteAlphabet;

/**
 * This class implements a simple discrete emission without any condition.
 * 
 * @author Jens Keilwagen, Michael Scharfe, Jan Grau
 */
public class DiscreteEmission extends AbstractConditionalDiscreteEmission {
	
	/**
	 * This is a simple constructor for a {@link DiscreteEmission} based on the equivalent sample size.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
	 * 
	 * @see #DiscreteEmission(AlphabetContainer, double[])
	 */
	public DiscreteEmission( AlphabetContainer con, double ess ) {
		super( con, 1, ess );
	}

	/**
	 * This is a simple constructor for a {@link DiscreteEmission} defining the individual hyper parameters.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param hyperParams the individual hyper parameters for each parameter
	 */
	public DiscreteEmission( AlphabetContainer con, double[] hyperParams ) {
		super( con, new double[][]{hyperParams} );
	}

	/**
	 * Creates a {@link DiscreteEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public DiscreteEmission( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected int getConditionIndex( boolean forward, int seqPos, Sequence seq ) {
		return 0;
	}

	@Override
	public String toString() {
		String res = "";
		DiscreteAlphabet abc = (DiscreteAlphabet) con.getAlphabetAt( 0 );
		for( int i = 0; i < probs[0].length; i++ ) {
			res += "P(X=" + abc.getSymbolAt( i ) + ") = " + probs[0][i] + "\t";
		}
		return res+"\n";
	}
}