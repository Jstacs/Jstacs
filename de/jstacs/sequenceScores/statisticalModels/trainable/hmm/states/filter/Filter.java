package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;

/**
 * A filter based on a {@link Sequence} and a position.
 * This simple filter might be applied in/for {@link State}s of an {@link AbstractHMM}.
 * 
 * @author Jens Keilwagen
 */
public interface Filter extends Storable {
	/**
	 * Returns <code>true</code> if the filter is passed.
	 * 
	 * @param anchor an anchor position
	 * @param seq the sequence
	 * 
	 * @return <code>true</code> if the filter is passed
	 */
	public boolean isAccepted( int anchor, Sequence seq );
}