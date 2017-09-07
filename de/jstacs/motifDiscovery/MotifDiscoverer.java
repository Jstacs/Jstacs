/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.motifDiscovery;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;

/**
 * This is the interface that any tool for de-novo motif discovery should
 * implement.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface MotifDiscoverer extends Cloneable, Storable {

	/**
	 * This <code>enum</code> can be used to determine which kind of profile
	 * should be returned.
	 * 
	 * @see MotifDiscoverer#getProfileOfScoresFor(int, int, Sequence, int,
	 *      KindOfProfile)
	 */
	public enum KindOfProfile {
		/**
		 * The profile should return
		 * <code>\log Q(seq,pos,component|class)</code>.
		 */
		UNNORMALIZED_JOINT,
		/**
		 * The profile should return
		 * <code>\log P(seq,pos|component,class)</code>.
		 */
		NORMALIZED_CONDITIONAL,
		/**
		 * The profile should return
		 * <code>\log Q(seq,pos|component,class)</code>.
		 */
		UNNORMALIZED_CONDITIONAL
	}
	
	/**
	 * This method returns a deep clone of the instance.
	 * 
	 * @return a deep clone of the instance
	 * 
	 * @throws CloneNotSupportedException if the instance could not be cloned
	 * 
	 * @see Cloneable
	 */
	public MotifDiscoverer clone() throws CloneNotSupportedException;

	/**
	 * This method returns the length of the motif with index <code>motif</code>
	 * .
	 * 
	 * @param motif
	 *            the index of the motif
	 * 
	 * @return the length of the motif with index <code>motif</code>
	 */
	public int getMotifLength(int motif);

	/**
	 * Returns the number of components in this {@link MotifDiscoverer}.
	 * 
	 * @return the number of components
	 */
	public int getNumberOfComponents();

	/**
	 * Returns the number of motifs for this {@link MotifDiscoverer}.
	 * 
	 * @return the number of motifs
	 */
	public int getNumberOfMotifs();

	/**
	 * Returns the number of motifs that are used in the component
	 * <code>component</code> of this {@link MotifDiscoverer}.
	 * 
	 * @param component
	 *            the component of the {@link MotifDiscoverer}
	 * 
	 * @return the number of motifs
	 */
	public int getNumberOfMotifsInComponent(int component);

	/**
	 * Returns the index of the component with the maximal score for the
	 * sequence <code>sequence</code>.
	 * 
	 * @param sequence
	 *            the given sequence
	 * 
	 * @return the index of the component with the maximal score for the given
	 *         sequence
	 * 
	 * @throws Exception
	 *             if the index could not be computed for any reasons
	 */
	public int getIndexOfMaximalComponentFor(Sequence sequence)
			throws Exception;

	/**
	 * Returns the global index of the <code>motif</code> used in
	 * <code>component</code>. The index returned must be at least 0 and less
	 * than {@link #getNumberOfMotifs()}.
	 * 
	 * @param component
	 *            the component index
	 * @param motif
	 *            the motif index in the component
	 * 
	 * @return the global index of the
	 *         <code>motif</code> in <code>component</code>
	 */
	public int getGlobalIndexOfMotifInComponent(int component, int motif);

	/**
	 * Returns the profile of the scores for component <code>component</code>
	 * and motif <code>motif</code> at all possible start positions of the motif
	 * in the sequence <code>sequence</code> beginning at <code>startpos</code>.
	 * This array should be of length <br>
	 * <code>sequence.length() - startpos - motifs[motif].getLength() + 1</code>.
	 * 
	 * <br>
	 * 
	 * A high score should encode for a probable start position.
	 * 
	 * @param component
	 *            the component index
	 * @param motif
	 *            the index of the motif in the component
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position in the sequence
	 * @param kind
	 *            indicates the kind of profile
	 * 
	 * @return the profile of scores
	 * 
	 * @throws Exception
	 *             if the score could not be computed for any reasons
	 */
	public double[] getProfileOfScoresFor(int component, int motif,
			Sequence sequence, int startpos, KindOfProfile kind)
			throws Exception;

	/**
	 * This method returns the probabilities of the strand orientations for a given subsequence if it is
	 * considered as site of the motif model in a specific component.
	 * 
	 * @param component
	 *            the component index
	 * @param motif
	 *            the index of the motif in the component
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position in the sequence
	 * 
	 * @return the probabilities of the strand orientations
	 * 
	 * @throws Exception
	 *             if the strand could not be computed for any reasons
	 */
	public double[] getStrandProbabilitiesFor(int component, int motif, Sequence sequence,
			int startpos) throws Exception;
}
