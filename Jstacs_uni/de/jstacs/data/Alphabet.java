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

package de.jstacs.data;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.parameters.InstanceParameterSet;

/**
 * Class for a set of symbols, i.e. an {@link Alphabet}. All {@link Alphabet}s
 * are immutable.
 * 
 * @author Jens Keilwagen
 */
public abstract class Alphabet implements Storable, InstantiableFromParameterSet, Comparable<Alphabet> {

	/**
	 * Checks if this {@link Alphabet} is consistent consistent with another
	 * {@link Alphabet}, i.e. both {@link Alphabet}s use the same symbols.
	 * 
	 * @param a
	 *            the second {@link Alphabet}
	 * 
	 * @return <code>true</code> if the {@link Alphabet}s are consistent,
	 *         <code>false</code> otherwise
	 * 
	 * @see Comparable#compareTo(Object)
	 */
	public final boolean checkConsistency( Alphabet a ) {
		return compareTo( a ) == 0;
	}

	/**
	 * Returns the minimal value of the{@link Alphabet}.
	 * 
	 * @return the minimal value of the {@link Alphabet}
	 */
	public abstract double getMin();

	/**
	 * Returns the length of the {@link Alphabet}.
	 * 
	 * @return the length of the {@link Alphabet}
	 */
	public abstract double length();

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public abstract String toString();

	/**
	 * The super class for the {@link InstanceParameterSet} of any
	 * {@link Alphabet}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static abstract class AlphabetParameterSet extends InstanceParameterSet {

		/**
		 * Creates a new {@link AlphabetParameterSet} from the class that can be
		 * instantiated using this {@link AlphabetParameterSet}.
		 * 
		 * @param instanceClass
		 *            the class
		 * 
		 * @see InstanceParameterSet#InstanceParameterSet(Class)
		 */
		public AlphabetParameterSet( Class<? extends Alphabet> instanceClass ) {
			super( instanceClass );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link AlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link AlphabetParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see InstanceParameterSet#InstanceParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public AlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ParameterSet#clone()
		 */
		@Override
		public AlphabetParameterSet clone() throws CloneNotSupportedException {
			return (AlphabetParameterSet)super.clone();
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return getInstanceClass().getSimpleName();
		}
	}
}