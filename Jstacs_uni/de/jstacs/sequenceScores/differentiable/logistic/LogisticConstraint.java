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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.sequenceScores.differentiable.logistic;

import de.jstacs.Storable;
import de.jstacs.data.Sequence;

/**
 * This interface defines the function {@latex.inline $f(\\underline{x})$} for sequence {@latex.inline $\\underline{x}$} which can be used in {@link LogisticDiffSS}.
 * 
 * @author Jens Keilwagen
 * 
 * @see LogisticDiffSS
 */
public interface LogisticConstraint extends Storable, Cloneable {

	/**
	 * This method returns the value f(seq) used in {@link LogisticDiffSS}
	 * 
	 * @param seq the sequence
	 * @param start the start position within the sequence
	 * 
	 * @return f(seq)
	 */
	public double getValue( Sequence seq, int start );
}
