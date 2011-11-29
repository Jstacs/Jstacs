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

package de.jstacs.utils;

/**
 * This class implements a {@link ProgressUpdater} doing nothing but forces a
 * crossvalidation that is used with an instance of this class to continue to
 * its end. <br>
 * May be used, if no progress update is wished during a crossvalidation.
 * 
 * @author Andre Gohr
 */
public class NullProgressUpdater implements ProgressUpdater {

	private static NullProgressUpdater NULLPROGRESSUPDATER = new NullProgressUpdater();

	private NullProgressUpdater() {}

	/**
	 * Returns a reference to the same {@link NullProgressUpdater} that is
	 * immutable.
	 * 
	 * @return always a reference to the same {@link NullProgressUpdater} that
	 *         is immutable since all methods except {@link #isCancelled()} are
	 *         not implemented. The method {@link #isCancelled()} always returns
	 *         <code>false</code> to force the crossvalidation to continue
	 */
	public static NullProgressUpdater getImmutableInstance() {
		return NullProgressUpdater.NULLPROGRESSUPDATER;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setMax(int)
	 */
	public void setMax( int max ) {}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setValue(int)
	 */
	public void setValue( int value ) {}

	/**
	 * Closes the supervision.
	 */
	public void close() {}

	/**
	 * After {@link NullProgressUpdater#setOffset()} is called the current value
	 * will be added to every value set by
	 * {@link NullProgressUpdater#setValue(int)}.
	 */
	public void setOffset() {}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#isCancelled()
	 */
	public boolean isCancelled() {
		return false;
	}

}
