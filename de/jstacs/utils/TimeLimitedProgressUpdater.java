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
 * This class is an extension of {@link DefaultProgressUpdater}. It prints the
 * percentage of iterations that is already done on the screen. Furthermore it
 * limits the time of the progress. If a progress calls the method
 * {@link #isCancelled()} and the time given in the constructor has passed
 * <code>false</code> is returned. Be aware: this class does not stop a progress
 * after a certain time, it only answers on request. The user has to stop the
 * progress!
 * 
 * @author Jens Keilwagen
 */
public class TimeLimitedProgressUpdater extends DefaultProgressUpdater {

	private long seconds;

	private boolean start;

	private Time t;

	/**
	 * Creates a new {@link TimeLimitedProgressUpdater}.
	 * 
	 * @param t
	 *            the current used time
	 * @param sec
	 *            seconds
	 * @param min
	 *            minutes
	 * @param hours
	 *            hours
	 * @param days
	 *            days
	 * 
	 * @throws IllegalArgumentException
	 *             if one of <code>sec</code>, <code>min</code>,
	 *             <code>hours</code> or <code>days</code> is less than 0
	 * 
	 * @see de.jstacs.utils.Time
	 * @see de.jstacs.utils.RealTime
	 * @see de.jstacs.utils.UserTime
	 */
	public TimeLimitedProgressUpdater( Time t, int sec, int min, int hours, int days ) throws IllegalArgumentException {
		super();
		if( sec < 0 || min < 0 || hours < 0 || days < 0 ) {
			throw new IllegalArgumentException( "The values of sec, min, hours and days have to be non-negative." );
		}
		this.seconds = sec + 60 * min + 3600 * hours + 86400 * days;
		this.t = t;
		start = false;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.DefaultProgressUpdater#setValue(int)
	 */
	@Override
	public void setValue( int value ) {
		if( start == false ) {
			t.reset();
			start = true;
		}
		super.setValue( value );
		if( value == max ) {
			start = false;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.DefaultProgressUpdater#isCancelled()
	 */
	@Override
	public boolean isCancelled() {
		return( t.getElapsedTime() >= seconds );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String erg;
		long e = seconds;
		if( e % 60 < 10 ) erg = ":0" + e % 60;
		else erg = ":" + e % 60;
		e /= 60;
		if( e % 60 < 10 ) erg = ":0" + e % 60 + erg;
		else erg = ":" + e % 60 + erg;
		e /= 60;
		erg = e % 24 + erg;
		if( e % 24 < 10 ) erg = "0" + erg;
		e /= 24;
		return "Timer: " + e + " days and " + erg;
	}
}