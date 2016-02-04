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
// from http://stackoverflow.com/questions/8506881/nice-label-algorithm-for-charts-with-minimum-ticks

/**
 * Class for creating nive tick marks on axes given a minimum and maximum value.
 * Greatly inspired by http://stackoverflow.com/questions/8506881/nice-label-algorithm-for-charts-with-minimum-ticks
 * 
 * @author Jan Grau
 *
 */
public class NiceScale {

	private double minPoint;
	private double maxPoint;
	private double maxTicks = 10;
	private double tickSpacing;



	private double range;
	private double niceMin;



	private double niceMax;

	/**
	 * Creates a {@link NiceScale} object for the given minimum and maximum
	 *
	 * @param min the minimum data point on the axis
	 * @param max the maximum data point on the axis
	 */
	public NiceScale(double min, double max) {
		this.minPoint = min;
		this.maxPoint = max;
		calculate();
	}

	/**
	 * Returns the spacing between the tick marks
	 * @return the spacing
	 */
	public double getTickSpacing() {
		return tickSpacing;
	}

	/**
	 * Returns the "nice" minimum value
	 * @return the minimum
	 */
	public double getNiceMin() {
		return niceMin;
	}

	/**
	 * Returns the "nice" maximum value
	 * @return the maximum
	 */
	public double getNiceMax() {
		return niceMax;
	}

	/**
	 * Calculate and update values for tick spacing and nice
	 * minimum and maximum data points on the axis.
	 */
	private void calculate() {
		this.range = niceNum(maxPoint - minPoint, true);
		this.tickSpacing = niceNum(range / (maxTicks - 1), true);
		this.niceMin =
				Math.ceil(minPoint / tickSpacing) * tickSpacing;
		this.niceMax =
				Math.ceil(maxPoint / tickSpacing) * tickSpacing;
	}

	/**
	 * Returns a "nice" number approximately equal to range. Rounds
	 * the number if round = true, Takes the ceiling if round = false.
	 *
	 * @param range the data range
	 * @param round whether to round the result
	 * @return a "nice" number to be used for the data range
	 */
	private double niceNum(double range, boolean round) {
		double exponent; /** exponent of range */
		double fraction; /** fractional part of range */
		double niceFraction; /** nice, rounded fraction */

		exponent = Math.floor(Math.log10(range));
		fraction = range / Math.pow(10, exponent);

		if (round) {
			if (fraction < 1.5)
				niceFraction = 1;
			else if (fraction < 3)
				niceFraction = 2;
			else if (fraction < 7)
				niceFraction = 5;
			else
				niceFraction = 10;
		} else {
			if (fraction <= 1)
				niceFraction = 1;
			else if (fraction <= 2)
				niceFraction = 2;
			else if (fraction <= 5)
				niceFraction = 5;
			else
				niceFraction = 10;
		}

		return niceFraction * Math.pow(10, exponent);
	}

	/**
	 * Sets the minimum and maximum data points for the axis.
	 *
	 * @param minPoint the minimum data point on the axis
	 * @param maxPoint the maximum data point on the axis
	 */
	public void setMinMaxPoints(double minPoint, double maxPoint) {
		this.minPoint = minPoint;
		this.maxPoint = maxPoint;
		calculate();
	}

	/**
	 * Sets maximum number of tick marks we're comfortable with
	 *
	 * @param maxTicks the maximum number of tick marks for the axis
	 */
	public void setMaxTicks(double maxTicks) {
		this.maxTicks = maxTicks;
		calculate();
	}
}