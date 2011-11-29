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
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sampling;


/**
 * This is the interface for sampling based on a sufficient statistic.
 * 
 * <br>
 * <br>
 * 
 * During the sampling the method {@link #drawParametersFromStatistic()} is used for
 * drawing the parameters from the posterior.
 * 
 * @author Jens Keilwagen, Michael Scharfe
 */
public interface SamplingFromStatistic extends SamplingComponent {

	/**
	 * This method draws the parameters using a sufficient statistic representing a posteriori
	 * density. It is recommended to write the parameters to a specific file using
	 * {@link SamplingComponent#acceptParameters()} so that they can later be parsed using the
	 * methods of the interface.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Before using this method the method {@link SamplingComponent#initForSampling(int)} should be
	 * called.
	 * 
	 * @throws Exception
	 *             if there is a problem with drawing the parameters, the initialization, ...
	 * 
	 * @see SamplingComponent#initForSampling(int)
	 * @see SamplingComponent#acceptParameters()
	 */
	public void drawParametersFromStatistic() throws Exception;
	
	/**
	* This method calculates the a-posteriori probability for the current statistics
	*
	* @return the logarithm of the a-posteriori probability
	*/ 
	public double getLogPosteriorFromStatistic();
}
