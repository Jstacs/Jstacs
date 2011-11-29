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

package de.jstacs.sampling;


import de.jstacs.data.Sample;

/**
 * This is the interface that any {@link de.jstacs.models.AbstractModel} has to implement if it
 * should be used in a sampling.
 * 
 * <br>
 * <br>
 * 
 * During the sampling the method {@link #drawParameters(Sample, double[] )} is used for
 * drawing the parameters from the posterior.
 * 
 * @author Berit Haldemann, Jens Keilwagen
 * 
 * @see de.jstacs.models.AbstractModel
 * @see de.jstacs.models.mixture.AbstractMixtureModel
 * @see de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING
 */
public interface GibbsSamplingModel extends SamplingComponent {

	/**
	 * This method draws the parameters of the model from the a posteriori
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
	 * @param data
	 *            a sample
	 * @param weights
	 *            the (non-negative) weights for each sequence of the sample
	 * 
	 * @throws Exception
	 *             if there is a problem with drawing the parameters, the model
	 *             is not initialized, ...
	 * 
	 * @see SamplingComponent#initForSampling(int)
	 * @see SamplingComponent#acceptParameters()
	 */
	public void drawParameters( Sample data, double[] weights ) throws Exception;
}
