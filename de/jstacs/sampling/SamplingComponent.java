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

import java.io.IOException;

/**
 * This interface defines methods that are used during a sampling.
 * 
 * Before starting a series of samplings the method
 * {@link SamplingComponent#initForSampling(int)} is called. <br>
 * 
 * Then before starting a specific sampling the method
 * {@link SamplingComponent#extendSampling(int, boolean)} is called. <br>
 * 
 * During the sampling the parameters have to be sampled from the posterior
 * by a specialized method that is defined in one of the sub-interfaces.
 * The sampled parameters can than be written to a file using the method
 * {@link SamplingComponent#acceptParameters()}. <br>
 * 
 * After finishing the sampling the method
 * {@link SamplingComponent#samplingStopped()} is called. <br>
 * 
 * After a sampling the methods
 * {@link SamplingComponent#parseParameterSet(int, int)} and
 * {@link SamplingComponent#parseNextParameterSet()} will be used for
 * computing for instance the (log-) likelihoods. <br>
 * 
 * <br>
 * <br>
 * 
 * The method {@link SamplingComponent#isInSamplingMode()} can be used to
 * check whether an object is currently used for sampling. If the object is in
 * sampling mode, it should not support any other method for changing the
 * parameters than the specialized sampling method of the sub-interfaces and
 * {@link SamplingComponent#parseParameterSet(int, int)}. Furthermore it is
 * legal to throw an {@link Exception} when the object is in sampling mode and a
 * method for saving, cloning or training the parameters is called.
 * 
 * 
 * @author Berit Haldemann, Jens Keilwagen, Michael Scharfe
 */
public interface SamplingComponent {

	/**
	 * This method allows the user to parse the set of parameters with index
	 * <code>n</code> of a certain <code>sampling</code> (from a file). The
	 * internal numbering should start with 0. The parameter set with index 0 is
	 * the initial (random) parameter set. It is recommended that a series of
	 * parameter sets is accessed by the following lines:
	 * 
	 * <p>
	 * <code>
	 * for( sampling = 0; sampling < numSampling; sampling++ )<br>
	 * {
	 * <dir>
	 *     boolean b = parseParameterSet( sampling, n );<br>
	 *     while( b )<br>
	 *     {<br>
	 *        //do something<br>
	 *        b = parseNextParameterSet();<br>
	 *     }
	 * </dir>
	 * }<br>
	 * </code>
	 * </p>
	 * 
	 * @param sampling
	 *            the index of the sampling
	 * @param n
	 *            the index of the parameter set
	 * 
	 * @return <code>true</code> if the parameter set could be parsed
	 * 
	 * @throws Exception
	 *             if there is a problem with parsing the parameters
	 * 
	 * @see SamplingComponent#parseNextParameterSet()
	 */
	public abstract boolean parseParameterSet( int sampling, int n ) throws Exception;

	/**
	 * This method allows the user to parse the next set of parameters (from a
	 * file).
	 * 
	 * @return <code>true</code> if the parameters could be parsed, otherwise
	 *         <code>false</code>
	 * 
	 * @see SamplingComponent#parseParameterSet(int, int)
	 */
	public abstract boolean parseNextParameterSet();

	/**
	 * This method initializes the instance for the sampling. For instance this
	 * method can be used to create new files where all parameter sets are
	 * stored.
	 * 
	 * @param starts
	 *            the number of different sampling starts that will be done
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see java.io.File#createTempFile(String, String, java.io.File )
	 */
	public abstract void initForSampling( int starts ) throws IOException;

	/**
	 * This method allows to extend a sampling.
	 * 
	 * @param sampling
	 *            the index of the sampling
	 * @param append
	 *            whether to append the sampled parameters to an existing file
	 *            or to overwrite the file
	 * 
	 * @throws IOException
	 *             if the file could not be handled correctly
	 */
	public abstract void extendSampling( int sampling, boolean append ) throws IOException;

	/**
	 * This method is the opposite of the method
	 * {@link SamplingComponent#extendSampling(int, boolean)}. It can be
	 * used for closing any streams of writer, ...
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see SamplingComponent#extendSampling(int, boolean)
	 */
	public abstract void samplingStopped() throws IOException;

	/**
	 * This method returns <code>true</code> if the object is currently used in
	 * a sampling, otherwise <code>false</code>.
	 * 
	 * @return <code>true</code> if the object is currently used in a sampling,
	 *         otherwise <code>false</code>
	 */
	public abstract boolean isInSamplingMode();

    /**
	 * This methods accepts the drawn parameters.
	 * Internally the drawn parameters should be saved (to a file).
	 * 
	 * @throws IOException
	 *             if the file could not be handled correctly
	 */
	public void acceptParameters() throws IOException;
}